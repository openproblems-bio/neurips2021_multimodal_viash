cat(">> Loading dependencies\n")

requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(testthat, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("DropletUtils", quietly = TRUE)

options(tidyverse.quiet = TRUE)
library(tidyverse)

## VIASH START
par <- list(
  input_mod1 = "output/task1_datasets/10x_pbmc_granulocyte_sorted_3k_rna/10x_pbmc_granulocyte_sorted_3k_rna.prepare_task1_dataset.output_mod1.h5ad",
  input_mod2 = "output/task1_datasets/10x_pbmc_granulocyte_sorted_3k_rna/10x_pbmc_granulocyte_sorted_3k_rna.prepare_task1_dataset.output_mod2.h5ad",
  output = "output_babel.h5ad",
  data_id = "granulocyte_multiome10x"
)
if (!dir.exists(dirname(par$input_mod1))) {
  stop(
    "Directory '", dirname(par$input_mod1), "' does not exist.\n",
    "You probably need to run `aws s3 sync s3://neurips2021-multimodal-public-datasets/ output/` first."
  )
}
## VIASH END

babel_location <- "../babel/bin/"      # location of babel executables

cat(">> Reading h5ad files\n")
ad1 <- anndata::read_h5ad(par$input_mod1)
if (is.null(ad1$var$gene_ids)) ad1$var$gene_ids <- colnames(ad1)
ad2 <- anndata::read_h5ad(par$input_mod2)
if (is.null(ad2$var$gene_ids)) ad2$var$gene_ids <- colnames(ad2)

mod1 <- unique(ad1$var$feature_types)

# subset train and test
ad1_train <- ad1[ad1$obs$group == "train", ]
ad1_test <- ad1[ad1$obs$group == "test", ]

# multiome_matrix for export to Babel's input format
multiome_matrix <- cbind(ad1_train$X, ad2$X)

# generate multiome anndata objects
ad_babel <- anndata::AnnData(
  X = multiome_matrix,
  var = bind_rows(ad1_train$var, ad2$var),
  obs = ad1_train$obs
)

# setting up babel dirs
tmpdir <- tempfile(pattern = "babel_temp", fileext = "/")
cat(">> Setting up directories for babel at ", tmpdir, "\n", sep = "")
dir.create(tmpdir)
on.exit(unlink(tmpdir, recursive = TRUE))

dir_data <- paste0(tmpdir, "data/")     # location of input files
dir.create(dir_data)
dir_model <- paste0(tmpdir, "model/")   # location of babel model
dir_pred <- paste0(tmpdir, "pred/")     # location of predictions

feature_type_map <- c(
  "GEX" = "Gene Expression",
  "ADT" = "Peaks", # try to make it run on ADT data as well
  "ATAC" = "Peaks"
)

cat(">> Writing train dataset as 10x-CellRanger H5 format\n")
DropletUtils::write10xCounts(
  paste0(dir_data, "train_input.h5"),
  t(ad_babel$X),
  gene.id = ad_babel$var$gene_ids,
  gene.symbol = colnames(ad_babel),
  barcodes = rownames(ad_babel),
  type = "HDF5",
  version = "3",
  genome = "GRCh38",
  gene.type = unname(feature_type_map[ad_babel$var$feature_types]),
  overwrite = TRUE
)

cat(">> Writing test dataset as 10x-CellRanger H5 format\n")
DropletUtils::write10xCounts(
  paste0(dir_data, "test_input.h5"),
  t(ad1_test$X),
  gene.id = ad1_test$var$gene_ids,
  gene.symbol = colnames(ad1_test),
  barcodes = rownames(ad1_test),
  type = "HDF5",
  version = "3",
  genome = "GRCh38",
  gene.type = unname(feature_type_map[ad1_test$var$feature_types]),
  overwrite = TRUE
)

cat(">> Babel: train model\n")
babel_train_cmd <- paste0(
  "/opt/conda/bin/conda run -n babel ",
  "python ", babel_location, "train_model.py ",
  "--data ", dir_data, "train_input.h5 ",
  "--outdir ", dir_model
)
out1 <- system(babel_train_cmd)

# check whether training succeeded
expect_equal(out1, 0, info = paste0("Model training failed with exit code ", out1))

cat(">> Babel: predict from model\n")
babel_pred_cmd <- paste0(
  "/opt/conda/bin/conda run -n babel ",
  "python ", babel_location, "predict_model.py ",
  "--checkpoint ", dir_model, " ",
  "--data ", dir_data, "test_input.h5 ",
  "--outdir ", dir_pred
)
out2 <- system(babel_pred_cmd)

# check whether prediction succeeded
expect_equal(out2, 0, info = paste0("Prediction failed with exit code ", out1))

cat(">> Read predictions\n")
pred <- anndata::read_h5ad(paste0(dir_pred, "/rna_atac_adata.h5ad"))

cat(">> Creating output object\n")
pred_long <-
  summary(pred$X) %>%
  mutate(
    i = match(rownames(pred), rownames(ad1_test))[i],
    j = match(colnames(pred), colnames(ad2))[j]
  )
pred_expanded <- 
  Matrix::sparseMatrix(
    i = pred_long$i,
    j = pred_long$j,
    x = pred_long$x,
    dims = c(nrow(ad1_test), ncol(ad2)),
    dimnames = list(rownames(ad1_test), colnames(ad2))
  ) %>%
  as("CsparseMatrix")

out <- anndata::AnnData(
  X = pred_expanded,
  uns = list(
    method_id = "babel"
  )
)

cat(">> Write predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")
