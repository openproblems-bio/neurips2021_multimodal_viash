cat(">> Loading dependencies\n")

requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(testthat, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("DropletUtils", quietly = TRUE)

options(tidyverse.quiet = TRUE)
library(tidyverse)

babel_location <- "/babel/bin/"
conda_bin <- "/opt/conda/bin/conda"

## VIASH START
# par <- list(
#   input_train_mod1 = "resources_test/predict_modality/test_resource.train_mod1.h5ad",
#   input_train_mod2 = "resources_test/predict_modality/test_resource.train_mod2.h5ad",
#   input_test_mod1 = "resources_test/predict_modality/test_resource.test_mod1.h5ad",
#   output = "output.h5ad"
# )
par <- list(
  input_train_mod1 = "output/public_datasets/predict_modality/totalvi_spleen_lymph_111_rna/totalvi_spleen_lymph_111_rna.censor_dataset.output_train_mod1.h5ad",
  input_train_mod2 = "output/public_datasets/predict_modality/totalvi_spleen_lymph_111_rna/totalvi_spleen_lymph_111_rna.censor_dataset.output_train_mod2.h5ad",
  input_test_mod1 = "output/public_datasets/predict_modality/totalvi_spleen_lymph_111_rna/totalvi_spleen_lymph_111_rna.censor_dataset.output_test_mod1.h5ad",
  output = "output.h5ad"
)
conda_bin <- "conda"
babel_location <- "../babel/bin/"
## VIASH END

cat(">> Reading h5ad files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)
if (is.null(input_train_mod1$var$gene_ids)) input_train_mod1$var$gene_ids <- colnames(input_train_mod1)
if (is.null(input_train_mod2$var$gene_ids)) input_train_mod2$var$gene_ids <- colnames(input_train_mod2)
if (is.null(input_test_mod1$var$gene_ids)) input_test_mod1$var$gene_ids <- colnames(input_test_mod1)

mod1 <- as.character(unique(input_train_mod1$var$feature_types))
mod2 <- as.character(unique(input_train_mod2$var$feature_types))
expect_equal(mod1, "GEX", info = "babel only runs on GEX→ATAC")
expect_equal(mod2, "ATAC", info = "babel only runs on GEX→ATAC")

# multiome_matrix for export to Babel's input format
multiome_matrix <- cbind(input_train_mod1$X, input_train_mod2$X)

# generate multiome anndata object for training
ad_babel <- anndata::AnnData(
  X = multiome_matrix,
  var = bind_rows(input_train_mod1$var, input_train_mod2$var)
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
  t(input_test_mod1$X),
  gene.id = input_test_mod1$var$gene_ids,
  gene.symbol = colnames(input_test_mod1),
  barcodes = rownames(input_test_mod1),
  type = "HDF5",
  version = "3",
  genome = "GRCh38",
  gene.type = unname(feature_type_map[input_test_mod1$var$feature_types]),
  overwrite = TRUE
)


cat(">> Babel: train model\n")
babel_train_cmd <- paste0(
  conda_bin, " run -n babel ",
  "python ", babel_location, "train_model.py ",
  "--data ", dir_data, "train_input.h5 ",
  "--outdir ", dir_model
)
# stringent filtering causes babel to almost always fail
# reason: https://github.com/wukevin/babel/blob/main/babel/sc_data_loaders.py#L168-L190

out1 <- system(babel_train_cmd)

# check whether training succeeded
expect_equal(out1, 0, info = paste0("Model training failed with exit code ", out1))

cat(">> Babel: predict from model\n")
babel_pred_cmd <- paste0(
  conda_bin, " run -n babel ",
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
