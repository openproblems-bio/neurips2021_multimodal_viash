cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(assertthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input_rna = "resources_test/common/pbmc_1k_protein_v3.normalize.output_rna.h5ad",
  input_mod2 = "resources_test/common/pbmc_1k_protein_v3.normalize.output_mod2.h5ad",
  output_mod1 = "resources_test/task1/pbmc_1k_protein_v3.output_mod1.h5ad",
  output_mod2 = "resources_test/task1/pbmc_1k_protein_v3.output_mod2.h5ad",
  output_solution = "resources_test/task1/pbmc_1k_protein_v3.solution.h5ad",
  mod1 = "RNA",
  seed = 1L
)
## VIASH END

cat("Reading input data\n")
rna_mod1 <- tolower(par$mod1) == "rna"
ad1_path <- if (rna_mod1) { par$input_rna } else { par$input_mod2 }
ad2_path <- if (rna_mod1) { par$input_mod2 } else { par$input_rna }

ad1_raw <- anndata::read_h5ad(ad1_path)
ad2_raw <- anndata::read_h5ad(ad2_path)

cat("Determining train/test split\n")
split <- 
  if (!is.null(ad1_raw$obs[["experiment"]]) && all(ad1_raw$obs[["experiment"]] %in% c("train", "test"))) {
    ad1_raw$obs[["experiment"]]
  } else {
    set.seed(par$seed)
    ix <- sample.int(
      nrow(ad1_raw), 
      size = nrow(ad1_raw) * 0.66,
      replace = FALSE
    )
    ifelse(seq_len(nrow(ad1_raw)) %in% ix, "train", "test")
  }
splor <- order(split)
ad1_raw <- ad1_raw[splor, ]
ad2_raw <- ad2_raw[splor, ]
split <- split[splor]

cat("Creating mod1 object\n")
out_mod1 <- anndata::AnnData(
  X = ad1_raw$X,
  obs = data.frame(
    row.names = rownames(ad1_raw),
    split = split
  ),
  uns = list(
    dataset_id = paste0(ad1_raw$uns[["dataset_id"]], "_task1"),
    modality = ad1_raw$uns[["modality"]]
  )
)

cat("Creating mod2 object\n")
out_mod2 <- anndata::AnnData(
  X = ad2_raw$X[split == "train", ],
  obs = data.frame(
    row.names = rownames(ad2_raw)[split == "train"],
    split = split[split == "train"]
  ),
  uns = list(
    dataset_id = paste0(ad2_raw$uns[["dataset_id"]], "_task1"),
    modality = ad2_raw$uns[["modality"]]
  )
)

cat("Create solution object\n")
out_solution <- anndata::AnnData(
  X = ad2_raw$X[split == "test", ],
  obs = data.frame(
    row.names = rownames(ad2_raw)[split == "test"],
    split = split[split == "test"]
  ),
  uns = list(
    dataset_id = paste0(ad2_raw$uns[["dataset_id"]], "_task1"),
    modality = ad2_raw$uns[["modality"]]
  )
)

cat("Saving output files as h5ad\n")
zzz <- out_mod1$write_h5ad(par$output_mod1, compression = "gzip")
zzz <- out_mod2$write_h5ad(par$output_mod2, compression = "gzip")
zzz <- out_solution$write_h5ad(par$output_solution, compression = "gzip")