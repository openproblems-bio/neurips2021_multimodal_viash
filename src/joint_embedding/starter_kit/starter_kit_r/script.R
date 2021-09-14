# Dependencies:
#   python: anndata
#   r: anndata, lmds
#
# R starter kit for the NeurIPS 2021 Single-Cell Competition.
# Parts with `TODO` are supposed to be changed by you.
#
# More documentation:
#
# https://viash.io/docs/creating_components/r/

cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE, quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(lmds, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
# Anything within this block will be removed by viash
# and will be replaced with the parameters as specified in
# your config.vsh.yaml.

dataset_path <- "sample_data/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter."
# dataset_path <- "output/datasets/joint_embedding/openproblems_bmmc_multiome_phase1/openproblems_bmmc_multiome_phase1.censor_dataset.output_"

par <- list(
  input_mod1 = paste0(dataset_path, "mod1.h5ad"),
  input_mod2 = paste0(dataset_path, "mod2.h5ad"),
  output = "output.h5ad",
  distance_method = "spearman",
  n_pcs = 50L
)
## VIASH END

method_id <- "r_starter_kit" # fill in the name of your method here

cat("Reading mod1 h5ad file\n")
input_mod1 <- read_h5ad(par$input_mod1)

cat("Performing DR on mod1\n")
dr_mod1 <- lmds(input_mod1$X, ndim = par$n_pcs / 2, distance_method = par$distance_method)
input_mod1_obs <- input_mod1$obs
input_mod1_uns <- input_mod1$uns

cat("Clearing mod1 from memory\n")
rm(input_mod1)
gc()

cat("Reading mod2 h5ad file\n")
input_mod2 <- read_h5ad(par$input_mod2)

cat("Performing DR on mod2\n")
dr_mod2 <- lmds(input_mod2$X, ndim = par$n_pcs / 2, distance_method = par$distance_method)

cat("Clearing mod2 from memory\n")
rm(input_mod2)
gc()

cat("Creating output AnnData\n")
dr <- cbind(dr_mod1, dr_mod2)
out <- anndata::AnnData(
  X = dr,
  uns = list(
    dataset_id = input_mod1_uns[["dataset_id"]],
    method_id = method_id
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")