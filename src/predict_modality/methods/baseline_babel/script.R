cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(DropletUtils)
# #https://rdrr.io/github/barkasn/nbHelpers/man/transpose_dgRMatrix.html
library(nbHelpers)

## VIASH START
par <- list(
  input_mod1 = "output_censored_rna.h5ad",
  input_mod2 = "output_censored_mod2.h5ad",
  output = "output_babel.h5ad",
  data_id = "granulocyte_multiome10x",
  mod1 = 'RNA'
)
## VIASH END

# folder names for Babel
babel_input_folder = 'export_train/' # location of temporary input files
babel_location = '../babel/bin/'     # location of babel executables
babel_model_output = 'babel_model'   # babel trained model location
babel_pred_output = 'babel_output'

cat("Reading h5ad files\n")
ad1 <- anndata::read_h5ad(par$input_mod1)
ad2 <- anndata::read_h5ad(par$input_mod2)

# specify feature type for each adata
if(par$mod1=='RNA'){
  ad1$var$feature_types = 'Gene Expression'
  ad2$var$feature_types = 'Peaks'
}else{
  ad1$var$feature_types = 'Peaks'
  ad2$var$feature_types = 'Gene Expression'
}

# subset train and merge both categories (required by Babel)
ad1_train <- ad1[ad1$obs$split == 'train',]

# multiome_matrix for export to Babel's input format
multiome_matrix = cbind(ad1_train$X, ad2$X)

# generate multiome anndata object
ad_babel <- anndata::AnnData(
  X = multiome_matrix,
  var = rbind(ad1_train$var, ad2$var),
  obs = ad1_train$obs
)

# write adata objects as 10x-CellRanger H5 format
# train object: test + train splits
write10xCounts(paste(babel_input_folder,'train_input_babel.h5', sep=""), t(ad_babel$X),
							gene.id = row.names(ad_babel$var), gene.symbol=row.names(ad_babel$var),
							barcodes = row.names(ad_babel$obs), type='HDF5',
							version='3', genome='GRCh38',
							gene.type=ad_babel$var$feature_types)

# test objet:
ad1_test = ad1[ad1$obs$split == 'test',]
x_test = transpose_dgRMatrix(ad1_test$X)

write10xCounts( paste(babel_input_folder, 'test_input_babel.h5', sep=""), x_test,
							gene.id = row.names(ad1_test$var), gene.symbol=row.names(ad1_test$var),
							barcodes = row.names(ad1_test$obs), type='HDF5',
							version='3', genome='GRCh38',
							gene.type=ad1_test$var$feature_types)


# RUN Babel
# switch conda enviornment
# Babel requires addata v0.6
use_condaenv("babel")

# train
babel_output_dir = paste(babel_model_output,"_",par$data_id, sep="")

babel_train_cmd = paste("python ",babel_location,"train_model.py --data ",babel_input_folder,"train_input_babel.h5 --outdir ",babel_output_dir,sep="")
system(babel_train_cmd)
# test
babel_prediction = paste(babel_pred_output,"_",par$data_id, sep="")
babel_pred_cmd = paste("python ",babel_location,"predict_model.py --checkpoint ",babel_output_dir," --data ",babel_input_folder,"test_input_babel.h5 --outdir ",babel_prediction, sep="")
system(babel_pred_cmd)

# Babel generated a model folder with h5 files for the training set
# After prediction of the test data, the final h5 file is in ./babel_model_output
# format [input:RNA, prediction: ATAC]
out = anndata::read_h5ad(paste(babel_prediction, '/rna_atac_adata.h5ad', sep=""))

# final h5 for prediction of test data
# add meta.data for the output
out$uns = list(method_id = 'babel')

# write output file
zzz <- out$write_h5ad(par$output, compression = "gzip")
