#!/bin/bash

set -ex
export PATH="${CONDA_BIN}:${PATH}"

cd $(dirname $0)

echo "installing root env:"
cat /tmp/environment.yml
conda install -c conda-forge mamba
mamba env update -n root  -f /tmp/environment.yml
jupyter serverextension enable --sys-prefix jupyter_server_proxy
jupyter serverextension enable dask_labextension --sys-prefix
jupyter serverextension enable --py jupyterlab_code_formatter --sys-prefix

# jupyter labextension install @bokeh/jupyter_bokeh
# jupyter labextension install @jupyter-widgets/jupyterlab-manager
# jupyter labextension install dask-labextension@3.0.0
# jupyter labextension install @ryantam626/jupyterlab_code_formatter
# jupyter labextension install @pyviz/jupyterlab_pyviz
# jupyter labextension install jupyterlab-execute-time
# jupyter labextension install @telamonian/theme-darcula
# jupyter labextension install jupyterlab-python-file

conda clean -afy
jupyter lab clean
jlpm cache clean
npm cache clean --force
find ${CONDA_DIR}/ -type f,l -name '*.pyc' -delete
find ${CONDA_DIR}/ -type f,l -name '*.a' -delete
find ${CONDA_DIR}/ -type f,l -name '*.js.map' -delete
rm -rf $HOME/.node-gyp
rm -rf $HOME/.local

conda create -n saturn
