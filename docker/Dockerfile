from dataintuitive/viash:latest

RUN rm -rf /etc/apt/sources.list.d/sbt*
RUN mkdir /opt/conda && \
    ln -s /opt/conda /srv/conda
ENV APP_BASE=/srv
ENV CONDA_DIR=/srv/conda
ENV CONDA_BIN=${CONDA_DIR}/bin
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -qq update && \
    apt-get -qq upgrade -y && \
    apt-get -qq install --yes --no-install-recommends \
        gettext-base \
        gnupg \
        locales \
        wget \
        curl \
        bzip2 \
        ca-certificates \
        ssh \
        file \
        git \
        awscli \
        rsync \
        sudo \
        openssh-server \
        htop \
        screen \
     > /dev/null && \
    apt-get -qq purge && \
    apt-get -qq clean && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir -p /run/sshd && \
    chmod 755 /run/sshd && \
    echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8
ENV SHELL /bin/bash
ENV NB_USER=jovyan
ENV NB_UID=1000
ENV USER=${NB_USER}
ENV HOME=/home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER} && \
    echo "$NB_USER ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/notebook && \
    mkdir -p ${APP_BASE} && \
    chown 1000:1000 -R /opt/conda && \
    chown -R $NB_USER:$NB_USER ${APP_BASE}
USER ${NB_USER}

COPY install-miniconda.bash /tmp/install-miniconda.bash
RUN bash /tmp/install-miniconda.bash
COPY install-jupyter.bash /tmp/install-jupyter.bash
COPY environment.yml /tmp/environment.yml
RUN  bash /tmp/install-jupyter.bash && \
    echo '' > ${CONDA_DIR}/conda-meta/history && \
    ${CONDA_BIN}/conda config --system --add channels conda-forge && \
    ${CONDA_BIN}/conda config --system --set auto_update_conda false && \
    ${CONDA_BIN}/conda config --system --set show_channel_urls true
ENV NPM_DIR ${APP_BASE}/npm
ENV NB_PYTHON_PREFIX ${CONDA_DIR}/envs/saturn
ENV PATH ${NB_PYTHON_PREFIX}/bin:${CONDA_BIN}:${NPM_DIR}/bin:${HOME}/.local/bin:${PATH}
WORKDIR ${HOME}

USER root
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin


COPY requirements.txt /tmp/requirements.txt
RUN /opt/conda/bin/conda install -n saturn python=3.7 pip ipykernel
RUN /opt/conda/envs/saturn/bin/pip install --no-cache-dir -r /tmp/requirements.txt
# install AWS CLI
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip -q awscliv2.zip
RUN ./aws/install
RUN /opt/conda/envs/saturn/bin/python -m ipykernel install --name python3 --display-name 'saturn (Python 3)' --prefix=/opt/conda

USER ${NB_USER}
