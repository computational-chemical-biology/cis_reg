FROM continuumio/anaconda3
MAINTAINER Ricardo R. da Silva <ridasilva@usp.br>

RUN conda create -n cis_reg -c bioconda meme 
RUN /opt/conda/bin/activate cis_reg 
RUN echo "source activate cis_reg" > ~/.bashrc
ENV PATH /opt/conda/envs/cis_reg/bin:$PATH

# install the notebook package
RUN pip install --no-cache --upgrade pip && \
    pip install --no-cache notebook jupyterlab

# create user with a home directory
ARG NB_USER=cis_reg
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}
USER ${USER}

COPY . .
RUN conda config --add channels conda-forge   
RUN conda config --set channel_priority strict      
RUN conda install -c conda-forge r-base=4.1.3
RUN Rscript install install_packages.R
RUN conda install -c bioconda/label/cf201901 homer

