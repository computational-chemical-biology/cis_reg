FROM continuumio/anaconda3
MAINTAINER Ricardo R. da Silva <ridasilva@usp.br>

RUN  apt-get update
RUN  apt-get -y install build-essential
RUN  apt-get -y install kallisto

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

#RUN conda create -n cis_reg -c bioconda homer meme
RUN conda env create -f environment.yml
#RUN conda activate cis_reg 
RUN echo "source activate cis_reg" > ~/.bashrc
ENV PATH /home/cis_reg/.conda/envs/cis_reg/bin:$PATH

# install the notebook package
RUN pip install --no-cache --upgrade pip && \
    pip install --no-cache notebook jupyterlab

RUN conda config --add channels conda-forge   
RUN conda config --set channel_priority strict      
#RUN conda install -c conda-forge r-base=4.1.3
#RUN conda install -c pkgs/r r-base=4.2.0 --force
RUN /home/cis_reg/.conda/envs/cis_reg/bin/Rscript install_packages.R
#RUN conda install -c bioconda kallisto

