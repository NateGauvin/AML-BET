# build from shiny directory

#FROM openanalytics/r-base
FROM rocker/r-ubuntu:20.04

MAINTAINER Nathaniel M. Gauvin

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    #libcurl4-gnutls-dev \
    libssl-dev \
    libsasl2-dev\
    libssh2-1-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libpq-dev \
    libnlopt-dev \ 
    gfortran \ 
    libblas-dev \
    liblapack-dev \
    unixodbc-dev \
    libmariadb-dev \
    libpng-dev
    
    

RUN  apt install -y --no-install-recommends software-properties-common dirmngr gnupg
RUN  wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc 

RUN  add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" 

RUN  apt-get clean \
    && rm -rf /var/lib/apt/lists/ 

# install required packages

RUN PATH=$PATH:/usr/lib

RUN add-apt-repository ppa:c2d4u.team/c2d4u4.0+
#RUN apt install --no-install-recommends -y r-cran-MASS
#RUN apt install --no-install-recommends -y r-cran-argparse
#RUN Rscript -e 'install.packages(c("ggplot2", "ggpubr" "argparse", "dplyr", "MASS", "tidyverse", "mongolite", "survival", "survminer", "bslib", "shiny", "shinyjs", "shinycssloaders"), repos="https://cloud.r-project.org")'
RUN apt install --no-install-recommends -y r-cran-curl
RUN apt install --no-install-recommends -y r-cran-forestploter
RUN apt install --no-install-recommends -y r-cran-grid
RUN apt install --no-install-recommends -y r-cran-modi
RUN install.r ggplot2 ggpubr argparse dplyr MASS mongolite survival survminer bslib shiny shinyjs shinycssloaders

#RUN R -e "install.packages('GGally', repos = 'https://cloud.r-project.org/')"
#RUN R -e "install.packages('shinyBS', repos = 'https://cloud.r-project.org/')"
#RUN R -e "install.packages('DT', repos = 'https://cloud.r-project.org/')"

#RUN rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# copy the app to the image
COPY . /AML-BET

EXPOSE 3838

CMD ["sudo", "R", "-e", "shiny::runApp('/AML-BET', host = '0.0.0.0', port = 3838)"]