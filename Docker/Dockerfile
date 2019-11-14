FROM continuumio/miniconda:4.7.12
MAINTAINER Irina Kuznetsova

RUN apt-get -y install ttf-dejavu

COPY conda.yml .
RUN \
   conda env update -n root -f conda.yml \
&& conda clean -a

RUN apt-get install -y procps
