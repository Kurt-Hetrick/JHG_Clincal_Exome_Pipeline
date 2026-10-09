FROM ubuntu:16.04

ENV ORAD_URL https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen-decompression/orad.2.7.0.linux.tar.gz
ENV ORAD_MD5 = "3aba41b14e48ec670e20a073edf7736a"

SHELL ["/bin/bash", "-c"]

RUN apt-get update -qq && \
	apt-get install -qq \
	bzip2 \
	gcc \
	g++ \
	make \
	python \
	zlib1g-dev \
	wget

RUN wget -nv $ORAD_URL && \
	tar -xzvf orad.2.7.0.linux.tar.gz && \
	cd orad.2.7.0.linux && \
	mv orad /usr/bin && \
	cd .. && \
	rm -rf orad.2.7.0.linux
