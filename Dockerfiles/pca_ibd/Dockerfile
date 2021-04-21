FROM ubuntu:16.04

ENV VCFTOOLS_URL https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz

ENV PLINK2_URL https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210416.zip

ENV KING_URL http://people.virginia.edu/~wc9c/KING/executables/Linux-king19.tar.gz

SHELL ["/bin/bash", "-c"]

RUN apt-get update -qq && \
	apt-get install -qq \
	wget \
	bzip2 \
	gcc \
	g++ \
	make \
	zlib1g-dev \
	pkg-config \
	zip

RUN mkdir /software
WORKDIR /software

# I'm not going to delete the source code directory for vcftools b/c of the perl scripts.

RUN wget -nv VCFTOOLS_URL && \
	tar -xzf vcftools-0.1.16.tar.gz && \
	cd vcftools-0.1.16 && \
	./configure && \
	make && \
	make install

RUN wget -nv $PLINK2_URL && \
	unzip plink_linux_x86_64_20210416.zip && \
	mv plink /usr/local/bin/ && \
	rm -rf LICENSE plink_linux_x86_64_20210416.zip prettify toy.map toy.ped

RUN wget -nv $KING_URL && \
	tar -xzf Linux-king19.tar.gz && \
	mv king /usr/local/bin/ && \
	rm -rf Linux-king19.tar.gz
