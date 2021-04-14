FROM rstudio/r-base:3.5.3-xenial

# tagged github releases do not go back to 2.0
ENV ANNOTSV_URL  https://lbgi.fr/AnnotSV/Sources/AnnotSV_2.0.tar.gz
ENV SAMTOOLS_URL https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2

RUN apt-get update \
 && apt-get install -y \
    libssl-dev \
    bzip2 \
    libncurses5-dev \
    bedtools=2.25.0-1 \
    python=2.7.12-1~16.04 \
    tcl=8.6.0+9 \
    datamash=1.0.7-1

RUN R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "BiocManager::install(c('Biostrings','IRanges','GenomicRanges','Rsamtools','GenomicAlignments'))" \
 && R -e "install.packages('ExomeDepth', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('optparse', dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN cd /tmp \
 && wget -nv $SAMTOOLS_URL \
 && tar -xvf samtools-1.11.tar.bz2 \
 &&	cd samtools-1.11 \
 &&	./configure --prefix=/usr \
 &&	make all all-htslib \
 &&	make install install-htslib \
 &&	cd .. \
 &&	rm -rf samtools-1.11*

RUN mkdir /annotsv \
 && cd /annotsv \
 && wget $ANNOTSV_URL \
 && tar -xzf AnnotSV_2.0.tar.gz \
 && rm AnnotSV_2.0.tar.gz

ENV ANNOTSV /annotsv/AnnotSV_2.0
ENV PATH    "$ANNOTSV/bin:$PATH"

RUN cd $ANNOTSV \
 && echo '-bedtools:     "/usr/bin/bedtools"' > configfile \
 && echo '-genomeBuild:  GRCh37'             >> configfile \
 && echo '-overlap:      70'                 >> configfile \
 && echo '-promoterSize: 500'                >> configfile \
 && echo '-SVinputInfo:  0'                  >> configfile \
 && echo '-SVminSize:	 50'                 >> configfile \
 && echo '-vcfPASS:      0'                  >> configfile
