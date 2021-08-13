FROM rstudio/r-base:3.5.3-xenial

RUN apt-get update \
 && apt-get install -y \
    libmagick++-dev=8:6.8.9.9-7ubuntu5.16 \
    libwebp-dev \
    cargo \
    software-properties-common

RUN add-apt-repository -y ppa:cran/ffmpeg-3 \
 && apt-get update

RUN add-apt-repository -y ppa:cran/poppler \
 && apt-get update

RUN apt-get install -y \
	libavfilter-dev \
	libpoppler-cpp-dev \
	libtesseract-dev \
	libleptonica-dev \
	tesseract-ocr-eng

RUN R -e "install.packages('magick', dependencies=TRUE, repos='http://cran.rstudio.com/')"
