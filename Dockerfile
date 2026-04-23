ARG R_VERSION=4.4
ARG MINC_VERSION=1.9.18
FROM rocker/r-ver:${R_VERSION}

ARG MINC_VERSION

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    freeglut3-dev \
    cmake \
    git \
    lsof \
    autoconf \
    automake \
    libtool \
    wget \
    python3 \
    python3-pip \
    python3-dev \
    imagemagick \
    libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

RUN wget -q \
    https://packages.bic.mni.mcgill.ca/minc-toolkit/Debian/minc-toolkit-1.9.18-20220625-Ubuntu_20.04-x86_64.deb \
    && dpkg -i minc-toolkit-*.deb || apt-get install -f -y \
    && rm minc-toolkit-*.deb

ENV MINC_PATH=/opt/minc/${MINC_VERSION}
ENV PATH="${MINC_PATH}/bin:${PATH}"
RUN echo "${MINC_PATH}/lib" > /etc/ld.so.conf.d/minc.conf && ldconfig

RUN pip3 install --no-cache-dir --break-system-packages pyminc

RUN git clone --recursive \
    https://github.com/Mouse-Imaging-Centre/minc-stuffs.git \
    /tmp/minc-stuffs \
    && cd /tmp/minc-stuffs \
    && ./autogen.sh \
    && ./configure --with-build-path=/opt/minc/${MINC_VERSION} --prefix=/usr/local \
    && make -j32 \
    && make install \
    && pip3 install --no-cache-dir --break-system-packages . \
    && rm -rf /tmp/minc-stuffs

RUN install2.r --error --skipinstalled \
    BiocManager \
    batchtools \
    dplyr \
    tidyr \
    readr \
    lme4 \
    purrr \
    shiny \
    Rcpp \
    Matrix \
    tibble \
    yaml \
    data.tree \
    visNetwork \
    rjson \
    DT \
    rlang \
    bigstatsr \
    testthat \
    rgl \
    plotrix \
    lmerTest \
    gridBase \
    igraph \
    devtools \
    && Rscript -e "BiocManager::install('qvalue', ask = FALSE, update = FALSE)" \
    && rm -rf /tmp/downloaded_packages

RUN mkdir -p /RMINC

COPY run-tests.sh /usr/local/bin/run-tests.sh
RUN chmod +x /usr/local/bin/run-tests.sh

WORKDIR /RMINC
ENTRYPOINT ["/usr/local/bin/run-tests.sh"]
