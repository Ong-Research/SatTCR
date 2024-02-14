FROM eddelbuettel/r2u:20.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    pandoc \
    pandoc-citeproc \
    curl \
    gdebi-core \
    && rm -rf /var/lib/apt/lists/*

RUN install.r \
    dada2 \
    jsonlite \
    tidyverse \
    ggplot2 \
    cowplot \
    htmltools \
    future \
    furrr \
    fuzzyjoin \
    qs \
    remotes \
    renv \
    knitr \
    rmarkdown \
    quarto

ENV QUARTO_PRINT_STACK=true
ENV XDG_RUNTIME_DIR=./tmp

ARG QUARTO_VERSION="1.4.549"
RUN curl -o quarto-linux-amd64.deb -L https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb
RUN gdebi --non-interactive quarto-linux-amd64.deb

CMD ["bash"]
