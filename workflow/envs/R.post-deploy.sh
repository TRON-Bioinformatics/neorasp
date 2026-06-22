#!/bin/bash

set -eou pipefail

TMP_DIR="$(mktemp -d)"
DEPENDENCY_SCRIPT="${TMP_DIR}/install_packages.R"

cat >$DEPENDENCY_SCRIPT <<__EOF__
#!/usr/bin/env Rscript

install.packages(c('argparse', 'remotes'), repos='http://cran.us.r-project.org', Ncpus=8)

library(remotes)
remotes::install_version("Matrix", version = "1.6-5", repos = "http://cran.us.r-project.org")
remotes::install_version("MASS", version = "7.3-60", repos = "http://cran.us.r-project.org")
remotes::install_git("https://github.com/TRON-Bioinformatics/splice2neo.git", ref = "v0.6.12", Ncpus=8)

__EOF__

Rscript --vanilla "$DEPENDENCY_SCRIPT"

rm -r "$TMP_DIR"
