#!/bin/bash

set -eou pipefail

TMP_DIR="$(mktemp -d)"
DEPENDENCY_SCRIPT="${TMP_DIR}/install_packages.R"

cat > $DEPENDENCY_SCRIPT << __EOF__
#!/usr/bin/env Rscript

install.packages('devtools', repos='http://cran.us.r-project.org', Ncpus=8)

library(devtools)
devtools::install_git('https://github.com/gagneurlab/FRASER', dependencies=TRUE, ref='1.99.4', Ncpus=8)

__EOF__

Rscript --vanilla "$DEPENDENCY_SCRIPT"

rm -r "$TMP_DIR"

