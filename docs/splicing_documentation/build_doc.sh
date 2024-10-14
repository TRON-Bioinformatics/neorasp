#!/bin/sh 

command -v Rscript || exit 1 
Rscript -e "library(bookdown); bookdown::render_book()"
