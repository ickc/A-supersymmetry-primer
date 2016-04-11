#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR
latexmk primerv7.tex
dvips primerv7.dvi
ps2pdf primerv7.ps