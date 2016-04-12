#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

cp primerv7.tex temp.tex
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7.pdf
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;