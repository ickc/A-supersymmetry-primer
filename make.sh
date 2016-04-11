#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

cp primerv7.tex temp.tex
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7.pdf
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;

cp primerv7.tex temp.tex
sed -i '' s/11pt/12pt/ temp.tex
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7_12pt.pdf
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;

cp primerv7.tex temp.tex
sed -i '' s/'%color '// temp.tex
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7_bluelink.pdf
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;

cp primerv7.tex temp.tex
sed -i '' s/11pt/12pt/ temp.tex
sed -i '' s/'%color '// temp.tex
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7_12pt_bluelink.pdf
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;

cp primerv7.tex temp.tex
sed -i '' s/'\\def\\signofmetric{1}'/'\\def\\signofmetric{0}'/ temp.tex
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7_altmetric.pdf
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;

cp primerv7.tex temp.tex
sed -i '' s/11pt/12pt/ temp.tex
sed -i '' s/'\\def\\signofmetric{1}'/'\\def\\signofmetric{0}'/ temp.tex
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7_12pt_altmetric.pdf
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;

cp primerv7.tex temp.tex
sed -i '' s/'%color '// temp.tex
sed -i '' s/'\\def\\signofmetric{1}'/'\\def\\signofmetric{0}'/ temp.tex
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7_bluelink_altmetric.pdf
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;

cp primerv7.tex temp.tex
sed -i '' s/11pt/12pt/ temp.tex
sed -i '' s/'%color '// temp.tex
sed -i '' s/'\\def\\signofmetric{1}'/'\\def\\signofmetric{0}'/ temp.tex
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7_12pt_bluelink_altmetric.pdf
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;