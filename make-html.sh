#!/bin/bash
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

CONTENT=$(<primerv7.tex)

# Turning non-standard LaTeX syntax into standard ones (e.g. def to newcommand, beq to begin eqnarray)
CONTENT=$(echo "$CONTENT" | sed \
	-e s/'\\def\\beq{\\begin{eqnarray}}'// \
	-e s/'\\def\\eeq{\\end{eqnarray}}'// \
	-e s/'\\def\\beq{\\begin{eqnarray*}}'// \
	-e s/'\\def\\eeq{\\end{eqnarray*}}'// \
	-e s/'\\beq'/'\\begin{eqnarray}'/g \
	-e s/'\\eeq'/'\\end{eqnarray}'/g \
	-e 's/\\def\([^{]*\)/\\newcommand{\1}/g' \
	-e s/'\\newcommand{\\thefootnote}'/'\\renewcommand{\\thefootnote}'/g \
	-e s/'\\newcommand{\\G}'/'\\renewcommand{\\G}'/g \
	-e s/'{\\centeron#1#2}'/'{\\centeron}[2]'/g \
	-e s/'{\\slashchar#1}'/'{\\slashchar}[1]'/g)

# ##############################################################################################################
# ## arXiv's hyperref setting
# sed -i '' '3 i\
# \\usepackage[dvips,hypertexnames=false,pdfpagelabels,bookmarksnumbered,linktocpage=true]{hyperref}
# ' temp.tex
# ## arXiv's version is 11pt
# sed -i '' s/12pt/11pt/ temp.tex
# ## LaTeX to PDF
# latexmk temp.tex
# dvips temp.dvi
# ps2pdf temp.ps
# mv temp.pdf primerv7-test.pdf
# mv temp.tex primerv7-test.tex
# ## Cleanup
# find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;
# ##############################################################################################################

# debug
# echo "$CONTENT" > primerv7-normalized.tex

# TeX to pandoc md
echo "$CONTENT" | pandoc -f latex -t markdown --normalize -s --wrap=none --atx-headers -o primerv7.md

# undo pandoc's replacement of eqnarray by aligned
sed -i '' s/'aligned'/'eqnarray'/g primerv7.md

# md to html
pandoc -S --base-header-level=1 --toc --toc-depth=6 -N --normalize -s --mathjax=https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_CHTML-full -t html5 -o primerv7.html primerv7.md