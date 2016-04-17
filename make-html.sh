#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

cp primerv7.tex temp.tex

# String substitution to replace beq and eeq
sed -i '' s/'\\def\\beq{\\begin{align}}'// temp.tex
sed -i '' s/'\\def\\eeq{\\end{align}}'// temp.tex

sed -i '' s/'\\beq'/'\\begin{eqnarray}'/g temp.tex
sed -i '' s/'\\eeq'/'\\end{eqnarray}'/g temp.tex

# TeX to pandoc md
pandoc -t markdown --normalize -s --wrap=none --atx-headers -o primerv7.md temp.tex

# String substitution: undo pandoc's replacement of eqnarray by aligned
sed -i '' s/'aligned'/'eqnarray'/g primerv7.md

# md to html
pandoc -f markdown+mmd_header_identifiers+abbreviations+autolink_bare_uris+mmd_link_attributes+markdown_attribute+mmd_title_block+tex_math_double_backslash -S --base-header-level=1 --toc --toc-depth=6 -N --normalize -s --mathjax=https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_CHTML-full -t html5 -o primerv7.html primerv7.md

# Cleanup
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;