#!/bin/bash
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# primerv7.pdf
cp primerv7.tex temp.tex
## arXiv's hyperref setting
sed -i '' '3 i\
\\usepackage[dvips,hypertexnames=false,pdfpagelabels,bookmarksnumbered,linktocpage=true]{hyperref}
' temp.tex
## arXiv's version is 11pt
sed -i '' s/12pt/11pt/ temp.tex
## LaTeX to PDF
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7.pdf
## Cleanup
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;

# primerv7_12pt.pdf
cp primerv7.tex temp.tex
## arXiv's hyperref setting
sed -i '' '3 i\
\\usepackage[dvips,hypertexnames=false,pdfpagelabels,bookmarksnumbered,linktocpage=true]{hyperref}
' temp.tex
## LaTeX to PDF
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7_12pt.pdf
## Cleanup
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;

# primerv7_bluelink.pdf
cp primerv7.tex temp.tex
## arXiv's hyperref setting plus blue links
sed -i '' '3 i\
\\usepackage[dvips,hypertexnames=false,pdfpagelabels,bookmarksnumbered,linktocpage=true,colorlinks=true,linkcolor=blue,citecolor=blue,filecolor=blue,urlcolor=blue]{hyperref}
' temp.tex
## arXiv's version is 11pt
sed -i '' s/12pt/11pt/ temp.tex
## LaTeX to PDF
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7_bluelink.pdf
## Cleanup
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;

# primerv7_12pt_bluelink.pdf
cp primerv7.tex temp.tex
## arXiv's hyperref setting plus blue links
sed -i '' '3 i\
\\usepackage[dvips,hypertexnames=false,pdfpagelabels,bookmarksnumbered,linktocpage=true,colorlinks=true,linkcolor=blue,citecolor=blue,filecolor=blue,urlcolor=blue]{hyperref}
' temp.tex
## LaTeX to PDF
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7_12pt_bluelink.pdf
## Cleanup
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;

# primerv7_altmetric.pdf
cp primerv7.tex temp.tex
## arXiv's hyperref setting
sed -i '' '3 i\
\\usepackage[dvips,hypertexnames=false,pdfpagelabels,bookmarksnumbered,linktocpage=true]{hyperref}
' temp.tex
## arXiv's version is 11pt
sed -i '' s/12pt/11pt/ temp.tex
## alternative metric (+---)
sed -i '' s/'\\def\\signofmetric{1}'/'\\def\\signofmetric{0}'/ temp.tex
## LaTeX to PDF
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7_altmetric.pdf
## Cleanup
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;

# primerv7_12pt_altmetric.pdf
cp primerv7.tex temp.tex
## arXiv's hyperref setting
sed -i '' '3 i\
\\usepackage[dvips,hypertexnames=false,pdfpagelabels,bookmarksnumbered,linktocpage=true]{hyperref}
' temp.tex
## alternative metric (+---)
sed -i '' s/'\\def\\signofmetric{1}'/'\\def\\signofmetric{0}'/ temp.tex
## LaTeX to PDF
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7_12pt_altmetric.pdf
## Cleanup
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;

# primerv7_bluelink_altmetric.pdf
cp primerv7.tex temp.tex
## arXiv's hyperref setting plus blue links
sed -i '' '3 i\
\\usepackage[dvips,hypertexnames=false,pdfpagelabels,bookmarksnumbered,linktocpage=true,colorlinks=true,linkcolor=blue,citecolor=blue,filecolor=blue,urlcolor=blue]{hyperref}
' temp.tex
## arXiv's version is 11pt
sed -i '' s/12pt/11pt/ temp.tex
## alternative metric (+---)
sed -i '' s/'\\def\\signofmetric{1}'/'\\def\\signofmetric{0}'/ temp.tex
## LaTeX to PDF
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7_bluelink_altmetric.pdf
## Cleanup
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;

# primerv7_12pt_bluelink_altmetric.pdf
cp primerv7.tex temp.tex
## arXiv's hyperref setting plus blue links
sed -i '' '3 i\
\\usepackage[dvips,hypertexnames=false,pdfpagelabels,bookmarksnumbered,linktocpage=true,colorlinks=true,linkcolor=blue,citecolor=blue,filecolor=blue,urlcolor=blue]{hyperref}
' temp.tex
## alternative metric (+---)
sed -i '' s/'\\def\\signofmetric{1}'/'\\def\\signofmetric{0}'/ temp.tex
## LaTeX to PDF
latexmk temp.tex
dvips temp.dvi
ps2pdf temp.ps
mv temp.pdf primerv7_12pt_bluelink_altmetric.pdf
## Cleanup
find . -iname 'temp.*' -exec mv '{}' ~/.Trash/ \;