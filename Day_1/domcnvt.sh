/bin/rm 01_MC_NVT.aux
/bin/rm 01_MC_NVT.bbl
/bin/rm 01_MC_NVT.blg
/bin/rm 01_MC_NVT.dvi
/bin/rm 01_MC_NVT.log
/bin/rm 01_MC_NVT.pdf
/bin/rm 01_MC_NVT.ps
latex -halt-on-error 01_MC_NVT.tex
#only continue if 01_MC_NVT.tex is error free
if (( $? == 0 )); then
	latex 01_MC_NVT.tex
	bibtex 01_MC_NVT
	bibtex 01_MC_NVT
	latex 01_MC_NVT
	latex 01_MC_NVT
	dvips -o 01_MC_NVT.ps 01_MC_NVT.dvi
#	ps2pdf 01_MC_NVT.ps
fi
