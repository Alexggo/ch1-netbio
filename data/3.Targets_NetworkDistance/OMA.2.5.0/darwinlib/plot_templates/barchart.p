set term postscript eps enhanced color "Helvetica" 18
set output "<DVAR name="plotfile">"
set title "<DVAR name="title">"
set ylabel "<DVAR name="ylabel">" 1,0 
set xrange [0.4:<DVAR name="nSeries">.6]
set boxwidth 0.70

#set xtics nomirror rotate ( "OMA" 1, "RoundUp strict" 2, "RoundUp sensitive" 
#3, "Inparanoid" 4, "ENSEMBL" 5, "orthoMCL" 6, "KOG" 7, "COG" 8, "Homologene" 
#9, "eggNOG" 10, "BBH 130" 11, "BBH 217" 12)

set xtics nomirror rotate (<DVAR name="xtics"> )
set format y "%3.2f"
set key Left reverse
#set bmargin 10
#set tmargin 10
#set lmargin 10
#set rmargin 10
set border 3

# each series has
plot "<DVAR name="datafile">" using 1:2 notit w boxes fill solid lt 1,\
     "<DVAR name="datafile">" using 1:2:3 notit w yerrorbars lt 0 pt -1
