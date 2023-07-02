set term postscript eps enhanced color "Helvetica" 24
set pointsize 2
set output "<DVAR name="plotfile">"
set title "<DVAR name="title">"
set xlabel "<DVAR name="xlabel">" 0,.2
set ylabel "<DVAR name="ylabel">" 1,0
set key top right

plot \
"<DVAR name="datafile">" using 1:($3-$5):(1.96*sqrt($4+$6)) title "bias" with errorbars lt 3,\
"<DVAR name="datafile">" using 1:(sqrt($7)):(1.96*sqrt($8)) title "std" with errorbars lt 1



