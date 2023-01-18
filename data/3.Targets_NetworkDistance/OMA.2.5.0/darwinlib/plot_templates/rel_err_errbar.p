set term postscript eps enhanced color "Helvetica" 24
set pointsize 2
set output "<DVAR name="plotfile">"
set title "<DVAR name="title">"
set xlabel "<DVAR name="xlabel">" 0,.2
set ylabel "<DVAR name="ylabel">" 1,0
set key top left

plot \
"<DVAR name="datafile">" using 1:3:(1.96*sqrt($4)) title "Rel. Error Var Matrix" with errorbars lt 1,\
"<DVAR name="datafile">" using 1:5:(1.96*sqrt($6)) title "Rel. Error with Var/Cov Matrix" with errorbars lt 3



