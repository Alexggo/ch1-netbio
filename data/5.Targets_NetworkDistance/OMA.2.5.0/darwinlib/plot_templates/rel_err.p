set term postscript eps enhanced color "<DVAR name="fontname" default="Helvetica">" <DVAR name="fontsize" default="24">
set pointsize <DVAR name="pointsize" default="2">
set output "<DVAR name="plotfile">"
set title "<DVAR name="title">"
set xlabel "<DVAR name="xlabel">" 0,.2
set ylabel "<DVAR name="ylabel">" 1,0
set key top left

plot \
"<DVAR name="datafile">" using 1:3 title "Rel. Error Var Matrix" with points lt 2,\
"<DVAR name="datafile">" using 1:5 title "Rel. Error with Var/Cov Matrix" with points lt 3



