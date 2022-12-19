set term postscript eps enhanced color "Helvetica" 24
set pointsize 2
set output "<DVAR name="plotfile">"
set title "<DVAR name="title">"
set xlabel "<DVAR name="xlabel">" 0,.2 
set ylabel "<DVAR name="ylabel">" 1,0 
set key top left

plot <DVAR name="series">

