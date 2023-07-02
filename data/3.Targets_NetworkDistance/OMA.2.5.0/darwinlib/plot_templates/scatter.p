set term postscript eps color "<DVAR name="font" default="Helvetica">" <DVAR name="fontsize" default="18">
set output "<DVAR name="plotfile" default="temp.ps">"
set title "<DVAR name="title" default="">"
set xlabel "<DVAR name="xlabel" default="">" 0,0 
set ylabel "<DVAR name="ylabel" default="">" 0,0
set xrange [<DVAR name="xmin" default="*">:<DVAR name="xmax" default="*">]
set yrange [<DVAR name="ymin" default="*">:<DVAR name="ymax" default="*">]


set key <DVAR name="legendPos" default="top left">

plot "<DVAR name="datafile">" u <DVAR name="x" default="1">:<DVAR name="y" default="2"> \
  notit w points lt 1 lw 1 ps <DVAR name="pointsize" default="1">

