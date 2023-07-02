set term postscript eps enhanced color "<DVAR name="fontname" default="Helvetica">" <DVAR name="fontsize" default="24">
set pointsize <DVAR name="pointsize" default="2">
set output "<DVAR name="plotfile">"
if ("<DVAR name="title" default="">" != "")  set title "<DVAR name="title" default="">" ; else unset title;
set xlabel "<DVAR name="xlabel">" 0,.2
set ylabel "<DVAR name="ylabel">" 1,0
set xrange [<DVAR name="xmin" default="*">:<DVAR name="xmax" default="*">]
set yrange [<DVAR name="ymin" default="*">:<DVAR name="ymax" default="*">]
#if ("<DVAR name="legend" default="">" != "") set key top left; else unset key;
set pm3d map
set palette <DVAR name="colorpalette" default="color">

splot "<DVAR name="datafile">" u 1:2:3 notit
