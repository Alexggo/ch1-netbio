set term postscript eps color "<DVAR name="font" default="Helvetica">" <DVAR name="fontsize" default="18">
set output "<DVAR name="plotfile" default="temp.ps">"
set title "<DVAR name="title" default="">"
set xlabel "<DVAR name="xlabel" default="">" 0,0 
set ylabel "<DVAR name="ylabel" default="">" 0,0 

<DLOOP name="series">
set style line <DLVAR name="ls"> lt rgb "#<DLVAR name="hexcolor">" lw 1 \
  pt <DLVAR name="pointtype" default="1"> ps <DLVAR name="pointsize" default="2">
</DLOOP>

<DBOOL name="legendPos" true="set key <DVAR name="legendPos" default="top left">" false="unset key">

plot \
<DLOOP name="series"> \
  "<DVAR name="datafile">" using <DLVAR name="x">:<DLVAR name="y"> \
    tit '<DLVAR name="label">' w points ls <DLVAR name="ls">, \
  <DBOOL name="xerr" true=" \ 
    '' using <DLVAR name="x">:<DLVAR name="y">:<DLVAR name="xerr"> w xerrorbars ls <DLVAR name="ls">,"> \
  <DBOOL name="yerr" true=" \
    '' using <DLVAR name="x">:<DLVAR name="y">:<DLVAR name="yerr"> w yerrorbars ls <DLVAR name="ls">,"> \
</DLOOP> \
  f(x)=0
