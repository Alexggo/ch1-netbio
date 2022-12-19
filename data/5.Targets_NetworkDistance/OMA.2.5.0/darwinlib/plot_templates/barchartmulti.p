set term postscript eps color "<DVAR name="font" default="Helvetica">" <DVAR name="fontsize" default="18">
set output "<DVAR name="plotfile" default="temp.ps">"
set title "<DVAR name="title">"
set ylabel "<DVAR name="ylabel">" 1,0 
set xrange [0.4:<DVAR name="nSeries">.6]
set boxwidth 1 

set xtics nomirror rotate (<DVAR name="xtics"> )
set format y "%3.2f"
set key Left reverse
set border 3

<DLOOP name="styles">
set style line <DLVAR name="nr" lt 1 lc rgb "<DLVAR name="color">" lw 2
</DLOOP>

set xtics nomirror rotate by -45 (<DVAR name="xtics"> )
set format y "%3.2f"
unset border
set key screen 0.5, screen 0.05 center horizontal
set label "<DVAR name="small_label">" at screen 0.99, screen 0.02  right font "Helvetica,9" textcolor rgbcolor "gray1"

<DLOOP name="mins">
set label "-" at <DLVAR name="offset">, graph 0.02 center font "Helvetica,32" front textcolor rgbcolor "white"
</DLOOP>

#series below
plot \
<DLOOP name="series"> \
    '<DVAR name="datafile">' i <DLVAR name="ind">:<DLVAR name="ind"> using 1:2 \
        notit w boxes fill solid lt <DLVAR name="lt">, \
    '<DVAR name="datafile">' i <DLVAR name="ind">:<DLVAR name="ind"> using 1:2:3
        notit w yerrorbars lt 0 pt -1, \
    '<DVAR name="datafile">' i <DLVAR name="ind">:<DLVAR name="ind"> using ($1+1):4 \
        notit w boxes fill solid lt <DLVAR name="lt2">, \
    '<DVAR name="datafile">' i <DLVAR name="ind">:<DLVAR name="ind"> using ($1+1):4:5 \
        notit w yerrorbars lt 0 pt -1, \
</DLOOP> \
      f(x) = 0
