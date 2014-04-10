#!/bin/sh

################################################
# Put the name of the output directory and plot here

# Whats the directory containing the datafiles?
outdir=output

# Whats the plot going to be called?
plotname=plot_comp

loc=2

################################################

for var in "$@"
do
    echo $var
    filename=$filename' '\"''$outdir'/'$var'.dat'\"' u 1:'$loc' w l lw 2 t '\"''$var''\"','
done

# Need to remove final character from filename
filename=${filename%?}

gnuplot << EOF

set term postscript landscape color enhanced 20
set output "$outdir/$plotname.eps"

set yr [-1.1:1.1]
set key spacing 1.5
set key top left
plot $filename

EOF

echo "Done plotting"
epstopdf --autorotate=All $outdir/$plotname.eps
echo "Conversion to PDF complete"
