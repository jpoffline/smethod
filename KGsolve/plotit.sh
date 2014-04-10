#!/bin/sh

####################################################
#
#   J. Pearson, Durham, March 2014
#
#
####################################################
#
# USEAGE
#
# To run for first time, need to change permissions
#  chmod +x plotit.sh
#
# To run normally, call with dir & fileID as arguments
#   (NOTE: dont need forward-slash on dir argument)
#
#  ./plotit.sh DIR FILEID
#
# NOTES
#
# The outputted figure will be named
# "DIR/FILEID_plot.eps"
#  - i.e. in same dir, with same fileID prefix
#
####################################################

echo "Start plotting"
gnuplot << EOF


# INPUT DATA FILE
inputdata="$1/$2.dat"

# OUTPUT IMAGE FILE
outputfig="$1/$2_plot.eps"

# Change column numbers here...

posloc=1
fldloc=2

# Do the plotting
set term postscript landscape color enhanced 20
set output outputfig


set xlabel "x"
set ylabel "field"


plot inputdata u posloc:fldloc  w l lw 2 notitle


EOF

echo "Done plotting"
epstopdf --autorotate=All $1/$2_plot.eps
echo "Conversion to PDF complete"

