#  File: roundoff.plt 
#
#  GNUplot plot file for roundoff.x output
#  
#  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
# 
#  Revision history
#   08-Feb-2019 Taken from Session 4 directory and modified
#

# Record timestamp.
set timestamp

# Titles and labels.
set title 'Demonstration of roundoff error'
set xlabel 'log10(N)'
set ylabel 'log10(relative error)'

# Move legend.
set key left

# Set axis ranges.
set xrange [1:9]
set yrange [-7:0]

# Fit line.
f1(x) = a1*x + b1
fit [4:7] f1(x) "roundoff_output.txt" using ($1):($2) via a1,b1 
fit_title1 = sprintf("%-+4.1f*x %-+4.1f",a1,b1)

# Set terminal type to X11.
set term x11

# Plot data with fit line.
plot "roundoff_output.txt" using ($1):($2) title 'Roundoff', a1*x + b1 title fit_title1

# Save an image.
set out "roundoff_output.ps"
set term postscript color enhanced
replot
