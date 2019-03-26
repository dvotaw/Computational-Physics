#  File: derivative_test.plt 
#
#  GNUplot plot file for derivative_test.x output.
#  
#  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
# 
#  Revision history
#   03/25/2019  Taken from Session 04 and modified.

# Put timestamp on plot.
set timestamp

# Titles and labels.
set title 'Test of Numerical Derivatives using exp(-x)'
set xlabel 'log10(mesh size)'
set ylabel 'relative error'

# Legend placement.
set key left

# Set scales.
set xrange [-10:0]
set yrange [-17:0]

# Fit functions.
f1(x) = a1*x + b1
fit [-8:-1] f1(x) "derivative_test.dat" using ($1):($2) via a1,b1 
fit_title1 = sprintf("%-+4.1f*x %-+4.1f",a1,b1)

f2(x) = a2*x + b2
fit [-5:-1] f2(x) "derivative_test.dat" using ($1):($3) via a2,b2 
fit_title2 = sprintf("%-+4.1f*x %-+4.1f",a2,b2)

f3(x) = a3*x + b3
fit [-2:-1] f3(x) "derivative_test.dat" using ($1):($4) via a3,b3 
fit_title3 = sprintf("%-+4.1f*x %-+4.1f",a3,b3)

f4(x) = a4*x + b4
fit [-1.2:0] f4(x) "derivative_test.dat" using ($1):($5) via a4,b4
fit_title4 = sprintf("%-+4.1f*x %-+4.1f",a4,b4)

# Set terminal to X11.
set term x11

# Plot data and fit lines.
plot "derivative_test.dat" using ($1):($2) title 'Forward difference', \
     a1*x + b1 title fit_title1, \
     "derivative_test.dat" using ($1):($3) title 'Central difference', \
     a2*x + b2 title fit_title2, \
     "derivative_test.dat" using ($1):($4) title 'First extrapolated difference', \
     a3*x + b3 title fit_title3, \
     "derivative_test.dat" using ($1):($5) title 'Second extrapolated difference', \
     a4*x + b4 title fit_title4

# Save the image.   
set out "derivative_test.ps"
set term postscript color enhanced
replot
