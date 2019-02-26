#  File: integration.plt 
#
#  GNUplot plot file for integration.x output.
#  
#  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
# 
#  Revision history:
#   25-Feb-2019: Copied from Session 4 and modified.
#

# Record time and date.
set timestamp

# Titles and labels.
set title 'Double precision numerical integration of exp(x)cos(x) on [0,1]'
set xlabel 'log10(Number of integration points)'
set ylabel 'log10(Relative error)'

# Legend placement.
set key right

# Set x and y scales.
set xrange [0:8]
set yrange [-16:0]

# Fit functions.
f1(x) = a1*x + b1
fit [0:3.5] f1(x) "integration_results.txt" using ($1):($2) via a1,b1 
fit_title1 = sprintf("%-+4.1f*x %-+4.1f",a1,b1)

f2(x) = a2*x + b2
fit [1.2:2.5] f2(x) "integration_results.txt" using ($3):($4) via a2,b2 
fit_title2 = sprintf("%-+4.1f*x %-+4.1f",a2,b2)

# Set terminal type to X11.
set term x11

# Plot data with fit lines.
plot "integration_results.txt" using ($1):($2) title 'Simpson', \
     a1*x + b1 title fit_title1, \
     "integration_results.txt" using ($3):($4) title 'Milne', \
     a2*x + b2 title fit_title2, \
     "integration_results.txt" using ($5):($6) title 'GSL QAGS'

# Postscript output.   
set out "integration_results.ps"
set term postscript color enhanced
replot
