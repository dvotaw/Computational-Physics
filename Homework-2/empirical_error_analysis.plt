#  File: integration.plt 
#
#  GNUplot plot file for integration.x output.
#  
#  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
# 
#  Revision history:
#   26-Feb-2019: Copied from integration.plt and modified.
#

# Record time and date.
set timestamp

# Titles and labels.
set title 'Empirical error analysis between Simpson and Milne'
set xlabel 'log10(Number of integration points)'
set ylabel 'log10(Relative error)'

# Legend placement.
set key right

# Set x and y scales.
set xrange [0:8]
set yrange [-16:0]

# Fit functions.
f1(x) = a1*x + b1
fit [1.2:3.5] f1(x) "integration_results.txt" using ($7):($8) via a1,b1 
fit_title1 = sprintf("%-+4.1f*x %-+4.1f",a1,b1)

f2(x) = a2*x + b2
fit [4:7] f2(x) "integration_results.txt" using ($7):($8) via a2,b2 
fit_title2 = sprintf("%-+4.1f*x %-+4.1f",a2,b2)

# Set terminal type to X11.
set term x11

# Plot data with fit lines.
plot "integration_results.txt" using ($7):($8) title 'Simpson vs. Milne', \
     a1*x + b1 title fit_title1, \
     a2*x + b2 title fit_title2

# Postscript output.   
set out "empirical_error_analysis.ps"
set term postscript color enhanced
replot
