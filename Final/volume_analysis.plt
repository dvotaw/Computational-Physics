#  File: volume_analysis.plt 
#
#  GNUplot plot file for analyzing how many points the annealing algorithm must calculate before convergence.
#  
#  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
# 
#  Revision history:
#   18-Mar-2019: Copied from a previous .plt file and modified.
#

# Record time and date.
set timestamp

# Titles and labels.
set title 'Volume analysis for simulated annealing'
set xlabel 'log10(Dimension)'
set ylabel 'log10(Fraction of space iterated over [%])'

# Legend placement.
set key right

# Set x and y scales.
set xrange [0.4:1]
set yrange [-1:9]

# Fit functions.
f1(x) = a1*x + b1
fit [0.4:1] f1(x) "results_avg.txt" using ($1):($2) via a1,b1 
fit_title1 = sprintf("%-+4.1f*x %-+4.1f",a1,b1)

a2 = 0
b2 = 2
fit_title2 = 'Brute force'

# Set terminal type to X11.
set term x11

# Plot data with fit lines.
plot "results_avg.txt" using ($1):($2) title 'Simulated annealing', \
     a1*x + b1 title fit_title1, \
     a2*x + b2 title fit_title2

# Postscript output.   
set out "volume_analysis.ps"
set term postscript color enhanced
replot
