#  File: Tmax_analysis.plt 
#
#  GNUplot plot file for analyzing the fidelity of SA as a function of the initial temperature.
#  
#  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
# 
#  Revision history:
#   29-Apr-2019: Copied from a previous .plt file and modified.
#

# Record time and date.
set timestamp

# Titles and labels.
set title 'Parameter analysis for simulated annealing'
set xlabel 'Tmax'
set ylabel 'Fidelity (fraction of runs that converged to global minimum)'

# Legend placement.
set key left

# Set terminal type to X11.
set term x11

# Plot data with fit lines.
plot "results_Tmax.txt" using ($1):($2) title "Simulated annealing" with lines

# Postscript output.   
set out "Tmax_analysis.ps"
set term postscript color enhanced
replot
