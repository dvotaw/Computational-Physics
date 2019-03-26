#  File: wavefunction.plt 
#
#  GNUplot plot file for eigen_basis.x output.
#  
#  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
# 
#  Revision history
#   03/26/2019  Taken from derivative_test.plt and modified.

# Put timestamp on plot.
set timestamp

# Titles and labels.
set title 'Wavefunction comparison'
set xlabel 'r'
set ylabel 'Wavefunction'

# Legend placement.
set key right

# Set scales.
set xrange [0:10]
set yrange [0:2]

# Set terminal to X11.
set term x11

# Plot exact and approximate wavefunctions.
plot 2*x*exp(-x) title 'Exact ground state wavefunction', \
     "wavefunction_N1_10b2.dat" using ($1):($2) title 'N = 1, b = 0.2', \
     "wavefunction_N5_10b2.dat" using ($1):($2) title 'N = 5, b = 0.2', \
     "wavefunction_N10_10b2.dat" using ($1):($2) title 'N = 10, b = 0.2', \
     "wavefunction_N20_10b2.dat" using ($1):(-$2) title 'N = 20, b = 0.2', \
     "wavefunction_N1_10b8.dat" using ($1):($2) title 'N = 1, b = 0.8', \
     "wavefunction_N5_10b8.dat" using ($1):($2) title 'N = 5, b = 0.8', \
     "wavefunction_N10_10b8.dat" using ($1):($2) title 'N = 10, b = 0.8', \
     "wavefunction_N20_10b8.dat" using ($1):($2) title 'N = 20, b = 0.8'

# Save the image.   
set out "wavefunction.ps"
set term postscript color enhanced
replot
