//  File: annealer_main.cpp
//
//  Main program for simulated annealing.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      14-Mar-2019  Wrote from scratch.
//      18-Mar-2019  Tuning parameters and testing for different dimensions.
//
//  Notes:  
//
//*********************************************************************//
#include <iostream>   // For cout().
#include <fstream>    // For file I/O.
#include <random>     // For Gaussian distribution.
#include <time.h>     // For clock().
#include "annealer.h" // For Annealer class.

using namespace std;

// PRNG seed.
extern const unsigned int SEED = 775975;

// This is the dimension of the space we're minimizing the cost function in.
extern const int DIM = 7;

// We'll discretize our steps in units of this.
extern const double LATTICE_SPACING = 1.0;

// Flip this Boolean to 1 to see some additional command line output.
const bool VERBOSE = 0;

// Create a RNG.
default_random_engine gen(SEED);

// This parametrizes the step size. It will be the sigma of the Gaussian step distribution.
const double STEP_SIZE = 2.*LATTICE_SPACING;
// Gaussian distribution with zero mean, will be used for stepping. Step will be a DIM-dimensional, uncorrelated Gaussian function. Other functions could be used instead.
normal_distribution<double> gaus(0.0, STEP_SIZE);

// Inline squaring function.
inline double sqr(double x) { return x*x; }
// Prototype for cost function, defined below.
double cost(double* x);
// Prototype for stepping function, defined below.
void step(const double* min, const double* max, double* x);

// Main function.
int main()
{
  /// Create an Annealer object and initialize its parameters.
  // Algorithm will terminate when we get below this temperature.
  const double Tmin = 1e-12;
  // Initial temperature.
  const double Tmax = 2600.;
  // Condition which must be met for cooling to occur.
  const int Nstays_init = 5*DIM;
  Annealer A(Tmin, Tmax, Nstays_init, VERBOSE);

  // Define the boundaries of the space we're minimizing the cost function on.
  double min[DIM];
  double max[DIM];

  // Set the boundaries.
  for(int i = 0; i < DIM; i++)
  {
    min[i] = -5.;
    max[i] = 5.;
  }

  // Declare an array to hold the position within the space.
  double x[DIM] = {};

  // Create a clock for timing.
  clock_t t;
  t = clock();

  // Call the simulated annealing function.
  int counter = A.anneal(min, max, x, cost, step);

  // Calculate the elapsed time.
  t = clock() - t;

  // Calculate how many lattice points there are in our space.
  int volume = 1;
  for(int i = 0; i < DIM; i++)
  {
    volume *= int((max[i] - min[i])/LATTICE_SPACING);
  }

  // Fraction of space we covered in the annealing process. In units of percent.
  double frac = 100.*double(counter)/double(volume);
  // How long the annealer took to converge. In units of seconds.
  double time = double(t)/double(CLOCKS_PER_SEC);

  // Some figures of merit.
  cout << "\n Looped over " << frac << "% of the search space in " << time << " seconds.\n";

  // Write to file for plotting.
  ofstream fout;
  fout.open("results.txt", ios::out | ios::app);
  fout << DIM << " " << frac << " " << time << "\n";

  // Return.
  return 0;
}

// This is the cost function. I'll use a Rosenbrock function in DIM dimensions. 
double cost(double* x)
{
  double ret = 0;
  for(int i = 0; i < DIM - 1; i++)
  {
    ret += 100.*sqr(x[i+1] - sqr(x[i])) + sqr(1. - x[i]);
  }
  return ret;
}

// This is the stepping function which will choose a new point to compare with the current point within the annealer.
void step(const double* min, const double* max, double* x)
{
  // Declare a variable which will hold random numbers drawn from the step distribution.
  double rnd;
  // Loop over dimensions of the space.
  for(int i = 0; i < DIM; i++)
  {
    // Grab a random number.
    rnd = int(gaus(gen))*LATTICE_SPACING;
    // Check whether this step lies within the boundaries.
    while(x[i] + rnd > max[i] || x[i] + rnd < min[i])
    {
      // Redraw until we get a valid step.
      rnd = int(gaus(gen))*LATTICE_SPACING;
    }
    // Update the position with the random step.
    x[i] += rnd;
  }
}
