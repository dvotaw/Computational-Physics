//  File: annealer.cpp
//
//  Function definitions for the annealer class.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      14-Mar-2019  Wrote from scratch.
//      19-Mar-2019  Added capability for adaptive accept/reject criteria. Hoping to improve convergence at low dimensions.
//
//  Notes:  
//
//*********************************************************************//
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "annealer.h"

using namespace std;

// We'll need these values here.
extern const unsigned int SEED;
extern const int DIM;
extern const double LATTICE_SPACING;

// Annealer class constructor.
Annealer::Annealer(const double Tmin, const double Tmax, const int Nstays_init, const bool verbosity)
{
  cout << "\n Initializing annealer.\n";
  T_min = Tmin;
  T = T_max = Tmax;
  N_stays = Nstays_init;
  VERBOSE = verbosity;
  count = 0;
  N_rej = 0;
}

// Annealer class destructor.
Annealer::~Annealer()
{
  cout << "\n Destroying annealer.\n";
}

// Annealer cooling schedule.
void Annealer::cooling_schedule(double T_old)
{
  // Implementing a linear cooling schedule with some cooling rate (must be less than 1).
  const double cooling_rate = 0.99;
  T = cooling_rate*T_old;
}

// Adapting the condition for cooling. We want to increase the required number of consecutive rejections as a function of Monte Carlo time.
void Annealer::adapt_rejection(int Nstays_old)
{
  const int adaptation = 0; // 1;
  N_stays = Nstays_old + adaptation;
}

// The annealing function.
double Annealer::anneal(const double* min, const double* max, double* x, double (*cost)(double*), void (*step)(const double*, const double*, double*))
{
  // Seed the uniform PRNG.
  srand(SEED);

  // Cost at the current point.
  double current_value;
  // Cost at the proposed point.
  double proposed_value;
  // Difference.
  double diff;
  // A temporary array to hold the current position while we step to the proposed position.
  double temp[DIM];

  // Initialize the position to some integers, distributed uniformly over the parameter space.
  for(int i = 0; i < DIM; i++)
  {
    x[i] = int(min[i] + (max[i] - min[i])*double(rand())/double(RAND_MAX));
  }

  /// Begin annealing.
  // Keep going until we reach the minimum temperature.
  while(T > T_min)
  {
    // Keep going until we reach the desired number of consecutive rejections.
    while(N_rej <= N_stays)
    {
      // Copy the current position to a temporary holder before we step.
      for(int i = 0; i < DIM; i++) { temp[i] = x[i]; }
      // Gaussian step to propose a new point.
      step(min, max, x);
      // Calculate the cost at the proposed point.
      proposed_value = cost(x);

      // Compare the proposed value to the current value.
      diff = proposed_value - current_value;

      // If we find a better value, take it.
      if(diff <= 0) { current_value = proposed_value; N_rej = 0; if(VERBOSE) { cout << "\n   Difference: " << diff << ". Found a better value, moving."; } }
      else
      {
        // If we find a worse value, take it anyway with a temperature-dependent probability.
        if(double(rand())/double(RAND_MAX) < exp(-diff/T)) { current_value = proposed_value; N_rej = 0; if(VERBOSE) { cout << "\n   Difference: " << diff << ". Found a worse value, moving."; } }
        // Or remain at the current point.
        else 
        {
          N_rej++;
          if(VERBOSE) { cout << "\n   Difference: " << diff << ". Found a worse value, staying."; }
          for(int i = 0; i < DIM; i++) { x[i] = temp[i]; }
        } 
      }

      // Show current temperature and position if verbosity is turned on.
      if(VERBOSE)
      {
        cout << " T: " << T << " x: <";
        for(int i = 0; i < DIM; i++)
        {
          if(i != DIM - 1) { cout << x[i] << ","; }
          else { cout << x[i] << ">"; }
        }
        cout << ".\n";
      }

      // Increment the counter.
      count++;
    }
    // Tell the user that cooling is occurring.
    cout << "\n  After " << count << " iterations, cooling from T = " << T << " to";
    // Cool.
    cooling_schedule(T);
    cout << " T = " << T << ", and adapting from Nstays = " << N_stays << " to";
    // Adapt rejection requirement.
    adapt_rejection(N_stays);
    cout << " Nstays = " << N_stays << ".\n";
    // Reset the number of consecutive rejections.
    N_rej = 0;
  }

  // Print final results. Should be the location of the global minimum.
  cout << "\n  Annealer converged to: <";
  for(int i = 0; i < DIM; i++)
  {
    if(i != DIM - 1) { cout << x[i] << ","; }
    else { cout << x[i] << ">"; }
  }
  cout << ", Minimum cost = " << cost(x) << ".\n";

  // Return this so we can see how many steps we needed.
  return count;
}
