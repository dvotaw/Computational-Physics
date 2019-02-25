//  File: integration_main.cpp
//
//  Main program for using various numerical integration methods.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      19-Feb-2019  Wrote from scratch.
//      25-Feb-2019  Added some comments.
//
//  Notes: 
//      Compile with: 
//        make -f make_integration 
//
//*********************************************************************// 
#include <iostream>                 // For cout().
#include <iomanip>                  // For setprecision().
#include <fstream>                  // For file I/O.
#include <cmath>                    // For fabs().
#include "integration_functions.h"  // For integration function declarations.

using namespace std;

// Function prototype of the integrand function.
double integrand(double x);

// Prototype of the integrand function with a pointer-to-void argument; needed for GSL.
double integrand_with_params(double x, void* params);

// Function prototype containing the exact antiderivative.
double exact_antiderivative(double x);

// Main function.
int main()
{
  // Define the maximum number of intervals we want to use to calculate the numerical integral.
  const int max = 1e8;

  // Define the output precision so it doesn't have to be changed in multiple places.
  const int OUTPUT_PRECISION = 16;

  // Declare a variable to hold the current number of intervals.
  int Npts;

  // Lower integration bound.
  const double xmin = 0.;
  // Upper integration bound.
  const double xmax = 1.;

  // Declare some variables to hold the exact result and the numerical result.
  double result, exact;

  // The exact result for the definite integral is the difference of the antiderivatives evaluated at the endpoints.
  exact = exact_antiderivative(xmax) - exact_antiderivative(xmin);

  // Output file.
  ofstream outfile;
  outfile.open("integration_results.txt");

  // Loop over numbers of integration intervals.
  for(Npts = 10; Npts <= max; Npts *= 10)
  {
    // Print to screen for user-friendliness.
    cout << " Starting " << Npts << " points...\n";

    outfile << setprecision(OUTPUT_PRECISION);

    // Write the number of points to file.
    outfile << fixed << log10(Npts+1);
    // Perform numerical integration using Simpson's rule.
    result = Simpson(Npts+1, xmin, xmax, &integrand);
    // Write to file.
    outfile << scientific << " " << log10(fabs(result - exact) / fabs(exact)) << " ";

    // Write the number of points to file.
    outfile << fixed << log10(Npts+5);
    // Perform numerical integral using Milne's rule.
    result = Milne(Npts+5, xmin, xmax, &integrand);
    // Write to file.
    outfile << scientific << " " << log10(fabs(result - exact) / fabs(exact)) << " ";

    // Write the number of points to file.
    outfile << fixed << log10(Npts);
    // Perform numerical integral using GSL QAGS.
    result = GSL_integration(Npts, xmin, xmax, &integrand_with_params);
    // Write to file.
    outfile << scientific << " " << log10(fabs(result - exact) / fabs(exact));

    // Newline for readability of output file.
    outfile << endl;
  }

  // Close the output file.
  outfile.close();

  // Return main().
  return 0;
}

// This function holds the function that we wish to integrate between x = xmin and x = xmax.
double integrand(double x)
{
  return exp(x) * cos(x);
}

// This function is just the same as the above, but with additional arguments. The arguments are necessary for GSL, but we will not use them.
double integrand_with_params(double x, void* params)
{
  return integrand(x);
}

// This function holds the exact antiderivative of the function in integrand(). We will use it to calculate the exact result, to compare with the numerical results.
double exact_antiderivative(double x)
{
  double ret = exp(x) * (cos(x) + sin(x)) / 2.;
  return ret;
}
