//  file: derivative_test.cpp
// 
//  Program to study the error in differentiation rules.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      03/25/2019  Copied from course website.
//
//  Notes:
//
//*****************************************************************
#include <iostream>       // For cout().
#include <iomanip>        // For manipulators.
#include <fstream>        // For file I/O.
using namespace std;
#include <gsl/gsl_math.h> // For GSL functions.
#include <gsl/gsl_diff.h> // For GSL differentiation.

// Function prototypes. 
double funct(double x, void* params_ptr);
double funct_deriv(double x, void* params_ptr);
double forward_diff(double x, double h, double (*f)(double x, void* params_ptr), void* params_ptr);
double central_diff(double x, double h, double (*f)(double x, void* params_ptr), void* params_ptr);
double extrap_diff(double x, double h, double (*f)(double x, void* params_ptr), void* params_ptr);

// Main function.
int main(void)
{
  void* params_ptr;		// Pointer to void containing function parameters. 

  const double hmin = 1.e-10;	// Minimum mesh size.
  double x = 1.;		// Point x. 
  double alpha = 1.;		// Function parameter. 
  double diff_cd, diff_fd;	// Central and forward differences. 
  double diff_extrap;		// Extrapolated difference. 
  double diff_gsl_cd;		// GSL adaptive central difference. 
  gsl_function My_F;		// GSL function.
  double abserr;                // Absolute error.

  ofstream out ("derivative_test.dat");	// Output file. 

  params_ptr = &alpha;		// Set pointer to the memory address of the function parameter. 

  double answer = funct_deriv(x, params_ptr);	// Exact answer for comparison.

  My_F.function = &funct;	// Create a GSL function with the function we wish to differentiate. 
  My_F.params = params_ptr;     // Set its parameters.
  gsl_diff_central(&My_F, x, &diff_gsl_cd, &abserr);	// GSL calculation.

  cout << "gsl_diff_central(" << x << ") = " << scientific
       << setprecision(16) << diff_gsl_cd << " +/- "
       << setprecision(6) << abserr << endl;
  cout << " actual relative error: " << setprecision(8)
       << fabs((diff_gsl_cd - answer) / answer) << endl;

  double h = 0.1;		// Initialize mesh spacing.
  while(h >= hmin)
    {
      // Calculate the derivative various ways.
      diff_fd = forward_diff(x, h, &funct, params_ptr);
      diff_cd = central_diff(x, h, &funct, params_ptr);
      diff_extrap = extrap_diff(x, h, &funct, params_ptr);

      // Print errors to file. 
      out << scientific << setprecision (8)
	  << log10(h) << "   "
	  << log10(fabs ((diff_fd - answer) / answer)) << "   "
	  << log10(fabs ((diff_cd - answer) / answer)) << "   "
	  << log10(fabs ((diff_extrap - answer) / answer)) << endl;

      h /= 2.;		// Reduce mesh size. 
    }

  out.close();         // Close output file.
  return 0;	       // Return main(). 
}

//************************** funct ***************************
double
funct (double x, void *params_ptr)
{
  double alpha;
  alpha = *(double *) params_ptr;

  return (exp (-alpha * x));
}

//************************** funct_deriv *********************
double
funct_deriv (double x, void *params_ptr)
{
  double alpha = *(double *) params_ptr;

  return (-alpha * exp (-alpha * x));
}

//************************** forward_diff *********************
double
forward_diff (double x, double h,
	      double (*f) (double x, void *params_ptr), void *params_ptr)
{
  return (f (x + h, params_ptr) - f (x, params_ptr)) / h;
}

//************************** central_diff *********************
double
central_diff (double x, double h,
	      double (*f) (double x, void *params_ptr), void *params_ptr)
{
  return (f (x + h / 2., params_ptr) - f (x - h / 2., params_ptr)) / h;
}

//************************** extrap_diff *********************
double
extrap_diff (double x, double h,
	     double (*f) (double x, void *params_ptr), void *params_ptr)
{
  /*
  return (8. * (f (x + h / 4., params_ptr) - f (x - h / 4., params_ptr))
	  - (f (x + h / 2., params_ptr) -
	     f (x - h / 2., params_ptr))) / (3. * h);
  */
  return ( 4.*central_diff (x, h/2., f, params_ptr) -
              central_diff (x, h, f, params_ptr) ) / 3.;	     
}
