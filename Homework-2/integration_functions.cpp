//  File: integration_functions.cpp
//
//  Function definitions for various numerical integration methods.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      19-Feb-2019  Wrote from scratch.
//      25-Feb-2019  Added some comments.
//
//  Notes:  
//
//*********************************************************************// 
#include "integration_functions.h" // For function declarations.
#include <gsl/gsl_integration.h>   // For GSL QAGS.

// Function to evaluate the integral of the function "integrand" between x = xmin and x = xmax numerically, with Npts integration points, using Simpson's rule.
double Simpson(int Npts, double xmin, double xmax, double (*integrand)(double x))
{
  // Step size.
  double h = ((xmax - xmin) / double(Npts - 1));
  // Declare some variables to hold the current x position and a running sum.
  double x, sum = 0.;		 

  // Sum over even points.   
  for(int i = 2; i < Npts; i += 2)
  {
    // Calculate x.
    x = xmin + h * double(i - 1);
    // Add the contribution from this point to the running sum.
    sum += (4./3.) * h * integrand(x);
  }
  // Sum over odd points.
  for(int i = 3; i < Npts; i += 2)
  {
    // Calculate x.
    x = xmin + h * double(i - 1);
    // Add the contribution from this point to the running sum.
    sum += (2./3.) * h * integrand(x);
  }    
  // Add the contributions from the endpoints. 
  sum += (1./3.) * h * (integrand(xmin) + integrand(xmax));

  // This is our numerical integration result.   
  return sum;
}

// Function to evaluate the integral of the function "integrand" between x = xmin and x = xmax numerically, with Npts integration points, using Milne's rule.
double Milne(int Npts, double xmin, double xmax, double (*integrand)(double x))
{
  double h = ((xmax - xmin) / double(Npts - 1));
  double x, sum = 0.;		 

  // Sum over even points.   
  for(int i = 2; i < Npts; i += 2)
  {
    // Calculate x.
    x = xmin + h * double(i - 1);
    // Add the contribution from this point to the running sum.
    sum += (64./45.) * h * integrand(x);
  }
  // Sum over first odd points.
  for(int i = 3; i < Npts; i += 4)
  {
    // Calculate x.
    x = xmin + h * double(i - 1);
    // Add the contribution from this point to the running sum.
    sum += (24./45.) * h * integrand(x);
  }
  // Sum over second odd points.
  for(int i = 5; i < Npts; i += 4)
  {
    // Calculate x.
    x = xmin + h * double(i - 1);
    // Add the contribution from this point to the running sum.
    sum += (28./45.) * h * integrand(x);
  }
  // Add the contributions from the endpoints.  
  sum += (14./45.) * h * (integrand(xmin) + integrand(xmax));

  // This is our numerical integration result.   
  return sum;
}

// Function to evaluate the integral of the function "integrand" between x = xmin and x = xmax numerically, with Npts integration points, using GSL adaptive integration with singularities.
double GSL_integration(int Npts, double xmin, double xmax, double (*integrand)(double x, void* params))
{
  // Create a pointer to a workspace, and allocate memory.
  gsl_integration_workspace* work_ptr = gsl_integration_workspace_alloc(Npts);

  // Define error limits.
  const double abs_error = 1.0e-8;
  const double rel_error = 1.0e-8;

  // Declare some variables to hold the result and its error.
  double result, error;	

  // Declare a GSL function.
  gsl_function my_func;
  // Set the GSL function equal to the integrand function from the arguments.
  my_func.function = integrand;

  // Perform the GSL QAGS integration.
  gsl_integration_qags(&my_func, xmin, xmax, abs_error, rel_error, Npts, work_ptr, &result, &error);

  // Free up the memomery allocated at the beginning.
  gsl_integration_workspace_free(work_ptr);

  // Return the result of the numerical integral.
  return result;

  // Return the absolute error of the numerical integral.
  //return error;
}
