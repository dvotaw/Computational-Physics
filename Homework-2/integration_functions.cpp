//  File: integration_functions.cpp
//
//  Function definitions for various numerical integration methods.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      19-Feb-2019  Wrote from scratch.
//
//  Notes:  
//
//*********************************************************************// 
#include "integration_functions.h"
#include <gsl/gsl_integration.h>

// Function to evaluate the integral of the function "integrand" between x = xmin and x = xmax numerically, with Npts integration points, using Simpson's rule.
double Simpson(int Npts, double xmin, double xmax, double (*integrand)(double x))
{
  double h = ((xmax - xmin) / double(Npts - 1));
  double x, sum = 0.;		 

  // Sum over even points.   
  for(int i = 2; i < Npts; i += 2)
  {
    x = xmin + h * double(i - 1);
    sum += (4./3.) * h * integrand(x);
  }
  // Sum over odd points.
  for(int i = 3; i < Npts; i += 2)
  {
    x = xmin + h * double(i - 1);
    sum += (2./3.) * h * integrand(x);
  }    
  // Endpoints. 
  sum += (h/3.) * (integrand(xmin) + integrand(xmax));
   
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
    x = xmin + h * double(i - 1);
    sum += (64./45.) * h * integrand(x);
  }
  // Sum over first odd points.
  for(int i = 3; i < Npts; i += 4)
  {
    x = xmin + h * double(i - 1);
    sum += (24./45.) * h * integrand(x);
  }
  // Sum over second odd points.
  for(int i = 5; i < Npts; i += 4)
  {
    x = xmin + h * double(i - 1);
    sum += (28./45.) * h * integrand(x);
  }
  // Endpoints. 
  sum += (14.*h/45.) * (integrand(xmin) + integrand(xmax));
   
  return sum;
}

// Function to evaluate the integral of the function "integrand" between x = xmin and x = xmax numerically, with Npts integration points, using GSL adaptive integration with singularities.
double GSL_integration(int Npts, double xmin, double xmax, double (*integrand)(double x, void* params))
{
  gsl_integration_workspace* work_ptr = gsl_integration_workspace_alloc(Npts);

  const double abs_error = 1.0e-8;
  const double rel_error = 1.0e-8;
  double result, error;	

  gsl_function my_func;
  my_func.function = integrand;

  gsl_integration_qags(&my_func, xmin, xmax, abs_error, rel_error, Npts, work_ptr, &result, &error);


  gsl_integration_workspace_free(work_ptr);

  return result;
}
