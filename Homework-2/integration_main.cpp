//  File: integration_main.cpp
//
//  Main program for using various numerical integration methods.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      19-Feb-2019  Wrote from scratch.
//
//  Notes: 
//         Compile with: 
//           make -f make_integration 
//
//*********************************************************************// 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "integration_functions.h"
using namespace std;

double integrand(double x);
double integrand_with_params(double x, void* params);

int main()
{
  const int Npts = 100;
  const double xmin = 0.;
  const double xmax = 1.;
  double result;

  result = Simpson(Npts, xmin, xmax, &integrand);
  cout << endl << result;

  result = Milne(Npts, xmin, xmax, &integrand);
  cout << endl << result;

  result = GSL_integration(Npts, xmin, xmax, &integrand_with_params);
  cout << endl << result;

  cout << endl;
  return 0;
}

double integrand(double x)
{
  const int order = 1;
  double sum = 0.;
  
  for(int i = 0; i <= order; i++)
  {
    sum += pow(x, i);
  }

  return sum;
}

double integrand_with_params(double x, void* params)
{
  return integrand(x);
}
