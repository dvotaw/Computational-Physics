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
double exact_antiderivative(double x);

int main()
{
  const int max = 1e8;
  int Npts;
  const double xmin = 0.;
  const double xmax = 1.;
  double result, exact;

  exact = exact_antiderivative(xmax) - exact_antiderivative(xmin);

  ofstream outfile;
  outfile.open("integration_results.txt");

  for(Npts = 11; Npts < max; Npts *= 11)
  {
    cout << " Starting " << Npts << " points...\n";
    outfile << scientific << setprecision(20) << log10(Npts);

    result = Simpson(Npts, xmin, xmax, &integrand);
    outfile << " " << log10(fabs(result - exact) / fabs(exact));

    result = Milne(Npts, xmin, xmax, &integrand);
    outfile << " " << log10(fabs(result - exact) / fabs(exact));

    result = GSL_integration(Npts, xmin, xmax, &integrand_with_params);
    outfile << " " << log10(fabs(result - exact) / fabs(exact));

    outfile << endl;
  }

  outfile.close();

  cout << endl;
  return 0;
}

double integrand(double x)
{
  return exp(x) * cos(x);
}

double integrand_with_params(double x, void* params)
{
  return integrand(x);
}

double exact_antiderivative(double x)
{
  double ret = exp(x) * (cos(x) + sin(x)) / 2.;
  return ret;
}
