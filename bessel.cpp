//  file: bessel.cpp 
//
//  Spherical Bessel functions via up and down recursion, compared with GSL.      
//                                                                     
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      29-Jan-2019  Taken from course website.
//
//  Notes:  
//   * Compile with: g++ -Wall -o bessel.x bessel.cpp -lgsl -lgslcblas
//   * Data saved as: x, yDown, yUp, Relative Error, yGSL                        
//  
//************************************************************************

// Header files.
#include <iostream>            // For cout.
#include <iomanip>             // For setwidth() and setprecision().
#include <fstream>             // For file output.
#include <cmath>               // For fabs().
#include <gsl/gsl_sf_bessel.h> // For spherical Bessel function.

using namespace std;


// Global constants.
const double xmax = 100.0; // Maximum x. 
const double xmin = 0.1;   // Minimum x. Avoid singularity at x = 0.  
const double step = 0.1;   // Step size in x.
const int order = 10;      // Spherical Bessel function order.
const int start = 50;	   // Starting point for down recursion. 

//********************************************************************

// Down recursion function.
double down_recursion(double x, int n, int m)
{
  double j[start + 2];
  j[m + 1] = j[m] = 1.;	
  for (int k = m; k > 0; k--)
  {
    j[k - 1] = ((2. * double(k) + 1.) / x) * j[k] - j[k + 1];
  }
  double scale = (sin(x) / x) / j[0];
  return j[n] * scale;
}

// Up recursion function.
double up_recursion(double x, int n)
{
  double term_three = 0.;
  double term_one = sin(x) / x; 
  double term_two = (sin(x) - x * cos(x)) / (x * x);
  for (int k = 1; k < n; k += 1)
  {
    term_three = ((2. * double(k) + 1.) / x) * term_two - term_one;	       
    term_one = term_two;
    term_two = term_three;
  }
  return term_three;
}

int main()
{
  double ans_down, ans_up, relative_error, ans_gsl;

  // Create an output file.
  ofstream my_out("bessel.dat");

  my_out << "# Spherical Bessel functions via up and down recursion\n"; 

  // Step over x values.
  for(double x = xmin; x <= xmax; x += step)
  {
    ans_down = down_recursion(x, order, start);
    ans_up = up_recursion(x, order);
    relative_error = fabs(ans_up - ans_down)/(fabs(ans_up) + fabs(ans_down));
    ans_gsl = gsl_sf_bessel_jl(order, x);

    my_out << fixed << setprecision(15) << setw(18) << x << " "
    << scientific << setprecision(15)
    << setw(18) << ans_down << " "
    << setw(18) << ans_up << " "
    << setw(18) << relative_error << " "
    << setw(18) << ans_gsl
    << endl;
  }
  cout << "\nData stored in bessel.dat.\n\n";

  // Close the output file.
  my_out.close();
  return 0;
}
