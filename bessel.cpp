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

/// Answers to questions: ///
/*
  When the relative error is close to 1, it means that one or both of the recursion methods is failing horribly.
  By comparison with the GSL function, up recursion is failing for x < 2, and down recursion is failing for x > 40.
  Between x = 2 and x = 10, the error appears to obey a power law with a negative slope.
  From x = 10 until about x = 30, the error is on the same order of magnitude as the machine precision for doubles.
  From x = 30 until about x = 40, the error appears to obey a different power law, now with a positive slope.
*/

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
  for(int k = m; k > 0; k--)
  {
    j[k - 1] = ((2. * double(k) + 1.) / x) * j[k] - j[k + 1];
  }
  double scale = (sin(x) / x) / j[0];
  return (j[n] * scale);
}

// Up recursion function.
double up_recursion(double x, int n)
{
  double term_three = 0.;
  double term_one = sin(x) / x; 
  double term_two = (sin(x) - x * cos(x)) / (x * x);
  for(int k = 1; k < n; k++)
  {
    term_three = ((2. * double(k) + 1.) / x) * term_two - term_one;	       
    term_one = term_two;
    term_two = term_three;
  }
  return term_three;
}

int main()
{
  // Define some doubles to hold the values of the Bessel function using 3 different methods, and the relative error between the up and down methods.
  double ans_down, ans_up, relative_error, ans_gsl;

  // Create an output file.
  ofstream my_out("bessel.dat");

  // Write a description line to the output file.
  my_out << "# Spherical Bessel functions via up and down recursion, relative error, and function from GSL\n"; 

  // Step over x values.
  for(double x = xmin; x <= xmax; x += step)
  {
    // Calculate J(x) using down recursion.
    ans_down = down_recursion(x, order, start);
    // Calculate J(x) using up recursion.
    ans_up = up_recursion(x, order);
    // Calculate relative error between up and down methods.
    relative_error = fabs(ans_up - ans_down) / (fabs(ans_up) + fabs(ans_down));
    // Get the value of J(x) from the GSL function.
    ans_gsl = gsl_sf_bessel_jl(order, x);

    // Write everything to file.
    my_out << fixed << setprecision(10) << setw(13) << x << " "
    << scientific << setprecision(10)
    << setw(13) << ans_down << " "
    << setw(13) << ans_up << " "
    << setw(13) << relative_error << " "
    << setw(13) << ans_gsl
    << endl;
  }
  cout << "\nData stored in bessel.dat.\n\n";

  // Close the output file.
  my_out.close();
  return 0;
}
