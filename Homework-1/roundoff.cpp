//  File: roundoff.cpp
//
//  This program demonstrates roundoff errors in computational math.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      28-Jan-2019  Written from scratch.
//
//  Notes:  
//   * compile with:  "g++ -o roundoff.x roundoff.cpp"
//
//*********************************************************************// 

/// Answers to questions: ///
/*
  The relative error was zero at least up to N = 50, so the plot begins at N = 100, because log(0) is singular.

  For N < 10^2 the relative error is zero, then it fluctuates between N ~ 10^2 and 10^4, then it appears to obey a power law between N ~ 10^4 and 10^7. From a linear least-squares fit,
  the power appears to be ~ 2. For N > 10^7 the relative error saturates to a constant value because the terms in the sums become smaller than the machine precision,
  so continuing to add them doesn't change the sums.

  The downward sum is more precise because it sums up the small values first, then the larger ones. So there is no precision lost due to roundoff.
  The upward sum involves a loss of precision because at large N, it's summing numbers of very different orders of magnitude.
*/

// Header files.
#include <iostream> // For cout, cin.
#include <iomanip>  // For I/O manipulators.
#include <fstream>  // For file output.
#include <cmath>    // For fabs().

using namespace std;

// Function for summing up.
inline float sum_up(int N)
{
  float sum = 0.;
  for(int i = 1; i <= N; i++)
  {
    sum += 1. / float(i);
  }
  return sum;
}

// Function for summing down.
inline float sum_down(int N)
{
  float sum = 0.;
  for(int i = N; i >= 1; i--)
  {
    sum += 1. / float(i);
  }
  return sum;
}

int main()
{
  // Define some floats. "upsum" and "downsum" will hold results of summing up and down. "relative_error" will be the relative error.
  float upsum, downsum, relative_error;

  // Output file.
  ofstream outfile;
  outfile.open("roundoff_output.txt", ios::out);

  // Loop over many values of N.
  for(int n = 100; n <= 1e9; n *= 2)
  {
    // Sum up.
    upsum = sum_up(n);
    // Sum down.
    downsum = sum_down(n);

    // Calculate the relative error.
    relative_error = fabs(upsum - downsum) / (0.5 * (fabs(upsum) + fabs(downsum)));

    // Print to screen.
    cout << setprecision(20) << "\nN = " << n << ", Up sum = " << upsum << ", Down sum = " << downsum << ".";

    // Write to file.
    outfile << setprecision(20) << log10(n) << " " << log10f(relative_error) << "\n";
  }

  cout << endl;
  // Close output file.
  outfile.close();

  return 0;
}

//*********************************************************************//
