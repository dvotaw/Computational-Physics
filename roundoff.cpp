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
  for(int n = 100; n <= 1e8; n *= 2)
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
