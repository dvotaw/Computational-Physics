//  File: roundoff_distribution.cpp
//
//  This program generates a file with many roundoff errors so we can see how they're distributed.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      7-Feb-2019  Written from scratch.
//
//  Notes:  
//   * compile with:  "g++ -o roundoff_distribution.x roundoff_distribution.cpp"
//
//*********************************************************************// 

/*
    The distribution looks approximately Gaussian, but not quite.
*/

// Header files.
#include <iostream> // For cout, cin.
#include <iomanip>  // For I/O manipulators.
#include <fstream>  // For file output.
#include <cmath>

using namespace std;

int main()
{
  // Calculate the same thing with a float and with a double.
  float xf;
  double xd;
  
  double relative_error;
 
  // How high do you want to iterate?
  const int Nmax = 10000;

  // Output file.
  ofstream outfile;
  outfile.open("roundoff_dist.txt", ios::out);

  // Loop over many values of N.
  for(int i = 1; i <= Nmax; i++)
  {
    // Calculate the same thing. Let's try 1/sqrt(x).
    xf = 1./sqrt(float(i));
    xd = 1./sqrt(double(i));

    // Calculate the relative error. No absolute value.
    relative_error = (xd - xf)/xf;

    // Write to file.
    outfile << setprecision(20) << relative_error << "\n";
  }

  cout << endl;
  // Close output file.
  outfile.close();

  return 0;
}

//*********************************************************************//
