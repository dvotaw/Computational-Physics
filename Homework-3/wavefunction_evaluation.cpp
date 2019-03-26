//  File: wavefunction_evaluation.cpp
// 
//  Program to evaluate the approximate ground state wavefunction for a Coulomb potential as a function of matrix dimension.
//  The program calculates the sum over i of (|exact_i| - |approximate_i|)/|exact_i|, where i runs over the discrete steps in r.
//  The results seem to indicate that this quantity decreases linearly as the matrix dimension is increased.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      03/26/2019  Written from scratch.
///******************************************************************
// Header files.
#include <iostream> // For cin().
#include <cmath>    // For exp().
#include <fstream>  // For file I/O.
using namespace std;

// Calculates exact hydrogen atom ground state wavefunction.
inline double exact(double r) { return 2. * r * exp(-r); }

// Main function.
int main()
{
  // Doubles to hold the r coordinate, the approximate wavefunction, and the exact wavefunction.
  double r, approx, ex;

  // Variable holding a running sum of the relative deviation.
  double rel_dev;

  // Array holding the dimension sizes.
  int N[] = {1, 5, 10, 20};

  // Input file name.
  char fName[100] = {};

  // Input file.
  ifstream fin;

  // Output file.
  ofstream fout;
  fout.open("wavefunction_evaluation.dat");
 
  // Loop over dimensions. 
  for(int i = 0; i < 4; i++)
  {
    // Reset the running sum.
    rel_dev = 0.;

    // Open the input file.
    sprintf(fName, "wavefunction_N%i_10b8.dat", N[i]);
    fin.open(fName, ios::in);

    // Read values from the file.
    while(fin.good())
    {
      fin >> r >> approx;
      // Calculate exact result.
      ex = exact(r);
      // Compare approximate result to exact result. Omit nodes.
      if(r > 0) { rel_dev += (fabs(ex) - fabs(approx))/fabs(ex); }
    }

    // Write results to file.
    fout << N[i] << " " << rel_dev << endl;

    // Close input file.
    fin.close();
  }

  // Close output file.
  fout.close();

  // Return main().
  return 0;
}
