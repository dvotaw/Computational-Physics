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

float sum_up(int N)
{
  float sum = 0.;
  for(int i = 1; i <= N; i++)
  {
    sum += 1. / float(i);
  }
  return sum;
}

float sum_down(int N)
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
  float upsum, downsum, diff;

  ofstream outfile;
  outfile.open("roundoff_output.txt", ios::out);

  for(int n = 100; n <= 1e8; n *= 2)
  {
    upsum = sum_up(n);
    downsum = sum_down(n);

    diff = fabs(upsum - downsum) / (0.5 * (fabs(upsum) + fabs(downsum)));

    cout << setprecision(20) << "\nN = " << n << ", Up sum = " << upsum << ", Down sum = " << downsum << ".";

    outfile << setprecision(20) << log10(n) << " " << log10f(diff) << "\n";
  }

  cout << endl;
  outfile.close();

  return 0;
}

//*********************************************************************//
