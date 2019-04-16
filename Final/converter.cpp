//  File: converter.cpp
//
//  Standalone program for combining annealer results.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      16-Apr-2019  Wrote from scratch.
//
//  Notes:  
//	Compile with: g++ -o converter.x converter.cpp
//	Command line argument: the number of data files to combine.
//
//*********************************************************************//
#include <iostream> // For cout().
#include <stdlib.h> // For atoi().
#include <fstream>  // For file I/O.
#include <cmath>    // For log10().
#include <string>   // For output file name.  

using namespace std;

// Main function.
int main(int argc, char* argv[])
{
  // This is the number of results_x.txt files that we want to combine. Naming convention starts from x = 0 to (N - 1).
  int N = 0;

  // C-String buffer.
  char buf[100] = {};

  // I need exactly one argument.
  if(argc != 2) { cerr << "\n Usage: ./converter.x <number of data files>\n" << endl; return 1; }
  else { N = atoi(argv[1]); }

  cout << "\n Analyzing " << N << " data files..." << endl;

  // Going to assume that the input files are structured as: "Dimension Time Volume", with 5 rows of text. We could be smarter about this, but where's the fun in that?
  const int nRows = 5;
  int dimension[nRows];
  double average_time[nRows] = {};
  double average_volume[nRows] = {};

  // Temporary variables for reading from file.
  double tmp_time, tmp_vol;

  // Dynamically allocate some input file objects.
  ifstream* fin = new ifstream[N];

  // Loop over files and average their results. This is a Monte Carlo algorithm, so the results for each individual run are random variables.
  for(int i = 0; i < N; ++i)
  {
    // Generate the file name.
    sprintf(buf, "results_%i.txt", i);
    // Open the file.
    fin[i].open(buf, ios::in);

    // Check that the file exists.
    if(!fin[i].good()) { cerr << "\n I didn't find file " << i << "! Aborting calculation.\n"; return 1; }

    // Loop over the 5 lines that I assume it contains.
    for(int j = 0; j < nRows; ++j)
    {
      // Read each line. Just overwrite dimension[j] each time.
      fin[i] >> dimension[j] >> tmp_time >> tmp_vol;

      // Add to running averages.
      average_time[j] += tmp_time;
      average_volume[j] += tmp_vol;

      // If I'm in the last file, convert my running sum to an average.
      if(i == N - 1)
      {
        average_time[j] /= double(N);
        average_volume[j] /= double(N);
      }
    }
    // Close the file I'm finished with.
    fin[i].close();    
  }

  // Memory management.
  delete[] fin;

  // Write results to an output file.
  ofstream fout;
  const string fName = "results_avg.txt";
  fout.open(fName.c_str(), ios::out | ios::trunc);
  for(int i = 0; i < nRows; ++i)
  {
    fout << log10(dimension[i]) << " " << log10(average_time[i]) << " " << log10(average_volume[i]) << endl;
  }
  fout.close();

  cout << "\n Results written to " << fName << ".\n" << endl;

  // Return normally.
  return 0;
}
