//  File: area_new.cpp
//
//  This program calculates the area of a circle, given the radius. Input and output can be via command line or file.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      28-Jan-2019  Took original version from PHY 480 course website.
//      28-Jan-2019  Modified original to calculate pi from atan(), use a circle class, file I/O, to be usable multiple times, etc.
//
//  Notes:  
//   * compile with:  "g++ -o area_new.x area_new.cpp"
//
//*********************************************************************// 

// Header files.
#include <iostream> // For cout, cin.
#include <iomanip>  // For I/O manipulators.
#include <cmath>    // For atan().
#include <fstream>  // For file I/O.
#include <string>   // For string class.

using namespace std;

//*********************************************************************//

// Pi is a global constant.
const double pi = 4. * atan(1);

// Circle class. Member data is the radius. Member functions include a constructor and an area calculator.
class circle
{
  public:
  circle(double r);
  inline double calc_area();
  private:
  double radius;
}; 

// Circle class constructor. Checks for negative radii.
circle::circle(double r)
{
  if(r >= 0) // If the input is non-negative, set the data member to its value.
  {
    radius = r;
  }
  else // If the input is negative, warn the user, and make it positive.
  {
    cout << "\n Warning, negative radius! Converting to positive.\n";
    radius = fabs(r);
  }
}

// Area calculator function. A = pi * r^2.
inline double circle::calc_area()
{
  return pi * radius * radius;
}

int main()
{
  // User input to answer command line prompts.
  char answer;
  // User input to decide whether to continue calculations.
  char keep_going = '1';
  // Will be the radius of the circle that the user inputs.
  double radius;
  // Will be the area of the circle for output.
  double area;

  // Precision for output. It's defined here so that it doesn't have to be changed in multiple places.
  const int OUTPUT_PRECISION = 10;

  // Input file.
  ifstream infile;
  // Input file name.
  string ifName;

  // Output file.
  ofstream outfile;
  // Output file name.
  string ofName;

  // Keep going until user wants to stop.
  while(keep_going == '1')
  {
    // Ask user what they want.
    cout << "\nRead from command line (0) or file (1)? ";
    cin >> answer;

    if(answer == '0') // Read from command line.
    {
      cout << "\n Enter radius of the circle: ";
      cin >> radius;
    }
    else if(answer == '1') // Read from file.
    {
      cout << "\n Input the file name as a string: ";
      cin >> ifName;

      infile.open(ifName.c_str(), ios::in);
      infile >> radius;
      infile.close();
    }
    else
    {
      cout << "\n Invalid input! Try again.\n";
      continue;
    }

    // Define an object of the circle class, and initialize it with the user-input radius.
    circle c(radius);
    // Calculate the area of the circle.
    area = c.calc_area();

    // Ask user what they want.
    cout << "\nWrite to command line (0) or file (1)? ";
    cin >> answer;

    if(answer == '0') // Write to command line.
    {
      cout << endl << setprecision(OUTPUT_PRECISION) << "  Radius = " << radius << ", Area = " << area << endl;
    }
    else if(answer == '1') // Write to file.
    {
      cout << "\n Input the file name as a string: ";
      cin >> ofName;

      outfile.open(ofName.c_str(), ios::out | ios::app);
      outfile << setprecision(OUTPUT_PRECISION) << "Radius = " << radius << ", Area = " << area << endl;
      outfile.close();
    }
    else
    {
      cout << "\n Invalid input! Try again.\n";
      continue;
    }

    // Ask user whether to stop or keep going.
    cout << "\nStop (0) or keep going (1)? ";
    cin >> keep_going;
  }

  cout << endl;
  return 0;
}

//*********************************************************************//
