//  File: eigen_basis.cpp
// 
//  Program to find bound state eigenvalues for various potentials
//   by diagonalizing the Hamiltonian using the GSL eigenvalue/eigenvector 
//   routines in a truncated harmonic oscillator basis.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      03/26/2019  Taken from course website.
//
//  Notes:
//   * Based on the documentation for the GSL library under
//      "Eigensystems" and on Chap. 15 of "Computational Physics"
//      by Landau and Paez.
//   * Uses the GSL functions for computing eigenvalues
//      and eigenvectors of matrices.  The steps for this part are:
//       * define and allocate space for matrices and vectors we need
//       * load the matrix to be diagonalized (pointed to by Hmat_ptr)
//       * find the eigenvalues and eigenvectors with gsl_eigensymmv
//       * sort the results numerically 
//       * print out the results  
//   * As a convention (advocated in "Practical C++"), we'll append
//      "_ptr" to all pointers.
//   * When sorting eigenvalues,
//      GSL_EIGEN_SORT_VAL_ASC => ascending order in numerical value 
//      GSL_EIGEN_SORT_VAL_DESC => descending order in numerical value 
//      GSL_EIGEN_SORT_ABS_ASC => ascending order in magnitude 
//      GSL_EIGEN_SORT_ABS_DESC => descending order in magnitude
//   * We use gls_integration_qagiu for the integrals from
//      0 to Infinity (calculating matrix elements of H).
//   * Start with l=0 (and generalize later)
//
//  To do:
//   * Add the Morse potential (function is given but not incorporated)
//   * Generalize to l>0.
//   * Make potential selection less kludgy.
//   * Improve efficiency (reduce run time)
//   * Split into more files (?) or convert to classes
//
///******************************************************************

/*  /////////////////// Answers to questions: ///////////////////

In order to make the calculation converge, b should be chosen to be large enough to capture the important features of the true wavefunction. 
If it's chosen to be too small, the convergence rate as a function of N will be poor, and you won't be able to reproduce the tail of the true wavefunction.

N should be made as large as possible. If you include too few basis functions, you won't be able to get all of the features right in the approximate wavefunction.
If b is chosen correctly, then choosing too small of an N will mostly affect the large-r behavior of the approximate wavefunction.

Based on the plots, I think b = 0.8 is good enough to get a good approximation of the true hydrogen atom ground state, at least if N > 10. N = 20 does very well.
*/

// Header files.
#include <iostream>		// For cout().
#include <iomanip>		// For manipulators.
#include <cmath>                // For math.
#include <fstream>
using namespace std;

#include <gsl/gsl_eigen.h>	        // GSL eigenvectors.
#include <gsl/gsl_integration.h>	// GSL integration.

// Structures.
typedef struct			// Structure holding Hij parameters. 
{
  int i;			// 1st matrix index. 
  int j;			// 2nd matrix index. 
  double mass;			// Particle mass. 
  double b_ho;			// Harmonic oscillator parameter. 
  int potential_index;		// Indicates which potential to use. 
}
hij_parameters;

typedef struct			// Structure holding potential parameters. 
{
  double param1;		// Any three parameters. 
  double param2;
  double param3;
}
potential_parameters;

// Potentials. 
double V_coulomb(double r, potential_parameters* potl_params_ptr);
double V_square_well(double r, potential_parameters* potl_params_ptr);
double V_morse(double r, potential_parameters* potl_params_ptr);

// i'th-j'th matrix element of Hamiltonian in HO basis. 
double Hij(hij_parameters ho_parameters);
double Hij_integrand(double x, void* params_ptr);

// Harmonic oscillator routines from harmonic_oscillator.cpp. 
extern double ho_radial(int n, int l, double b_ho, double r);
extern double ho_eigenvalue(int n, int l, double b_ho, double mass);

// To square double precision numbers.
inline double sqr(double x)  {return x*x;}

//************************** Main program ***************************
int main()
{
  hij_parameters ho_parameters;  // Parameters for the Hamiltonian.

  // Will store output file name.
  char fName[100] = {};

  // Limits for plotting wavefunction.
  const double rmin = 0.;
  const double rmax = 10.;
  const double delta_r = 0.05;

  // Pick the potential based on the integer "answer". 
  int answer = 0;
  while(answer != 1 && answer != 2)	// Don't quit until 1 or 2!
    {
      cout << "Enter 1 for Coulomb or 2 for square well potential: ";
      cin >> answer;
    }
  ho_parameters.potential_index = answer;

  // Set up the harmonic oscillator basis. 
  double b_ho;			// HO length parameter. 
  cout << "Enter the oscillator parameter b: ";
  cin >> b_ho;

  double mass = 1;		 // Measure mass in convenient units. 
  ho_parameters.mass = mass;    
  ho_parameters.b_ho = b_ho;

  // Pick the dimension of the basis (matrix).
  int dimension;		// Dimension of the matrices and vectors. 
  cout << "Enter the dimension of the basis: ";
  cin >> dimension;

  // Decide which wavefunction to write to file.
  int whichPrint;
  cout << "Which wavefunction would you like to save to file? (0, 1, 2, ..., N) ";
  cin >> whichPrint;

  // Generate unique filename, and open output file. Use int(10*b) so we don't have to worry about extra digits in the filename.
  sprintf(fName, "wavefunction_N%i_10b%i.dat", dimension, int(10.*b_ho));
  ofstream fout;
  fout.open(fName, ios::trunc);
 
  //  Define and allocate space for the vectors, matrices, and workspace. 
                               // Original GSL matrix with Hamiltonian. 
  gsl_matrix* Hmat_ptr = gsl_matrix_alloc(dimension, dimension); 
                               // GSL vector with eigenvalues. 
  gsl_vector* Eigval_ptr = gsl_vector_alloc(dimension);	
                               // GSL matrix with eigenvectors. 
  gsl_matrix* Eigvec_ptr = gsl_matrix_alloc(dimension, dimension);	
                               // The workspace for GSL. 
  gsl_eigen_symmv_workspace* worksp = gsl_eigen_symmv_alloc(dimension);	

  // Load the Hamiltonian matrix pointed to by Hmat_ptr. 
  for(int i = 0; i < dimension; i++)
  {
    for(int j = 0; j < dimension; j++)
    {
      ho_parameters.i = i;
      ho_parameters.j = j;
      gsl_matrix_set(Hmat_ptr, i, j, Hij(ho_parameters));
    }
  }

  // Find the eigenvalues and eigenvectors of the real, symmetric
  //  matrix pointed to by Hmat_ptr.  It is partially destroyed
  //  in the process. The eigenvectors are pointed to by 
  //  Eigvec_ptr and the eigenvalues by Eigval_ptr.
  gsl_eigen_symmv(Hmat_ptr, Eigval_ptr, Eigvec_ptr, worksp);

  // Sort the eigenvalues and eigenvectors in ascending order. 
  gsl_eigen_symmv_sort(Eigval_ptr, Eigvec_ptr, GSL_EIGEN_SORT_VAL_ASC);

  // Print out the results.   
  // Allocate a pointer to one of the eigenvectors of the matrix. 
  gsl_vector* eigenvector_ptr = gsl_vector_alloc(dimension);	
  for(int i = 0; i < dimension; i++)
  {
    double eigenvalue = gsl_vector_get(Eigval_ptr, i);
    gsl_matrix_get_col(eigenvector_ptr, Eigvec_ptr, i);

    cout << "eigenvalue " << i+1 << " = " << scientific << eigenvalue << endl;

    // Get desired eigenvector.
    if(i == whichPrint)
    {
      // Dynamically allocate an array of doubles to hold the expansion coefficients.
      double* coef; coef = new double[dimension];

      // Define a temporary variable to hold the value of the wavefunction at a given r.
      double wf;
      
      // Get expansion coefficients from components of desired eigenvector.  
      for(int j = 0; j < dimension; j++)
      {
        coef[j] = gsl_vector_get(eigenvector_ptr, j);
      }

      // Loop over r and calculate the radial wavefunction at each point.
      for(double r = rmin; r <= rmax; r += delta_r)
      {
        wf = 0.;
        // Calculate the wavefunction by summing contributions from each of the HO basis functions.
        for(int k = 0; k < dimension; k++)
        {
          wf += coef[k] * ho_radial(k+1, 0, b_ho, r); 
        }
      // Write to file.
      fout << setprecision(12) << r << " " << wf << endl; 
      }
     // Memory management.
     delete[] coef; 
    }
  }

  // Free the space used by the vector and matrices and workspace. 
  gsl_matrix_free(Eigvec_ptr);
  gsl_vector_free(Eigval_ptr);
  gsl_matrix_free(Hmat_ptr);
  gsl_vector_free(eigenvector_ptr);
  gsl_eigen_symmv_free(worksp);

  // Close output file.
  fout.close();

  // Return main().
  return 0;
}

//************************************************************

//************************** Hij ***************************
//  
// Calculate the i'th-j'th matrix element of the Hamiltonian
//  in a Harmonic oscillator basis.  This routine just passes
//  the integrand Hij_integrand to a GSL integration routine
//  (gsl_integration_qagiu) that integrates it over r from 0
//  to infinity
//
// Take l=0 only for now 
//
//*************************************************************
double Hij(hij_parameters ho_parameters)
{
  gsl_integration_workspace* work = gsl_integration_workspace_alloc(1000);
  gsl_function F_integrand;

  double lower_limit = 0.;	// Start integral from 0 (to infinity). 
  double abs_error = 1.0e-8;	// To avoid round-off problems.
  double rel_error = 1.0e-8;	// The result will usually be much better.
  double result = 0.;		// The result from the integration. 
  double error = 0.;		// The estimated error from the integration.

  void* params_ptr;		// Void pointer passed to function. 

  params_ptr = &ho_parameters;	// We'll pass i, j, mass, b_ho.

  // Set up the integrand. 
  F_integrand.function = &Hij_integrand;
  F_integrand.params = params_ptr;

  // Carry out the integral over r from 0 to infinity. 
  gsl_integration_qagiu(&F_integrand, lower_limit, abs_error, rel_error, 1000, work, &result, &error);

  // Return the result of the integration.
  return result;
}

//************************** Hij_integrand ***************************
// 
// The integrand for the i'th-j'th matrix element of the 
//  Hamiltonian matrix.
//   * uses a harmonic oscillator basis
//   * the harmonic oscillator S-eqn was used to eliminate the
//      2nd derivative from the Hamiltonian in favor of the
//      HO energy and potential.  This was checked against an
//      explicit (but crude) 2nd derivative (now commented).
//
//************************************************************
double Hij_integrand(double x, void* params_ptr)
{
  potential_parameters potl_params;	// Parameters to pass to potential. 
  double Zesq;			// Ze^2 for Coulomb potential. 
  double R, V0;			// Radius and depth of square well. 
  int potential_index;		// Index 1,2,... for potental. 

  int l = 0;			// Orbital angular momentum. 
  int n_i, n_j;			// Principal quantum number (1,2,...). 
  double mass, b_ho;		// Local HO parameters. 
  double hbar = 1.;		// Units with hbar = 1. 
  double omega;			// Harmonic oscillator frequency. 
  double ho_pot;		// Value of HO potential at current x. 

  n_i = ((hij_parameters *)params_ptr)->i + 1;	// n starts at 1. 
  n_j = ((hij_parameters *)params_ptr)->j + 1;
  mass = ((hij_parameters *)params_ptr)->mass;
  b_ho = ((hij_parameters *)params_ptr)->b_ho;
  omega = hbar / (mass * b_ho * b_ho);
  ho_pot = (1. / 2.) * mass * (omega * omega) * (x * x); // HO otential.
  potential_index = ((hij_parameters *)params_ptr)->potential_index;

  // Set up the potential according to potential index. 
  switch(potential_index)
    {
    case 1:			// Coulomb. 
      Zesq = 1.;
      potl_params.param1 = Zesq;
      return (ho_radial(n_i, l, b_ho, x) * (ho_eigenvalue(n_j, l, b_ho, mass) - ho_pot + V_coulomb(x, &potl_params)) * ho_radial(n_j, l, b_ho, x));
      break;
    case 2:			// Square well. 
      R = 1.;
      V0 = 50.;
      potl_params.param1 = V0;
      potl_params.param2 = R;
      return (ho_radial(n_i, l, b_ho, x) * (ho_eigenvalue(n_j, l, b_ho, mass) - ho_pot + V_square_well(x, &potl_params)) * ho_radial(n_j, l, b_ho, x));
      break;
    default:
      cout << "Shouldn't get here!\n";
      return 1;
      break;
    }
}

//************************** Potentials *************************

//************************** V_coulomb ***************************
//
// Coulomb potential with charge Z:  Ze^2/r
//  --> hydrogen-like atom
//
//   Zesq stands for Ze^2.
//
//**************************************************************
double V_coulomb(double r, potential_parameters* potl_params_ptr)
{
  double Zesq = potl_params_ptr->param1;

  return (-Zesq / r);
}

//**************************************************************

//************************* V_square_well **********************
//
// Square well potential of radius R and depth V0.
//
//**************************************************************
double V_square_well(double r, potential_parameters* potl_params_ptr)
{
  double V0 = potl_params_ptr->param1;
  double R = potl_params_ptr->param2;

  if(r < R)
  {
    return (-V0);  // Inside the well of depth V0. 
  }
  else
  {
    return 0.;	   // Outside the well. 
  }
}

//************************************************************

//************************** V_morse ***************************
//
// Morse potential with equilibrium bond length r_eq and potential
//  energy for bond formation D_eq
//
//**************************************************************
double V_morse(double r, potential_parameters* potl_params_ptr)
{
  double D_eq = potl_params_ptr->param1;
  double r_eq = potl_params_ptr->param2;

  return (D_eq * sqr(1. - exp(-(r-r_eq))));
}

//**************************************************************
