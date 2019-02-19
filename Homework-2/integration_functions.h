//  File: integration_functions.h
//
//  Header file to store the functions for the various numerical integration methods.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      19-Feb-2019  Wrote from scratch.
//
//  Notes:  
//
//*********************************************************************// 

extern double Simpson(int Npts, double xmin, double xmax, double (*integrand)(double x));

extern double Milne(int Npts, double xmin, double xmax, double (*integrand)(double x));

extern double GSL_integration(int Npts, double xmin, double xmax, double (*integrand)(double x, void* params));
