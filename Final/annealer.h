//  File: annealer.h
//
//  Header file for the annealer class.
//
//  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
//
//  Revision history:
//      14-Mar-2019  Wrote from scratch.
//
//  Notes:  
//
//*********************************************************************//
#ifndef ANNEALER_H
#define ANNEALER_H

class Annealer
{
  public:
    // Constructor
    Annealer(const double Tmin, const double Tmax, const int Nstays_init, const bool verbosity);
    // Destructor
    ~Annealer();

    // The cooling schedule of the algorithm.
    void cooling_schedule(double T_old);

    void adapt_rejection(int Nstays_old);

    // Simulated annealing function.
    double anneal(const double* min, const double* max, double* x, double (*cost)(double*), void (*step)(const double*, const double*, double*));

  private:
    // Current temperature of annealer.
    double T;
    // Minimum temperature; sets when the algorithm terminates.
    double T_min;
    // Maximum temperature; this is the initial temperature.
    double T_max;
    // Number of times required to stay at the same point before reducing the temperature.
    int N_stays;
    // Current number of consecutive rejections.
    int N_rej;
    // Generic counter.
    int count;
    // Verbosity.
    bool VERBOSE;
};

#endif 
