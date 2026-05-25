#include <iostream>
#include <map>
#include <vector>
using namespace std;

#ifndef __RWALK1D_H__
#define __RWALK1D_H__

class Rwalk1D {
    public:
        Rwalk1D(int N=1, double x=0., double pL=0.5, double pR=0.5,double dt=1, double dx=1);

        ~Rwalk1D();

        // particle simulation
        void Run(int nsteps); // number of steps

        // getters
        const vector<double>& GetTrajectory(int n=1); // particle number
        double GetTimeStep();
        double GetSpaceStep();


    private:
        double x0; // init coo
        int N; // number of particles
        double pL, pR; // probabilities (left, same, righ)
        double dt, dx; // steps (time, space)
        map<int, vector<double> > mT; // trajectories
};

#endif