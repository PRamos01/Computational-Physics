#include <iostream>
#include <random>
using namespace std;

#include "Rwalk1D.h"

Rwalk1D::Rwalk1D(int num_part, double x, double prob_L, double prob_R,double tempo, double distance){
    N=num_part;
    x0=x;
    pL=prob_L;
    pR=prob_R;
    dt=tempo;
    dx=distance;
    for(int i=0; i<N; i++){
        mT[i].push_back(x0);
    }
}

Rwalk1D::~Rwalk1D(){};

void Rwalk1D::Run(int nsteps){
    minstd_rand R1(3765376);
    uniform_real_distribution<double> D1(0.,1.);
    for (int part=0;part<N;part++){
        for (int i=1; i<=nsteps; i++){
            double aleatorio=D1(R1);
            if (aleatorio<=pL){
                mT[part].push_back(mT[part][i-1]-dx);
            }
            else{
                mT[part].push_back(mT[part][i-1]+dx);
            }
        }
        for(int i=0; i<mT[part].size(); i++){
        }
    }
}

double Rwalk1D::GetTimeStep(){
    return dt;
}

double Rwalk1D::GetSpaceStep(){
    return dx;
}

const vector<double>& Rwalk1D::GetTrajectory(int n){
    return mT[n];
}