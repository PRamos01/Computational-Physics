#ifndef __FUNCTOR__
#define __FUNCTOR
using namespace std;

#include "TCanvas.h"

class Functor
{
    public:
        Functor(string name="");
        ~Functor() = default;

        virtual double operator()(double x);
        virtual void Draw(array<float,2> range, int npoints, array<string, 2> titles);

    protected:
        static TCanvas *c;
        string name;
};


#endif