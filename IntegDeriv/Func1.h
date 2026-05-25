#include <functional>
#include "Functor.h"

#ifndef __FUNC1__
#define __FUNC1__

class Func1 :public Functor
{
    public:
        //using FuncaoLambda = function<double(double)>;
        Func1();
        ~Func1()=default;

        double operator()(double x);
        void Draw(array<float,2> range, int npoits, array<string,2> titles);
    
    private:
    //FuncaoLambda f;
        
};


#endif