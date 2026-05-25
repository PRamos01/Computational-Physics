#include <iostream>
#include <vector>
#include <functional>
#include <map>
using namespace std;
#include "ODEpoint.h"

class ODEsolver
{
public:
    ODEsolver()=default;
    ODEsolver(vector<function<double(ODEpoint)>> f);
    ~ODEsolver()=default;
    void Euler(vector<ODEpoint>& V, ODEpoint P0, double DT, double step);
    void Trapezio(vector<ODEpoint>& V, ODEpoint P0, double DT, double step);
    void Leap_Frog(vector<ODEpoint>& V, ODEpoint P0, double DT,double step);
    void RK2(vector<ODEpoint>& V, ODEpoint P0, double DT, double step);
    void RK4(vector<ODEpoint>& V, ODEpoint P0, double DT, double step);
    //void Euler_Cromer(vector<ODEpoint>& V, ODEpoint P0, double DT, double step);
private:
    vector<function<double(ODEpoint)>> F;
    map<string, vector<ODEpoint> > MS; // key: "euler", "trapezoidal",
};
