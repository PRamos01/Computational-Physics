#include <iostream>
#include <vector>
using namespace std;

#ifndef __ODEPOINT_H__
#define __ODEPOINT_H__

class Xvar {
public:
    Xvar() = default;
    Xvar(int); // number of dependent variables
    Xvar(vector<double>); // passing vector
    // using initializer list to build object: X({1,2})
    ~Xvar()=default;
    Xvar(const Xvar&); // copy constructor
    Xvar& operator=(const Xvar&); // assignment operator
    Xvar& operator+(const Xvar&); // operator+
    double& operator[](int); // X[i]
    friend Xvar operator*(double, const Xvar&); // scalar*X
    friend ostream& operator<< (ostream&, const Xvar&);
    vector<double>& X(); // accessor to x

protected:
    vector<double> x;
    int n_var;
};



class ODEpoint: public Xvar{
public:
    ODEpoint();
    ODEpoint(double t_, Xvar a_);
    ODEpoint(double, vector<double>);
    ~ODEpoint()=default;

    //vector<double>& X();
    double& T();
    
    friend ostream& operator<< (ostream& s, ODEpoint& P);

    void SetODEpoint(double t_, Xvar& p);
    void SetODEpoint(double t_, std::vector<double> v);

private:
    double time;
    vector<double> var;
};
#endif