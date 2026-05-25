#include <iostream>
#include <vector>
using namespace std;
#include "ODEpoint.h"
//ATENCAO: DIFERENTE DOS OUTROS CODIGOS PORQUE O MEU RECEBE UM VETOR COMO ARGUMENTO

ODEpoint::ODEpoint():Xvar(){time=0;}


ODEpoint::ODEpoint(double t, vector<double> v): time(t), Xvar(v){var=v;}

//vector<double>& ODEpoint::X(){return var;}


ODEpoint::ODEpoint(double t_, Xvar a_) : time(t_), Xvar(a_) {;}

double& ODEpoint::T(){return time;}

ostream& operator<< (ostream& s, ODEpoint& P){
    s << "[" << P.time <<"; (";
    for(int i=0; i<int(P.X().size()-1); ++i){
        s<<P.X()[i]<<",";
    }
    s<<P.X().back()<<")]";
    return s;
}

void ODEpoint::SetODEpoint(double t_, Xvar& p){
    time=t_;
    Xvar::X() = p.X();
}

void ODEpoint::SetODEpoint(double t_, vector<double> v){
    time=t_;
    var=v;
    Xvar::X() = v;
}

Xvar::Xvar(int n){n_var=n;}

Xvar::Xvar(vector<double> variaveis){
    x=variaveis;
    n_var=variaveis.size();
}

Xvar::Xvar(const Xvar& objeto){
    x=objeto.x;
    n_var=objeto.n_var;   
}

Xvar& Xvar::operator=(const Xvar& objeto){
    x=objeto.x;
    n_var=objeto.n_var;
    return *this;
}

Xvar& Xvar::operator+(const Xvar& objeto){
    for (int i=0; i<x.size(); i++){
        x[i]+=objeto.x[i];
    }
    return *this;
}

double& Xvar::operator[](int n){
    return x[n];
}

Xvar operator*(double escalar, const Xvar& objeto){
    vector<double> variaveis;
    for(int i=0; i<objeto.x.size(); i++){
        variaveis[i]=objeto.x[i]*escalar;
    }
    return Xvar(variaveis);
}

ostream& operator<< (ostream& s, const Xvar& objeto){
    for(int i=0; i<int(objeto.x.size()-1); ++i){
        s<<objeto.x[i]<<",";
    }
    s<<objeto.x.back()<<")]";
    return s;
}

vector<double>& Xvar::X(){
    return x;
} 

