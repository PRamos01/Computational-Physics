#include <iostream>
#include "ODEsolver.h"
using namespace std;

//ATENCAO: DIFERENTE DOS OUTROS CODIGOS PORQUE O MEU RECEBE UM VETOR COMO ARGUMENTO


ODEsolver::ODEsolver(vector<function<double(ODEpoint)>> f){
    F=f;
};


void ODEsolver::Euler(vector<ODEpoint>& V, ODEpoint(P0), double DT, double step){
    V.push_back(P0);
    double t=0;
    while(t<DT){
        ODEpoint Pcur = V.back();
        vector<double> var_atualiz;
        for(int i=0; i<F.size();i++){
            var_atualiz.push_back(Pcur.X()[i]+step*F[i](Pcur));
        }
        V.push_back(ODEpoint(t+step, var_atualiz));
        t=V.back().T();
    }
    MS["Euler"]=V;
}

void ODEsolver::Trapezio(vector<ODEpoint>& V, ODEpoint(P0), double DT, double step){
    //é o predictor
    V.push_back(P0);
    double t=0;
    while(t<DT){
        ODEpoint Pcur = V.back();
        vector<double> var_next;
        for(int i=0; i<F.size(); i++){
            var_next.push_back(Pcur.X()[i]+step*F[i](Pcur));
        }
        ODEpoint Pnext(t+step,var_next);

        vector<double> var_incremento;
        for(int j=0; j<F.size();j++){
            var_incremento.push_back(Pcur.X()[j]+0.5*step*( F[j](Pcur) + F[j](Pnext)));
        }

        V.push_back(ODEpoint(t+step, var_incremento));
        t=V.back().T();
    }
    MS["Trapezoidal"]=V;
}

void ODEsolver::Leap_Frog(vector<ODEpoint>& V, ODEpoint P0, double DT, double step){
    double t=0;
    V.push_back(P0);
    vector<double> v_primeiro;
    for(int i=0; i<F.size(); i++){
        v_primeiro.push_back(P0.X()[i]+step*F[i](P0));
    }
    V.push_back(ODEpoint(t+step, v_primeiro));

    while(t<DT){
        ODEpoint Pcur = V.back();
        ODEpoint Pantes= V[V.size()-2];

        vector<double> var_next;
        for(int i=0; i<F.size(); i++){
            var_next.push_back(Pantes.X()[i]+2*step*F[i](Pcur));
        }
        V.push_back(ODEpoint(t+step, var_next));
        t=V.back().T();
    }
    MS["LeapFrog"]=V;
}

void ODEsolver::RK2(vector<ODEpoint>& V, ODEpoint P0, double DT, double step){
    double t=0;
    V.push_back(P0);
    while(t<DT){
        ODEpoint Pcur = V.back();
        vector<double> v_int;
        for(int i=0; i<F.size(); i++){
            v_int.push_back(Pcur.X()[i]+(step/2)*F[i](Pcur));
        }
        ODEpoint Pint (t+step/2, v_int);

        vector<double> v_final;
        for(int i=0; i<F.size(); i++){
            v_final.push_back(Pcur.X()[i]+(step)*F[i](Pint));
        }
        V.push_back(ODEpoint(t+step,v_final));
        t=V.back().T();
    }
    MS["RK2"]=V;
}

void ODEsolver::RK4(vector<ODEpoint>& V, ODEpoint P0, double DT, double step){
    double t=0;
    V.push_back(P0);
    while(t<DT){
        ODEpoint Pcur = V.back();

        vector<double> v_int1;
        for(int i=0; i<F.size(); i++){
            v_int1.push_back(Pcur.X()[i]+(step/2)*F[i](Pcur));
        }
        ODEpoint Pint1(t+step/2, v_int1);

        vector<double> v_int2;
        for(int i=0; i<F.size(); i++){
            v_int2.push_back(Pcur.X()[i]+(step/2)*F[i](Pint1));
        }
        ODEpoint Pint2(t+step/2,v_int2);

        vector<double> v_int3;
        for(int i=0; i<F.size(); i++){
            v_int3.push_back(Pcur.X()[i]+(step)*F[i](Pint2));
        }
        ODEpoint Pint3(t+step, v_int3);

        vector<double>  v_Pfinal;
        for(int i=0;i<F.size(); i++){
            v_Pfinal.push_back(Pcur.X()[i] + (step/6)*(F[i](Pcur)+2*F[i](Pint1)+2*F[i](Pint2)+F[i](Pint3)));
        }

        V.push_back(ODEpoint(t+step, v_Pfinal));
        t=V.back().T();
    }
    MS["RK4"]=V;
}

/*void ODEsolver::Euler_Cromer(vector<ODEpoint>& V, ODEpoint P0, double DT,double step){
    V.push_back(P0);
    double t=0;
    while(t<DT){
        ODEpoint Pcur = V.back();
        double v_seguinte = Pcur.X()[1]+step*F[1](Pcur);
        V.push_back(ODEpoint(t+step, { Pcur.X()[0]+step*v_seguinte, v_seguinte}));
        t=V.back().T();
    }
    MS["Euler_Cromer"]=V;
}*/
