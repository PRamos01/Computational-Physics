#include <iostream>
#include <functional>
using namespace std;
#include "Func1.h"

#ifndef __INTEGDERIV__
#define __INTEGDERIV__
class IntegDeriv {
    public:
        IntegDeriv(Functor&);
        ~IntegDeriv() = default;

        //INTEGRAIS
        void TrapezoidalRule(double xi, double xf, double& Integral, double& Error);
        void TrapezoidalRule_adptive_scheme(double xi, double xf, double& Integral, double& Error);
        //Falta definir void Romberg_Integration(double xi, double xf, double& Integral, double& Error);
        void simpsonRule(double xi, double xf, double& Integral, double& Error);
        void midpointRule(double xi, double xf, double& Integral, double& Error);

        //DERIVADAS
        double first_derivative_five_points(double x, double h=0.01);
        double second_derivative(double x, double h);

        double forward_first_derivative(double x, double h);
        double central_first_derivative(double x, double h);
        double backward_first_derivative(double x, double h);

        //SAO AS FOWARDS DERIVATIVES. VER SE E NECESSARIOS IMPLEMENTAR AS OUTRAS 
        double third_derivarive(double x, double h);
        double fourth_derivative(double x, double h);


        //Zeros funcoes
        double bisection_method(double a, double b, double tolerancia_max=0.001); //é necessário ter f(a)*f(b)<0
        double regula_falsi(double a, double b, double tolerancia_max=0.001); 
        double newton_raphson(double a, double tolerancia_max=0.001);
        double secant_method(double a, double tolerancia_max=0.001); 


        //Integração Monte Carlo
        void integrateMC(double xi, double xf, double& Integral, double& Relative_Error, int& Ngen);
        void integrateMC_IS(double xi, double xf,function<double(double)> px, function<double(double)> xy,
                double& Integral, double& Relative_Error, int& Ngen);
        void integrate_von_Neumann(double xi, double xf, double& Integral, double& Relative_Error, int& Ngen);


    private:
        Functor& F;
};

#endif