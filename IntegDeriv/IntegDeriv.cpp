#include <iostream>
#include<cmath>
#include <vector>
#include <random>
#include <functional>
using namespace std;
#include "IntegDeriv.h"

IntegDeriv::IntegDeriv(Functor& f) : F(f) {}


double formula_simpson(double xi, double xf,int n, Functor& f){
    double h = (xf- xi)/n;
    vector<double> valores_x; 
    for (int i = 0; i <= n; i++){
        valores_x.push_back(xi + h*i);
    }

    double primeiro_somatorio = 0;
    for (int j=1; j<=n/2 -1; j++){
        primeiro_somatorio+=f(valores_x[2*j]);
    }

    double segundo_somatorio=0;
    for(int j=1; j<=n/2; j++){
        segundo_somatorio+=f(valores_x[2*j - 1]);
    }
    return (h)*(f(valores_x[0])+ 2*primeiro_somatorio + 4* segundo_somatorio + f(valores_x[n])) / 3;
}

void IntegDeriv::simpsonRule(double xi, double xf, double& Integral, double& Error){
    //CORRIGIR ERRO
    double max_erro = 0.001;
    double erro = 1+max_erro;
    int n=2;
    double integral1=0;
    double integral2=1;
    while (fabs(erro)>max_erro)
    {
        integral1 =formula_simpson(xi, xf, n,F);
        //integral2 = formula_simpson(xi,xf,n*2, F);
        //erro = integral2-integral1;
        double passo=(xf-xi)/n;
        erro=((xf-xi)*n*n*n*n*fourth_derivative((xf+xi)/2, 0.1))/180;
        n=n*2;
    }
    Error=erro;
    Integral = integral1;
}


void IntegDeriv::TrapezoidalRule(double xi, double xf, double& Integral, double& Error){
    double max_erro=0.001;
    double erro = 1+ max_erro;
    int n=2;
    double integral1=0;
    while (fabs(erro)>max_erro)
    {
        double h=(xf-xi)/n;
        integral1=0;
        double somatorio_int=0;
        double somatorio_der=0;
        for(int i=0; i<n; i++){
            double x1=xi+i*h;
            double x2=x1+h; 
            somatorio_int+=(F(x1)+F(x2));

            //ERRO_SEGUNDA_DERIVADA
            somatorio_der+=second_derivative((x1+x2)/2, h/2);
        }
        integral1=somatorio_int*(h/2);
        double media_derivadas=somatorio_der/n;
        erro=((h*h)*(xf-xi)*abs(media_derivadas))/12;
        n=n*2;
    }
    Error=erro;
    Integral=integral1;
}


void IntegDeriv::TrapezoidalRule_adptive_scheme(double xi, double xf, double& Integral, double& Error){
    //FALTA CORRIGIR O ERRO, FAZER AS DERIVADAS
    double max_erro = 0.001;
    double erro = 1+max_erro;
    int n=2;
    double h = xf-xi;
    double integral1=(h/2)*(F(xi)+F(xf));
    double integral2=1;
    int n_integrais=1;
    while(fabs(erro)>max_erro){
        n_integrais+=1;
        integral2=integral1;
        integral1 = integral1/2;
        double somatorio=0;
        for (int i=1; i<=pow(2, n_integrais-2);i++){
            somatorio+=F(xi+(2*i-1)*(h/pow(2, n_integrais-1)));
        }
        integral1+=(h/pow(2, n_integrais-1))*somatorio;
        //NAO HA EXPRESSOES DE DERIVADAS PARA ESTE, POR ISSO FICA ASSIM
        erro=integral2-integral1;
    }
    Error=erro;
    Integral = integral1;
}

void IntegDeriv::midpointRule(double xi, double xf, double& Integral, double& Error){
    double erro_max=0.001;
    double erro=1+erro_max;
    int n=2;
    double integral1=0;
    while(fabs(erro)>erro_max){
        double h = (xf-xi)/n; //parcela
        integral1=0;
        double somatorio_derivadas=0;
        for(int i=0; i<n; i++){
            double x1=xi+i*h;
            double x2=x1+h; 
            integral1+=h*F((x1+x2)/2);

            //ERRO_SEGUNDA_DERIVADA
            somatorio_derivadas+=second_derivative((x1+x2)/2, h/2);
        }
        double media_derivadas=somatorio_derivadas/n;
        erro=((h*h)*(xf-xi)*abs(media_derivadas)) / 24;
        n=n*2;
    }
    Error=erro;
    Integral=integral1;
}


/*void IntegDeriv::simpsonRule(double xi, double xf, double& Integral, double& Error){
    //perceber onde entre a quarta derivada
    int n=2;
    double erro_max=0.0001;
    double erro=1+erro_max;
    double integral1=0;
    while (fabs(erro)>erro_max)
    {   
        double integral2=integral1;
        double h=(xf-xi)/n;
        integral1=F(xi)+F(xf);
        for(int i=1; i<n-1; i++){
            if(i%2!=0){
                integral1+=4*F(xi+i*h);
            }
            else{
                integral1+=2*F(xi+i*h);
            }
        }
        integral1=integral1*(h/3);
        erro=integral1-integral2;
        n=n*2;
    }
    Integral=integral1;
    Error=erro;
}
*/

//###################################################################
//DERIVADAS

double IntegDeriv::first_derivative_five_points(double x, double h){
    double resultado = F(x-2*h)-8*F(x-h)+8*F(x+h)-F(x+2*h);
    return (1./(12*h))*resultado;
}

double IntegDeriv::second_derivative(double x, double h){
    double resultado = -F(x-2*h)+16*F(x-h)-30*F(x)+16*F(x+h)-F(x+2*h);
    return (1./(12*h*h))*resultado;
}

double IntegDeriv::forward_first_derivative(double x, double h){
    double diference = F(x+h)-F(x);
    return diference/h;
}

 double IntegDeriv::central_first_derivative(double x, double h){
    double diference = F(x+h)-F(x-h);
    return diference/(2*h);

 }

 double IntegDeriv::backward_first_derivative(double x, double h){
    double diference = F(x)-F(x-h);
    return diference/h;
 }

 double IntegDeriv::third_derivarive(double x, double h){
    double diference = second_derivative(x+h,h) - second_derivative(x,h);
    return diference/h;
 }

 double IntegDeriv::fourth_derivative(double x, double h){
    double diference = third_derivarive(x+h, h) - third_derivarive(x+h,h);
    return diference/h;
 }

/*double IntegDeriv::third_derivarive(double x, double h){
  double result = (1/(12*h*h))*(-first_derivative_five_points(x-2*h) +16*first_derivative_five_points(x-h) 
  -30*first_derivative_five_points(x) + 16*first_derivative_five_points(x+h) - first_derivative_five_points(x+2*h));
  return result;
}*/

/*double IntegDeriv::fourth_derivative(double x, double h){
  double result = (1/(12*h*h))*(-second_derivative(x-2*h, h) +16*second_derivative(x-h, h) //aumenta a precisao
  -30*second_derivative(x, h) + 16*second_derivative(x+h, h) - second_derivative(x+2*h, h));
  return result;
}*/




//################################################################################
//Zeros Funcoes
double IntegDeriv::bisection_method(double a, double b, double tolerancia_max){
    if (F(a)*F(b)>0){
        cout << "Impossível calcular"<<endl;
        return 999;
    }

    double x1=a;
    double x2=b;
    double x_meio=(x2+x1)/2.;

    //se eu multiplicar por x_meio, quanod aquilo se aproxima de 0, o while explode
    while(fabs(F(x_meio))>tolerancia_max){
        x_meio=(x1+x2)/2.;
        double f1=F(x1);
        double f2=F(x2);
        double fr = F(x_meio);
        if(f1*F(x_meio)<=0){
            x2=x_meio;
        }
        else{
            x1=x_meio;
        }
    }
    return x_meio;
}


double IntegDeriv::regula_falsi(double a, double b, double tolerancia_max){
    if (F(a)*F(b)>0){
        cout << "Impossível calcular"<<endl;
        return 999;
    }

    double x1=a;
    double xr = x1 - F(x1)*((b-x1)/(F(b)-F(x1)));
    while(fabs(F(xr))>=tolerancia_max){
        x1=xr;
        xr = x1 - F(x1)*((b-x1)/(F(b)-F(x1)));
    }
    return xr;
}

double IntegDeriv::newton_raphson(double a, double tolerancia_max){
    double x1 = a;
    double xr = x1 -(F(x1)/first_derivative_five_points(x1, 0.001));

    while(fabs(F(xr))>= tolerancia_max){
        x1 = xr;
        xr = x1 -(F(x1)/first_derivative_five_points(x1, 0.001));
    }
    return xr;

}

double IntegDeriv::secant_method(double a, double tolerancia_max){
    double x1=a;
    double x2 = a - F(a)/first_derivative_five_points(a, 0.001);
    double xr = x2;

    while(fabs(F(xr))>=tolerancia_max){
        xr = x2 - (F(x2)*(x2-x1)/(F(x2)-F(x1)));
        x1=x2;
        x2=xr;
    }   
    return xr;
}


//#####################################################################
//MÉTODO DE INTEGRACAO DE MONTE CARLO

void IntegDeriv::integrateMC(double xi, double xf, double& Integral, double& Relative_Error, int& Ngen){
    double F1 = 0; //para calcular a soma da funcao nos varios pontos
    double f1_var = 0; //para calcular a variancia
    minstd_rand R1(time(0));
    uniform_real_distribution D1(0.,1.);

    for (int i=0; i<Ngen; i++){
        double x = xi + (xf-xi)*D1(R1);
        //podia incluir o "TRandom3.h" e usar gRandom->Uniform()

        F1+=F(x);
        f1_var += F(x)*F(x);
    }

    double f1media = F1/Ngen;  //calcula a <f(x)> (valor medio)

    //Calcular integral
    Integral = f1media*(xf-xi);

    //Calcular variancia
    double Var1 = f1_var/Ngen-f1media*f1media;

    //Calcular o erro
    double Erro = (xf-xi)/sqrt(Ngen)*sqrt(Var1);

    Relative_Error = Erro/Integral;
}


void IntegDeriv::integrateMC_IS(double xi, double xf,function<double(double)> px, function<double(double)> xy,double& Integral, double& Relative_Error, int& Ngen){
    //Não percebo este método
}


//Escolher máximo no von_Neumann
double maximo(function<double(double)> F, double xi, double xf, double N){
    double maximo1 = 0;
    double h = (xf-xi)/N;
    for (double i=xi; i<xf; i+=xf){
        if(F(i)>maximo1){
            maximo1=F(i);
        }
    }
    return maximo1+1; //tem de ser maior que o maximo da funcao
}

void IntegDeriv::integrate_von_Neumann(double xi, double xf, double& Integral, double& Relative_Error, int& Ngen){
    double F1=0;
    double fmax = maximo(F, xi, xf, 10000);

    minstd_rand R1(time(0));
    minstd_rand R2(time(0)+1);
    uniform_real_distribution D1(xi,xf); //para gerar o valor de x
    uniform_real_distribution D2(0., fmax); //para gerar o valot de y
    double A_retangulo = (xf-xi)*fmax;

    double Nr_dentro=0; //numero de pontos que vao cair dentro da area
    for(int i=0; i<Ngen;i++){
        double x=D1(R1);
        double y = D2(R2);
        if(y<=F(x)){
            Nr_dentro+=1;
        }
    }

    Integral = (xf-xi)*fmax*(Nr_dentro/Ngen);

    double Erro = (A_retangulo / Ngen) * sqrt(Nr_dentro*(1-Nr_dentro/Ngen));
    Relative_Error = Erro/Integral;
}