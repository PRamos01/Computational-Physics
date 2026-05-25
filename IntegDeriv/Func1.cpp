#include <cmath>
#include "Func1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TAxis.h"

//Func1::Func1(FuncaoLambda F) : Functor("func1") {f=F;}

Func1::Func1() : Functor("func1") {}
double Func1::operator()(double x){
    return 1./(1+x*x);
}

void Func1::Draw(array<float,2> range, int npoints, array<string, 2> titles){
    vector<double> x, y;
    double intervalo = range[1] - range[0];
    double dx = intervalo/npoints;
    for (double k = range[0]; k <= range[1]; k+= dx){
        x.push_back(k);
        y.push_back(operator()(k));
    }

    TGraph G(x.size(), x.data(), y.data());
    G.SetMarkerStyle(kFullCircle);
    G.SetMarkerColor(kBlue);
    G.SetTitle(name.c_str());
    G.GetXaxis()->SetTitle(titles[0].data());
    G.GetYaxis()->SetTitle(titles[1].data());

    G.Draw("AP");
    c->SetGrid();
    c->Update();

    c->SaveAs("FIG_GRAFICO.pdf");
    gSystem ->ProcessEvents(); //global variables

    c ->WaitPrimitive();
}