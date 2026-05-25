#include <vector>
#include <string>
#include "Functor.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TAxis.h"

TCanvas* Functor::c = nullptr; 

Functor::Functor(string name_): name(name_){
        c=new TCanvas("c","canvas Functor", 1000, 1000);
}


double Functor::operator()(double x){return 0;}

void Functor::Draw(array<float,2> range, int npoints, array<string, 2> titles){

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

    //TApplication A("A", 0,0);
    G.Draw("AP");
    c->SetGrid();
    c->Update();
    G.SetTitle(name.c_str());
    G.GetXaxis()->SetTitle(titles[0].data());
    G.GetYaxis()->SetTitle(titles[1].data());
    c->SaveAs("FIG_GRAFICO.pdf");

    gSystem ->ProcessEvents(); //global variables
    c ->WaitPrimitive();
}