#include <iostream>
#include <cmath>
using namespace std;
#include "mapas_de_luz.h"

#ifndef mapas_de_luz_cpp
#define mapas_de_luz_cpp

lightmap::lightmap(array<int,2> ncell, array<float,2> size){
    n_cell_x=ncell[0];
    n_cell_y=ncell[1];

    tamanho_x=size[0];
    tamanho_y=size[1];

    //distancia de cada celula nas duas direcões
    dx = tamanho_x/n_cell_x;
    dy = tamanho_y/n_cell_y;

    //criar as celulas e atribuir à grelha
    for (int i=0; i<n_cell_x; i++){
        vector<cell> vet;
        for(int j=0; j<n_cell_y; j++){
            cell celula;
            celula.area=dx*dy; //[cm²]
            celula.center_coo[0]=(i+0.5)*dx; //[cm]
            celula.center_coo[1]=(j+0.5)*dy; //[cm]
            vet.push_back(celula);
        }
        GRID.push_back(vet);
    }
}

lightmap::lightmap(lightsource S, array<int,2> ncell, array<float,2> size){
    n_cell_x=ncell[0];
    n_cell_y=ncell[1];
    
    //adicinar a lightsource ao vetor com as varias luzes
    fotoes.push_back(S);

    tamanho_x=size[0];
    tamanho_y=size[1];

    //distancia de cada celula nas duas direcões
    dx = tamanho_x/n_cell_x;
    dy = tamanho_y/n_cell_y;

    GRID.resize(n_cell_x, vector<cell>(n_cell_y));
    //criar as células e atribuir à grelha
    for (int i=0; i<n_cell_x; i++){
        for(int j=0; j<n_cell_y; j++){
            //definir as características das células
            cell celula;
            celula.area = dx*dy;  //[cm²]
            celula.center_coo[0]=(i+0.5)*dx;  //[cm]
            celula.center_coo[1]=(j+0.5)*dy;  //[cm]
            
            //Calcular a potência [W] provenientes desta luz para todas as células
            //Potência = irradiância * área
            celula.power=0;
            celula.power+= Irradiance(celula.center_coo) * celula.area; // [W]
            GRID[i][j]=celula;  //Acrescentar a célula à  grelha
        }
    }
}

pair<int,int> lightmap::GetCellIndex(float x, float y) const{
    pair<int, int> index;

    //Determinar o index
    int index_x= floor(x/dx);
    int index_y=floor(y/dy);
    index.first=index_x;
    index.second=index_y;
    return index;
}

pair<float,float> lightmap::GetCellCoo(int index_x, int index_y) const{
    pair<float, float> coo;
    coo.first= GRID[index_x][index_y].center_coo[0];
    coo.second = GRID[index_x][index_y].center_coo[1];
    return coo;
}


float lightmap::GetCellPower(int index_x, int index_y) const{
    return GRID[index_x][index_y].power; //[W]
}

float lightmap::GetCellPower(float x, float y) const{
    pair<int, int> index = GetCellIndex(x,y);
    return GRID[index.first][index.second].power;  //[W]
}

int lightmap::GetCellNx() const{
    return n_cell_x;
}

int lightmap::GetCellNy() const{
    return n_cell_y;
}

const vector<vector<cell>>& lightmap::GetCells() const{
    return GRID;
}

array<float,3> lightmap::GetNormal(){
    return {0.,0.,1.};
}

void lightmap::AddLightSource(lightsource S){
    //Adicionar a luz ao vetor das luzes
    fotoes.push_back(S);
    
    //Calcular a potência [W] proveniente da luz para todas as células
    //Potência = irradiânica * área
    for (int i=0; i<GRID.size(); i++){
        for (int j=0; j<GRID[0].size(); j++){
            GRID[i][j].power+=this->Irradiance(GRID[i][j].center_coo)*GRID[i][j].area; //[W]
        }
    }
}


const cell& lightmap::GetMaxCell() const{
    //determinar a céula com potência máxima da grelha
    float max_power=0.; //inteiro para avaliar se existe uma potência maior
    pair<int, int> index;  //par de valores relativos ao índice da célula com maior potência
    for (int i=0; i<GRID.size();i++){
        for(int j=0; j<GRID[0].size();j++){
            if (GRID[i][j].power>max_power){
                max_power=GRID[i][j].power;
                index.first=i;
                index.second=j;
            }
        }
    }
    return GRID[index.first][index.second];
}

float lightmap::distance2cell(array<float,3> SourceCoo, array<float,3> PointCoo){
    return sqrt(pow((SourceCoo[0]-PointCoo[0]),2)+pow((SourceCoo[1]-PointCoo[1]),2)+pow((SourceCoo[2]-PointCoo[2]),2));
}

float lightmap::Irradiance(array<float,3> PointCoo){
    //vetor entre a luz e a celula
    array<float, 3> r = {PointCoo[0]-fotoes[fotoes.size()-1].coo[0], PointCoo[1]-fotoes[fotoes.size()-1].coo[1], PointCoo[2]-fotoes[fotoes.size()-1].coo[2]};
    
    //distanica entre a luz e a celula
    float distancia=distance2cell(fotoes[fotoes.size()-1].coo, PointCoo);

    float cos_alfa = r[2]/ distancia; //As coordenada x e y são 0
    float irrad = (fotoes[fotoes.size()-1].power/(4*M_PI)) * (cos_alfa/(distancia*distancia)); //[W/cm²]
    return irrad;
}

float lightmap::CellPower(const cell& C){
    return C.power;
}

#endif