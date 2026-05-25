#include <iostream>
#include <vector>
#include <array>
using namespace std;

#ifndef mapas_de_luz_h
#define mapas_de_luz_h


struct cell {
    array<float,3> center_coo {0,0,100}; // cell center coordinates (cm)↪→
    float area=0; // cell area (cmˆ2)
    float power=0; // cell incident power (W)
};

struct lightsource {
    array<float,3> coo{0.,0.,0.}; // source coordinates (cm)
    float power=100; // source power (W)
};

class lightmap
{
private:
    vector<vector<cell> > GRID; // cell grid
    int n_cell_x;
    int n_cell_y;
    float tamanho_x;
    float tamanho_y;
    float dx; //tamanho de cada celula em x
    float dy; //tamanho de cada celula em y
    vector<lightsource> fotoes;

public:
    lightmap() = default;
    lightmap(array<int,2> ncell, array<float,2> size);
    lightmap(lightsource S, array<int,2> ncell, array<float,2> size);
    pair<int,int> GetCellIndex(float x, float y) const; // return cell indices
    pair<float,float> GetCellCoo(int index_x, int index_y) const; //return cell coordinates
    float GetCellPower(int index_x, int index_y) const; // return cell power Watts
    float GetCellPower(float x, float y) const;
    int GetCellNx() const; // get number of cells along x
    int GetCellNy() const;

    const vector<vector<cell>>& GetCells() const; // return cells grid
    array<float,3> GetNormal(); // return normal to plane

    void AddLightSource(lightsource S); // add sources
    const cell& GetMaxCell() const; // get cell with max power
    float distance2cell(array<float,3> SourceCoo, array<float,3> PointCoo);
    float Irradiance(array<float,3> PointCoo);
    float CellPower(const cell& C);
};

#endif

