#include <iostream>
#include <Eigen/Core>
using namespace std;
using namespace Eigen;

#ifndef FCmatrixAlgo_h
#define FCmatrixAlgo_h

class FCmatrixAlgo {
    public:
    FCmatrixAlgo() = default;
    ~FCmatrixAlgo() = default;

    /*
    Implements Gauss elimination
    */
    
    void GaussElimination(Matrix<double, Dynamic, Dynamic>&, Matrix<double, Dynamic,1>& ); //no pivoting

    void GaussEliminationPivot(Matrix<double,Dynamic, Dynamic>&, Matrix<double, Dynamic,1>&, Matrix<double,Dynamic,1>&); //make pivoting
    
    /*
    Implements LU decomposition (Doolitle)
    */
    //CASO CORRA MAL
    //void LUdecomposition(Matrix<double, Dynamic, Dynamic>&, Matrix<int,Dynamic,1>&, bool bpivot=false);

    void LUdecomposition(Matrix<double, Dynamic, Dynamic>&, Matrix<double,Dynamic,1>&, bool bpivot=false);

};

#endif