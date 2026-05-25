#include <Eigen/Dense>
#include <iostream>
#include "FCmatrixAlgo.h"
using namespace std;
using namespace Eigen;

#ifndef __EqSolver_h__
#define __EqSolver_h__

class EqSolver {
    public:
        // constructors and destructor
        EqSolver();
        EqSolver(const Matrix<double,Dynamic,Dynamic>&, const Matrix<double,Dynamic,1>&);
        ~EqSolver() = default;

        // output (optional)
        //friend ostream& operator<<(ostream&, const EqSolver&);

        // solvers
        const Matrix<double,Dynamic,1>& GaussSolver(bool pivot=false);
        const Matrix<double,Dynamic,1>& LUSolver(bool pivot=false);
        void IterativeJacobiSolver(Matrix<double, Dynamic,1>&, int& itmax, double tol=1.E-3); // tolerance on convergence
        void IterativeGaussSeidelSolver(Matrix<double,Dynamic,1>&,int& itmax, double tol=1.E-3);
        
    private:
        Matrix<double,Dynamic,Dynamic> M; //coefficients matrix
        Matrix<double,Dynamic,1> b; // constants vector

        Matrix<double,Dynamic,Dynamic> M_copia; //copy coefficients matrix
        Matrix<double,Dynamic,1> b_copia; // copy constants vector
        //FCmatrixAlgo F;
};

#endif