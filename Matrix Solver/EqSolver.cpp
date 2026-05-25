#include <iostream>
#include "EqSolver.h"
 //DEPOIS METER O .h
using namespace Eigen;


EqSolver::EqSolver(){}

EqSolver::EqSolver(const Matrix<double, Dynamic, Dynamic>& matriz, const Matrix<double, Dynamic, 1>& const_vector) {
    M = matriz;
    b = const_vector;

    //Guarda uma cópia da matriz e do vetor para poder alterar os outros, mantendo a matriz e o vetor originais
    M_copia = M;
    b_copia = b;
}

const Matrix<double, Dynamic, 1>& EqSolver::GaussSolver(bool pivot) {
    FCmatrixAlgo F;

    //Iguala a matriz e o vetor aos originais
    M = M_copia;
    b = b_copia;

    Matrix<double, Dynamic, 1> solucao(M.rows());
    solucao.setZero();

    if (!pivot) {
        F.GaussElimination(M, b);     
    } 

    else {
        Matrix<double, Dynamic, 1> r_i(M.rows());
        for (int i = 0; i < M.rows(); i++) {
            r_i(i) = i;
        }
        F.GaussEliminationPivot(M, b, r_i);
    }

    // Solve
        for (int k=M.rows()-1;k>=0;k--){
            double somatorio=0.;
            for (int j=k+1;j<M.cols();j++){
                somatorio+=M(k,j)*solucao(j);
            }
            solucao(k)=(b(k)-somatorio)/M(k,k);
        }
    
    b = solucao;
    cout << "A solucao através de gauss:\n" << b << endl;
    return b;
}

//DEPOIS VER SE NÃO VOLTO ATRÁS
/*
const Matrix<double,Dynamic,1>& EqSolver::LUSolver(bool pivot){
    FCmatrixAlgo F;
    //Criação da matriz para a ordem das linhas
    Matrix<int, Dynamic, 1> r_i(M.rows());
    for (int i = 0; i < M.rows(); i++) {
        r_i(i) = i;
    }

    //Sem pivot
    if (!pivot) {
        F.LUdecomposition(M, r_i);   
    } 

    //Com pivot
    else {
        F.LUdecomposition(M,r_i,true);
    }
    //trocar o valor das constantes
    for (int ordem = 0; ordem < r_i.rows(); ordem++){
    if (r_i(ordem) > ordem) {
        int temp = b(ordem);
        b(ordem) = b(r_i(ordem));
        b(r_i(ordem)) = temp;
        }
    }


    //Resolver, ir descobrir a primeira matriz(y)
    Matrix<double, Dynamic, 1> matriz_y(M.rows());
    matriz_y.setZero();
    for(int k=0;k<matriz_y.rows();k++){
        double somatorio=0.;
        for(int j=0;j<k;j++){
            somatorio+=M(k,j)*matriz_y(j);
        }
        matriz_y(k)=b(k)-somatorio;
    }

    //Resolucao 2 matriz(x)
    Matrix<double, Dynamic, 1> solucao(M.rows());
    solucao.setZero();
    for(int k=solucao.rows()-1;k>=0;k--){
        double somatorio=0.;
        for(int j=k+1;j<solucao.rows();j++){
            somatorio+= M(k,j)*solucao(j);
        }
        solucao(k)=(matriz_y(k)-somatorio)/M(k,k);
    }
    cout << "a solucao do LU SOLVER é:\n"<< solucao<< endl;
    b = solucao;
    return b;
}*/

const Matrix<double,Dynamic,1>& EqSolver::LUSolver(bool pivot){
    FCmatrixAlgo F;

    //Iguala a matriz e o vetor aos originais
    M = M_copia;
    b = b_copia;


    //Sem pivot
    if (!pivot) {
        F.LUdecomposition(M, b);   
    } 

    //Com pivot
    else {
        F.LUdecomposition(M,b,true);
    }


    //Resolver, ir descobrir a primeira matriz(y)
    Matrix<double, Dynamic, 1> matriz_y(M.rows());
    matriz_y.setZero();
    for(int k=0;k<matriz_y.rows();k++){
        double somatorio=0.;
        for(int j=0;j<k;j++){
            somatorio+=M(k,j)*matriz_y(j);
        }
        matriz_y(k)=b(k)-somatorio;
    }

    //Resolucao 2 matriz(x)
    Matrix<double, Dynamic, 1> solucao(M.rows());
    solucao.setZero();
    for(int k=solucao.rows()-1;k>=0;k--){
        double somatorio=0.;
        for(int j=k+1;j<solucao.rows();j++){
            somatorio+= M(k,j)*solucao(j);
        }
        solucao(k)=(matriz_y(k)-somatorio)/M(k,k);
    }
    cout << "a solucao do LU SOLVER é:\n"<< solucao<< endl;
    b = solucao;
    return b;
}




void EqSolver::IterativeJacobiSolver(Matrix<double, Dynamic,1>& x, int& itmax, double tol){
    //Iguala a matriz e o vetor aos originais
    M = M_copia;
    b = b_copia;

    
    // Verificar se a matriz é diagonalmente dominante, caso não o seja gera um aviso de que a solução pode não ser convergente
    // Para isso, verifica-se se a entrada da diagonal de cada linha é maior do que a soma das restantes entradas dessa linha
    for (int row = 0; row<M.rows(); row++)
    {
        double soma = 0.;
        for (int col = 0; col<M.cols(); col++)
        {
            if (row!=col) soma += fabs(M(row, col));
        }

        if (fabs(M(row, row))>=soma) continue;
        else
        {
            cout << "Matrix não é diagonalmente dominante, logo pode não convergir" << endl;
        }
    }


    /*VectorXd x_new = x;
    bool end = false;
    double error;
    for (int it = 1; it<=itmax; it++)
    {
        for (int i = 0; i < M.rows(); i++)
        {
            double soma = 0.;
            for (int j = 0; j < M.cols(); j++)
            {
                if (j != i) soma += M(i, j) * x(j);
            }
            x_new(i) = (b(i) - soma) /M(i, i);
        }
        error = (x_new -x).norm();

        if (error<tol){
            printf("Converged in %i iterations.\n", it);
            end = true;
            break;
        }
        x = x_new;
    }
    if (!end) cout << "Did not converge within the maximum number of iterations." << endl;
    cout << " A solucão é do jacobi:\n " << x<<endl;*/



        // linear system of m unknowns
        VectorXd x_aux=VectorXd::Zero(M.rows()); //zero's
        bool btol = false;
        int it = 0.;
        double eps = 1.E-4; //tolerance
        while (!btol && (it++ < itmax)) {
            x_aux = x;
            for (int i=0; i<M.rows(); i++) {
            x(i) = 0.;
            for (int j=0; j<M.rows(); j++){
                if (i != j) x(i) += -M(i,j)*x_aux(j);
            }
            x(i) += b(i);
            x(i) /= M(i,i);

            //guarantee that all vector entries are converging equally
            for(int i=0; i<M.rows(); i++){
                if (abs(x(i)-x_aux(i))<eps){
                    btol=true;
                }
                else{
                    btol=false;
                    break;
                }
            }
            }

        }
        cout << " A solucão do jacobi é:\n " << x<<endl;
    }
    

void EqSolver::IterativeGaussSeidelSolver(Matrix<double,Dynamic,1>& x,int& itmax, double tol){

    //Iguala a matriz e o vetor aos originais
    M = M_copia;
    b = b_copia;
    
    // Verificar se a matriz é diagonalmente dominante, caso não o seja gera um aviso de que a solução pode não ser convergente
    // Para isso, verifica-se se a entrada da diagonal de cada linha é maior do que a soma das restantes entradas dessa linha
    for (int row = 0; row<M.rows(); row++)
    {
        double soma = 0.;
        for (int col = 0; col<M.cols(); col++)
        {
            if (row!=col) soma += fabs(M(row, col));
        }

        if (fabs(M(row, row))>=soma) continue;
        else
        {
            cout << "Matrix não é diagonalmente dominante, logo pode não convergir" << endl;
        }
    }

    VectorXd x_aux=VectorXd::Zero(M.rows()); //zero's
    bool btol = false;
    int it = 0.;
    double eps = 1.E-4; //tolerance
    while (!btol && (it++ < itmax)) {
        x_aux = x;
        for (int i=0; i<M.rows(); i++) {
            x(i) = 0.;
            for (int j=0; j<M.rows(); j++){
            if (i != j) x(i) += -M(i,j)*x(j);
            }
            x(i) += b(i);
            x(i) /= M(i,i);
            //guarantee that all vector entries are converging equally
            /*if (fabs(x(i)-x_aux(i)) < eps){
                btol = true;
                }
            else btol = false;*/
            for(int i=0; i<M.rows(); i++){
                if (abs(x(i)-x_aux(i))<eps){
                    btol=true;
                }
                else{
                    btol=false;
                    break;
                }
            }
        }
    }
    cout << " A solucão do seidel é:\n " << x<<endl;
}
    