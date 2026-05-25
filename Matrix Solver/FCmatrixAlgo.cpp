#include <iostream>
#include "FCmatrixAlgo.h"
using namespace std;
using namespace Eigen;

void FCmatrixAlgo::GaussElimination(Matrix<double, Dynamic, Dynamic>& coeficientes, Matrix<double, Dynamic, 1>& const_vector) {
    int cols = coeficientes.cols();
    int rows = coeficientes.rows();
    int linhas = 0;
    int colunas = 0;

    while (linhas < rows && colunas<cols) {
        //cout << coeficientes(linhas, colunas) << endl;

        if (coeficientes(linhas, colunas) == 0) {
            bool troca = false;

            for (int k = linhas + 1; k < rows; k++) {
                if (coeficientes(k, colunas) != 0) {
                    coeficientes.row(linhas).swap(coeficientes.row(k));
                    swap(const_vector(linhas), const_vector(k));
                    troca = true;
                    break;
                }
            }

            if (!troca) {
                // Todos os elementos são zero nesta coluna, avance para a próxima coluna
                colunas++;
                continue;
            }
        }

        //cout << "MATRIX:\n" << coeficientes << endl;

        if (coeficientes(linhas, colunas) != 0) {
            for (int j = linhas + 1; j < rows; j++) {
                double lambda = coeficientes(j, colunas) / coeficientes(linhas, colunas);
                //cout << "LAMBDA: " << lambda << endl;
                coeficientes.row(j) -= lambda * coeficientes.row(linhas);
                const_vector(j) -= lambda * const_vector(linhas);
                //cout << coeficientes << endl;
            }

            linhas++;
        }

        colunas++;
    }

    //cout << "A matriz é:\n" << coeficientes << endl;
    //cout << "vetor solucoes\n" << const_vector << endl;
}


void FCmatrixAlgo::GaussEliminationPivot(Matrix<double,Dynamic, Dynamic>& coeficientes, Matrix<double, Dynamic,1>& const_vector, Matrix<double,Dynamic,1>& row_indexing){
    int cols = coeficientes.cols();  //numero de colunas
    int rows = coeficientes.rows();  //numero de linhas
    int linhas = 0;  //iterar sobre as linhas
    int colunas = 0;  //iterar sobre as colunas

    //Vetor que dá o valor máximo de cada linha
    VectorXd escala(rows);
    //escala=coeficientes.rowwise().maxCoeff();
    for(int i=0; i<rows; i++){
        int maximo =0;
        for (int j=0; j<cols; j++){
            if(maximo<abs(coeficientes(i,j))){
                maximo = abs(coeficientes(i,j));
            }
        }
        escala(i)=maximo;
    }
    

    while (linhas < rows && colunas<cols){
        
        //Vetor que dá o valor relativo de cada linha
        VectorXd ordem(rows);
        for (int i = 0; i < coeficientes.rows(); i++) {
            ordem(i) = abs(coeficientes(i, colunas) / escala(i));
        }

        //iteração para descobrir a linha com maior valor relativo
        int pos=linhas;
        int max=ordem(linhas);
        for (int k=linhas; k<rows; k++){
            if (ordem(k)>ordem(linhas)){
                pos=k;
                max=ordem(k);
            }
        }

        //trocar a ordem das linhas, se necessário
        if(pos!=linhas){
            coeficientes.row(linhas).swap(coeficientes.row(pos));
            swap(ordem(linhas), ordem(pos));
            swap(const_vector(linhas), const_vector(pos));
            swap(row_indexing(linhas), row_indexing(pos));

        }

        //se todos os valores da coluna do pivot forem 0, avanca só a coluna
        if (coeficientes(linhas, colunas) != 0) {
            for (int j = linhas + 1; j < rows; j++) {
                double lambda = coeficientes(j, colunas) / coeficientes(linhas, colunas);
                coeficientes.row(j) -= lambda * coeficientes.row(linhas);
                const_vector(j) -= lambda * const_vector(linhas);
                
            }
            linhas++;
        }
        colunas++;
    }
}

//SE TIVER DE VOLTAR ATŔAS
/*void FCmatrixAlgo::LUdecomposition(Matrix<double, Dynamic, Dynamic>& coeficientes, Matrix<int,Dynamic,1>& row_indexing, bool bpivot){
    int cols = coeficientes.cols();
    int rows = coeficientes.rows();
    int linhas = 0;
    int colunas = 0;
    MatrixXd L_M(rows,cols);
    L_M.setZero();

    //SEM PIVOTING
    if(bpivot==false){
        while (linhas < rows && colunas<cols) {
            //cout << coeficientes(linhas, colunas) << endl;

            if (coeficientes(linhas, colunas) == 0) {

                for (int k = linhas + 1; k < rows; k++) {
                    if (coeficientes(k, colunas) != 0) {
                        coeficientes.row(linhas).swap(coeficientes.row(k));
                        swap(row_indexing(linhas), row_indexing(k));
                        break;
                    }
                }
            }


            //cout << "MATRIX:\n" << coeficientes << endl;

            if (coeficientes(linhas, colunas) != 0) {
                for (int j = linhas + 1; j < rows; j++) {
                    double lambda = coeficientes(j, colunas) / coeficientes(linhas, colunas);
                    //cout << "LAMBDA: " << lambda << endl;
                    coeficientes.row(j) -= lambda * coeficientes.row(linhas);
                    L_M(j,colunas)= lambda;
                }

                linhas++;
            }
            colunas++;
        }
    }

    //COM PIVOTING
    else{
        //Vetor que dá o valor máximo de cada linha
        VectorXd escala(rows);
        escala=coeficientes.rowwise().maxCoeff();

        while (linhas < rows && colunas<cols){
            //Vetor que dá o valor relativo de cada linha
            VectorXd ordem(rows);
            for (int i = 0; i < coeficientes.rows(); i++) {
                ordem(i) = abs(coeficientes(i, colunas) / escala(i));
            }

            //iteração para descobrir a linha com maior valor relativo
            int pos=linhas;
            int max=ordem(linhas);
            for (int k=linhas; k<rows; k++){
                if (ordem(k)>ordem(linhas)){
                    pos=k;
                    max=ordem(k);
                }
            }

            //trocar a ordem das linhas, se necessário
            if(pos!=linhas){
                coeficientes.row(linhas).swap(coeficientes.row(pos));
                swap(ordem(linhas), ordem(pos));
                swap(row_indexing(linhas), row_indexing(pos));
                swap(L_M(linhas), L_M(pos));
            }

            //se todos os valores da coluna do pivot forem 0, avanca só a coluna
            if (coeficientes(linhas, colunas) != 0) {
                for (int j = linhas + 1; j < rows; j++) {
                    double lambda = coeficientes(j, colunas) / coeficientes(linhas, colunas);
                    coeficientes.row(j) -= lambda * coeficientes.row(linhas);
                    L_M(j,colunas)=lambda;
                }
                linhas++;
            }
            colunas++;
        }
    }
    coeficientes+=L_M;
}
*/




void FCmatrixAlgo::LUdecomposition(Matrix<double, Dynamic, Dynamic>& coeficientes, Matrix<double,Dynamic,1>& b, bool bpivot){
    int cols = coeficientes.cols();
    int rows = coeficientes.rows();
    int linhas = 0;
    int colunas = 0;
    MatrixXd L_M(rows,cols);
    L_M.setZero();

    //SEM PIVOTING
    if(bpivot==false){
        while (linhas < rows && colunas<cols) {
            //cout << coeficientes(linhas, colunas) << endl;

            if (coeficientes(linhas, colunas) == 0) {

                for (int k = linhas + 1; k < rows; k++) {
                    if (coeficientes(k, colunas) != 0) {
                        coeficientes.row(linhas).swap(coeficientes.row(k));
                        swap(b(linhas), b(k));
                        break;
                    }
                }
            }


            //cout << "MATRIX:\n" << coeficientes << endl;

            if (coeficientes(linhas, colunas) != 0) {
                for (int j = linhas + 1; j < rows; j++) {
                    double lambda = coeficientes(j, colunas) / coeficientes(linhas, colunas);
                    //cout << "LAMBDA: " << lambda << endl;
                    coeficientes.row(j) -= lambda * coeficientes.row(linhas);
                    L_M(j,colunas)= lambda;
                }

                linhas++;
            }
            colunas++;
        }
    }

    //COM PIVOTING
    else{
        //Vetor que dá o valor máximo de cada linha
        VectorXd escala(rows);
        escala=coeficientes.rowwise().maxCoeff();

        while (linhas < rows && colunas<cols){
            //Vetor que dá o valor relativo de cada linha
            VectorXd ordem(rows);
            for (int i = 0; i < coeficientes.rows(); i++) {
                ordem(i) = abs(coeficientes(i, colunas) / escala(i));
            }

            //iteração para descobrir a linha com maior valor relativo
            int pos=linhas;
            int max=ordem(linhas);
            for (int k=linhas; k<rows; k++){
                if (ordem(k)>ordem(linhas)){
                    pos=k;
                    max=ordem(k);
                }
            }

            //trocar a ordem das linhas, se necessário
            if(pos!=linhas){
                coeficientes.row(linhas).swap(coeficientes.row(pos));
                swap(ordem(linhas), ordem(pos));
                swap(b(linhas), b(pos));
                swap(L_M(linhas), L_M(pos));
            }

            //se todos os valores da coluna do pivot forem 0, avanca só a coluna
            if (coeficientes(linhas, colunas) != 0) {
                for (int j = linhas + 1; j < rows; j++) {
                    double lambda = coeficientes(j, colunas) / coeficientes(linhas, colunas);
                    coeficientes.row(j) -= lambda * coeficientes.row(linhas);
                    L_M(j,colunas)=lambda;
                }
                linhas++;
            }
            colunas++;
        }
    }
    coeficientes+=L_M;
}