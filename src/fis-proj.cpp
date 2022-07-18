#include <iostream>
#include "matrix.hpp"
#include "gmres.hpp"
#include <random>

int main(int argc, char* argv[]){

    assert(argc>1 && "Path has to be passed as first argument");

    std::string path = argv[1]; 

    mat matrix = read_mat(path);
    

    std::vector<double> v;
    v.resize(matrix.M,0);
    v[0] = 1.0;

    std::cout << "Matrix Size " << matrix.N << ", number of non-zero elements " << matrix.nnz << std::endl;
    std::vector<double> res = vecprod(matrix, v);

    std::default_random_engine generator;   
    
    std::uniform_real_distribution<double> unif(0,1);
    double r = unif(generator);

    std::vector<double> b(matrix.N);
    std::vector<double> x0(matrix.M);
    std::vector<double> r0(matrix.M);

    for (size_t i = 0; i < matrix.M; i++) x0[i] = unif(generator);
    for (size_t i = 0; i < matrix.N; i++) b[i] = unif(generator);

    r0 = vecprod(matrix, x0);


    for (size_t i = 0; i < matrix.N; i++) r0[i] = b[i]-r0[i];

    int m = 5;
    krylov K = gramSchmidt(matrix, m, r0);

    printf("H_bar \n");
    for (size_t i = 0; i < m; i++)
    {
        printf("i= ");
        for (size_t j = 0; j <= i+1; j++)
        {
            printf("% .3e ",K.H_bar[i][j]);
        }
        printf("\n");
    }

    gmres(matrix,x0,b,200);


    return 0;
}