#include "matrix.hpp"
#include <vector>

struct krylov_old
{
    int m;
    mat H_bar;
    mat H;
    std::vector<std::vector<double>> V;
};


struct krylov
{

    std::vector<std::vector<double>> H_bar;
    std::vector<std::vector<double>> V;
};

std::vector<double> getKrylov(mat A, std::vector<std::vector<double>> &V)
{
    std::cout << V[0][0] << std::endl;
    int sz = A.N;
    assert(sz == A.M);
    for (std::vector<double>v : V) assert(sz == v.size()); 

    //std::cout << "hallo" << std::endl;
    int j = V.size();

    std::vector<double> w(sz);
    w = vecprod(A,V[j-1]);
    
    //std::cout << "hallo 2" << std::endl;
    std::vector<double> h(j+1);


    for (size_t i = 0; i < j; i++)
    {
        h[i] = dot(V[i], w);

        //std::cout << i << " " << j << std::endl;
        for (size_t k = 0; k < sz; k++)
        {
            //if (k <= 10) printf("%d %d %.3f %.3f %.3f \n", (int) i, (int) k , w[k], h[i], V[i][k]);
            w[k] -= h[i]*V[i][k];
        }
    }
    h[j] = norm(w);
    //printf("% .3e \n", h[j]);
    std::vector<double> v(sz);

    for (size_t k = 0; k < sz; k++) v[k]=w[k]/h[j];

    V.push_back(v);
    
    return h;
}

krylov gramSchmidt(mat A, int m, std::vector<double> r0){

    krylov K;
    double r = norm(r0);
    int sz = r0.size();

    std::vector<double> v(sz);
    std::vector<std::vector<double>> H_bar;
    std::vector<std::vector<double>> V;

    for (size_t i = 0; i < sz; i++) v[i] = r0[i]/r;
    V.push_back(v);


    //std::vector<double> h(m+1);
    for (size_t i = 0; i < m; i++)
    {
        std::vector<double> h_temp(i+2);
        h_temp = getKrylov(A, V);
        H_bar.push_back(h_temp);
        //for (size_t j = 0; j <= j+1; j++) h[j]=h_temp;
    }


    for (size_t i = 0; i < m; i++)
    {
        printf("i= ");
        for (size_t j = 0; j < m; j++)
        {
            printf("% .3e ",dot(V[i],V[j]));
        }
        printf("\n");
    }

    K.H_bar  = H_bar;
    K.V = V;
    return K;
}

void gmres(mat A, std::vector<double> x0, std::vector<double> b, int m)
{
    int sz = A.N;

    //impelemnt checks for dimenstions

    std::vector<double> v(sz);
    std::vector<double> r0(sz);
    std::vector<double> g;
    r0 = vecprod(A,b);
    double r = norm(r0);
    for (size_t i = 0; i < sz; i++) r0[i] = b[i]-r0[i];
    for (size_t i = 0; i < sz; i++) v[i] = r0[i]/r;

    g.resize(m+1,0);
    g[0]=r;

    krylov K;
    std::vector<double> h;
    std::vector<std::vector<double>> H_bar;
    std::vector<std::vector<double>> V;
    V.push_back(v);

    std::vector<double> c(m);
    std::vector<double> s(m);
    for (size_t j = 0; j < m; j++)
    {
        printf("m  %d \n", (int) j);
        h.resize(j+2,0.0);
        h = getKrylov(A, V);
        H_bar.push_back(h);
        for (size_t k = 1; k <= j; k++)
        {
            double temp_h  =  c[k-1]*H_bar[j][k-1] + s[k-1]*H_bar[j][k];
            H_bar[j][k] = -s[k-1]*H_bar[j][k-1] + c[k-1]*H_bar[j][k];
            H_bar[j][k-1] = temp_h;
        }
        //printf("Hallo 2\n");

        double temp_sc = sqrt(pow(H_bar[j][j], 2)+pow(H_bar[j][j+1], 2));
        c[j] = H_bar[j][j]  /temp_sc;
        s[j] = H_bar[j][j+1]/temp_sc;
        //printf("Hallo 3\n");
        H_bar[j][j] = c[j]*H_bar[j][j] + s[j]*H_bar[j][j+1];
        g[j+1]= -s[j]*g[j];
        printf("err % .3e \n",abs(g[j+1])/r);
        g[j] = c[j]*g[j];
    }


}