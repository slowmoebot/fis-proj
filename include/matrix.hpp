#pragma once
#include <iostream>
#include <string>
#include <cassert>
#include <fstream>
#include <vector>
#include <math.h>

struct mat
{
    std::string mat_type;
    int N;
    int M;
    int nnz;
    std::vector<int> I;
    std::vector<int> J;
    std::vector<double> val;
};

struct mat_coo
{
    int N;
    int M;
    std::vector<int> I;
    std::vector<int> J;
    std::vector<double> val;
};

struct entry
{
    int row;
    int col;
    double val;
};


mat read_mat(std::string path){
    mat matrix;
    std::ifstream fin(path);
    fin >> matrix.mat_type;
    fin >> matrix.N >> matrix.nnz;

    matrix.M = matrix.N;

    matrix.I.resize(matrix.N+1, 0);
    for (size_t i = 0; i < matrix.N + 1; i++)
    {
        fin >> matrix.I[i];
        matrix.I[i]--;
    }

    matrix.J.resize(matrix.nnz, 0);
    matrix.val.resize(matrix.nnz, 0);
    for (size_t i = 0; i < matrix.nnz; i++)
    {
        fin >> matrix.J[i] >> matrix.val[i];
        matrix.J[i]--;
    }

    fin.close();
    return matrix;
}

mat coo2mat(mat_coo inp){
    mat matrix;
    matrix.mat_type = "n";
    matrix.N = inp.N;
    matrix.M = inp.M;
    matrix.nnz = inp.I.size();
    matrix.I.resize(matrix.N+1,0);
    matrix.J.resize(matrix.nnz,0);
    matrix.val.resize(matrix.nnz,0.0);

    std::vector<std::vector<entry>> temp_mat(matrix.N);

    for (size_t i = 0; i < matrix.nnz; i++)
    {
        entry temp_entry;
        temp_entry.row = inp.I[i];
        temp_entry.col = inp.J[i];
        temp_entry.val = inp.val[i];
        temp_mat[temp_entry.row].push_back(temp_entry);
    }
    
    int k = 0;
    for (size_t i = 0; i < matrix.N; i++)
    {
        matrix.I[i]=k;
        for (entry e : temp_mat[i])
        {
            matrix.J[k] = e.col;
            matrix.val[k] = e.val;
            k++;
        }    
    }
    matrix.I[matrix.N] = matrix.nnz;

    return matrix;
    
}

mat coo2mat_old(mat_coo inp){


    int nnz = inp.I.size();

    mat matrix;

    matrix.N = inp.N;
    matrix.M = inp.M;
    matrix.mat_type = "n";
    matrix.nnz = nnz;
    matrix.I.resize(inp.M+1,0);
    matrix.J.resize(nnz,0);
    matrix.val.resize(nnz,0.0);


    std::vector<entry> temp_mat[inp.M];

    assert(nnz==inp.J.size() && nnz==inp.val.size());

    for (int i = 0; i < nnz; i++)
    {
        entry temp_entry;
        temp_entry.row = inp.I[i];
        temp_entry.col = inp.J[i];
        temp_entry.val = inp.val[i];
        temp_mat[temp_entry.row].push_back(temp_entry);
    }

    int k=0;
    for (int i=0; i<inp.M; i++)
    {
        std::cout << i << " " << matrix.I.size() << std::endl;

        matrix.I[i]=k;
        for(entry e: temp_mat[i])
        {
            matrix.J[k]=e.col;
            matrix.val[k]=e.val;
            k++;
        }
    }
    matrix.I[inp.M] = nnz;
    return matrix;
}

double dot(std::vector<double> a, std::vector<double> b){
    int sz = a.size();
    assert(b.size()==sz);
    double res = 0;
    for (size_t i = 0; i < sz; i++)
    {
        res += a[i]*b[i];
    }
    return res;
}

double norm(std::vector<double> x){
    int sz = x.size();
    double res = 0;
    for (size_t i = 0; i < sz; i++)
    {
        res += pow(x[i],2.0);
    }
    return sqrt(res);
}

std::vector<double> vecprod_csr(mat matrix, std::vector<double> vec){
    int sz = vec.size();
    assert(matrix.N == sz);

    std::vector<double> res;
    res.resize(matrix.M,0);

    for (size_t i = 0; i < matrix.M; i++)
    {
        int i1 = matrix.I[i];
        int i2 = matrix.I[i+1];
        int temp_len = i2 - i1; 
        std::vector<double> temp1(temp_len);
        std::vector<double> temp2(temp_len);
        for (size_t j = i1; j < i2; j++)
        {
            temp1[j - i1] = matrix.val[j];
            temp2[j - i1] = vec[matrix.J[j]];
        }
        double temp_res = dot(temp1, temp2);
        res[i] = temp_res;
    }

    return res;
}

std::vector<double> vecprod_csc_t(mat matrix, std::vector<double> vec)
{
    size_t sz=vec.size();

    assert(sz==matrix.N);

    std::vector<double> res;
    res.resize(matrix.M,0);

    for (size_t i = 0; i < matrix.M; i++)
    {
        int i1 = matrix.I[i];
        int i2 = matrix.I[i+1];
        for (size_t j = i1; j < i2; j++)
        {
            if (i != matrix.J[j]) res[matrix.J[j]] += matrix.val[j]*vec[i];
        }
    }
    return res;
}

std::vector<double> vecprod(mat matrix, std::vector<double> vec)
{
    // Get result from vector product in compressed sparse row format
    
    std::vector<double> res = vecprod_csr(matrix, vec);

    // If matrix is symmetric the transpose of the matrix needs to be multiplied as well
    if (matrix.mat_type == "s") {

        // Get result from transposed vector product without diagonal
        std::vector<double> temp = vecprod_csc_t(matrix, vec);


        // Add results
        for (size_t i = 0; i < matrix.M; i++) res[i] += temp[i];     
    }

    //Return final result
    return res;
}