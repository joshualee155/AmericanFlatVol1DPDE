//
//  clsAmericanFlatVol1DPDE.cpp
//  AmericanFlatVol1DPDE
//
//  Created by Xiangyu Joshua Li on 19/09/2015.
//  Copyright © 2015 Xiangyu Joshua Li. All rights reserved.
//

#include "clsAmericanFlatVol1DPDE.hpp"
#include <iostream>
#include <cmath>

double max(double, double);
vector<double> triSolver(double, double, double, vector<double>);

clsAmericanFlatVol1DPDE::clsAmericanFlatVol1DPDE(double _S, double _K, double _r, double _T, double _sig, double _div, bool _isCall, int _nsteps, int _tsteps, double _theta, int _nStdev):
S0(_S), K(_K), T(_T), r(_r), vol(_sig), div(_div), isCall(_isCall), Ssteps(_nsteps), Tsteps(_tsteps), theta(_theta), nStdev(_nStdev)
{
    C.resize(Tsteps+1);
    for (int i = 0; i < Tsteps+1; ++i )
    {
        C[i].resize(Ssteps+1);
    }
    isDone = false;
}

clsAmericanFlatVol1DPDE::~clsAmericanFlatVol1DPDE()
{}

void clsAmericanFlatVol1DPDE::run()
{
    double d_t = T/Tsteps;
    double d_x = 2.0*vol*nStdev/Ssteps;
    double beta = 0.5*d_t*vol*vol/d_x/d_x;
    double alpha = 0.25*d_t*vol*vol/d_x;
    
    vector<double> x(Ssteps), S(Ssteps);
    for ( int j = 0; j < Ssteps+1; ++j )
    {
        x[j] = - nStdev*vol + j*d_x;
        S[j] = S0*exp(x[j] + (r - div)*T);
    }
    
    for (int j = 0; j < Ssteps+1; ++j)
    {
        C[Tsteps][j] = isCall ? max(S[j] - K, 0.0) : max(K - S[j], 0.0);
    }
    
    for (int i = Tsteps - 1; i >= 0; --i)
    {
        // Fill in the boundary conditions
        C[i][0] = isCall ? 0.0 : K*exp(-r*(T - i*d_t)) - S[0]*exp(-(r - div)*(T-i*d_t));
        C[i][Ssteps] = isCall ? S[Ssteps]*exp(-(r - div)*(T-i*d_t)) - K*exp(-r*(T - i*d_t)) : 0.0;
        
        // Correct the sub-boundary points
        C[i+1][1] += (alpha+beta)*C[i][0];
        C[i+1][Ssteps-1] -= (alpha - beta)*C[i][Ssteps];
        
        vector<double> tmp_Sol = triSolver(1.0 + 2.0*beta, -(alpha+beta), alpha-beta, C[i+1]);
        
        for ( int j = 1 ; j < Ssteps; ++j)
        {
            C[i][j] = max(tmp_Sol[j-1], isCall ? max( S[j]*exp(-(r - div)*(T-i*d_t)) - K, 0 ) : max( K - S[j]*exp(-(r - div)*(T-i*d_t)), 0 ));
        }
    }
    isDone = true;
}

double clsAmericanFlatVol1DPDE::get_price(double S_t)
{
    if (!isDone)
    {
        std::cerr << "Call run() before query for price!" << std::endl;
        return -1;
    }
    double d_x = 2.0*vol*nStdev/Ssteps;
    vector<double> x(Ssteps), S(Ssteps);
    
    int j;
    for ( j = 0; j < Ssteps+1; ++j )
    {
        x[j] = - nStdev*vol + j*d_x;
        S[j] = S0*exp(x[j]);
        if (S_t < S[j]) break;
    }
    if (j == 0 || j == Ssteps)
    {
        return C[0][j];
    }
    else
    {
        return C[0][j-1] + (C[0][j]-C[0][j-1])/(S[j]-S[j-1])*(S_t-S[j-1]);
    }
}

double max(double a, double b)
{
    return (a>b) ? a : b;
}

vector<double> triSolver(double b, double a, double c, vector<double> y)
{
    int M = static_cast<int>(y.size()-2);
    vector<double> Sol(M), c_vec(M-1);
    c_vec[0] = c/b;
    y[1] /= b;
    for (int i = 1; i < M - 1 ; ++i)
    {
        c_vec[i] = c/(b - a*c_vec[i-1]);
        y[i+1] = (y[i+1] - a*y[i])/(b - a*c_vec[i-1]);
    }
    y[M] = (y[M] - a*y[M-1])/(b - a*c_vec[M-2]);
    Sol[M] = y[M];
    for (int i = M-1; i >= 0; --i)
    {
        Sol[i] = y[i+1]-c_vec[i]*Sol[i+1];
    }
    return Sol;
}

