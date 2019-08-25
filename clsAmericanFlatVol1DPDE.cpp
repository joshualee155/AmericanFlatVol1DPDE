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
vector<double> tri_solver(double, double, double, const vector<double> &);

clsAmericanFlatVol1DPDE::clsAmericanFlatVol1DPDE(double _S, double _K, double _r, double _T, double _sig, double _div, bool _isCall, int _nsteps, int _tsteps, double _theta, int _nStdev):
S0(_S), K(_K), T(_T), r(_r), vol(_sig), div(_div), isCall(_isCall), Ssteps(_nsteps), Tsteps(_tsteps), theta(_theta), nStdev(_nStdev)
{
    C.resize(Tsteps+1);
    for (auto i = 0; i < Tsteps+1; ++i )
    {
        C[i].resize(Ssteps+1);
    }
    isDone = false;
}

clsAmericanFlatVol1DPDE::~clsAmericanFlatVol1DPDE()
{
    for (auto i = 0; i < Tsteps+1; ++i )
    {
        C[i].resize(0);
    }
    C.resize(0);
}

void clsAmericanFlatVol1DPDE::run()
{
    auto const d_t = T/Tsteps;
    auto const d_x = 2.0*vol*nStdev/Ssteps;
    auto const beta = 0.5*d_t*vol*vol/d_x/d_x;
    auto const alpha = 0.25*d_t*vol*vol/d_x;
    
    vector<double> x(Ssteps+1), s(Ssteps+1);
    for (auto j = 0; j < Ssteps+1; ++j)
    {
        x[j] = - nStdev*vol + j*d_x;
        s[j] = S0*exp(x[j] + (r - div)*T);
    }
    
    for (auto j = 0; j < Ssteps+1; ++j)
    {
        C[Tsteps][j] = isCall ? max(s[j] - K, 0.0) : max(K - s[j], 0.0);
    }
    
    for (auto i = Tsteps - 1; i >= 0; --i)
    {
        // Fill in the boundary conditions
        C[i][0] = isCall ? 0.0 : K*exp(-r*(T - i*d_t)) - s[0]*exp(-(r - div)*(T-i*d_t));
        C[i][Ssteps] = isCall ? s[Ssteps]*exp(-(r - div)*(T-i*d_t)) - K*exp(-r*(T - i*d_t)) : 0.0;
        
        // Correct the sub-boundary points
        C[i+1][1] += (alpha+beta)*C[i][0];
        C[i+1][Ssteps-1] -= (alpha - beta)*C[i][Ssteps];
        
        auto tmp_sol = tri_solver(1.0 + 2.0*beta, -(alpha+beta), alpha-beta, C[i+1]);
        
        for ( auto j = 1 ; j < Ssteps; ++j)
        {
            C[i][j] = max(tmp_sol[j-1], isCall ? max( s[j]*exp(-(r - div)*(T-i*d_t)) - K, 0 ) : max( K - s[j]*exp(-(r - div)*(T-i*d_t)), 0 ));
        }
    }
    isDone = true;
}

double clsAmericanFlatVol1DPDE::get_price(const double s_t)
{
    if (!isDone)
    {
        std::cerr << "Call run() before query for price!" << std::endl;
        return -1;
    }
    auto const d_x = 2.0*vol*nStdev/Ssteps;
    vector<double> x(Ssteps), s(Ssteps);
    
    int j;
    for ( j = 0; j < Ssteps+1; ++j )
    {
        x[j] = - nStdev*vol + j*d_x;
        s[j] = S0*exp(x[j]);
        if (s_t < s[j]) break;
    }
    if (j == 0 || j == Ssteps)
    {
        return C[0][j];
    }
    else
    {
        return C[0][j-1] + (C[0][j]-C[0][j-1])/(s[j]-s[j-1])*(s_t-s[j-1]);
    }
}

double max(double a, double b)
{
    return (a>b) ? a : b;
}

vector<double> tri_solver(double b, double a, double c, const vector<double> & y_)
{
    auto const m = static_cast<int>(y_.size()-2);
    //vector<double> y(y_);
    auto y(y_);

    vector<double> sol(m), c_vec(m-1);
    c_vec[0] = c/b;
    y[1] /= b;
    for (auto i = 1; i < m - 1 ; ++i)
    {
        c_vec[i] = c/(b - a*c_vec[i-1]);
        y[i+1] = (y[i+1] - a*y[i])/(b - a*c_vec[i-1]);
    }
    y[m] = (y[m] - a*y[m-1])/(b - a*c_vec[m-2]);
    sol[m-1] = y[m];
    for (auto i = m-2; i >= 0; --i)
    {
        sol[i] = y[i+1]-c_vec[i]*sol[i+1];
    }
    return sol;
}

