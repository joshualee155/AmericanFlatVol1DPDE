//
//  clsAmericanFlatVol1DPDE.hpp
//  AmericanFlatVol1DPDE
//
//  Created by Xiangyu Joshua Li on 19/09/2015.
//  Copyright © 2015 Xiangyu Joshua Li. All rights reserved.
//

#ifndef clsAmericanFlatVol1DPDE_hpp
#define clsAmericanFlatVol1DPDE_hpp

#include <vector>

using std::vector;

class clsAmericanFlatVol1DPDE
{
public:
    clsAmericanFlatVol1DPDE();
    clsAmericanFlatVol1DPDE(double, double, double, double, double, double, bool,  int = 50, int = 20, double = 0.0, int = 4);
    ~clsAmericanFlatVol1DPDE();
    
    void run();
    double get_price(double);
    
private:
    
    double S0, K, T, r, vol, div;
    bool isCall;
    double theta;
    int Tsteps;
    int Ssteps;
    int nStdev;
    bool isDone;
    vector <vector <double>> C;
    
    
};





#endif /* clsAmericanFlatVol1DPDE_hpp */
