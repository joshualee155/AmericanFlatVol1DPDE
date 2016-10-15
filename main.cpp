//
//  main.cpp
//  AmericanFlatVol1DPDE
//
//  Created by Xiangyu Joshua Li on 18/09/2015.
//  Copyright © 2015 Xiangyu Joshua Li. All rights reserved.
//

#include <iostream>
#include "clsAmericanFlatVol1DPDE.hpp"

int main(int argc, const char * argv[]) {
    
    clsAmericanFlatVol1DPDE my_option(100, 100, 0.01, 0.5, 0.20, 0.005, true,
                                      400, 200, 0.0, 8);
    my_option.run();
    
    std::cout << my_option.get_price(100) << std::endl;
    
    return 0;
}
