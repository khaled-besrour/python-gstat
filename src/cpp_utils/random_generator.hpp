//
//  random_generator.hpp
//  kriging
//
//  Created by khaled besrour on 19/10/2022.
//

#ifndef random_generator_hpp
#define random_generator_hpp


#include <iostream>
#include <random>

class RandomGenerator {
    std::mt19937 gen32;
    
    std::uniform_real_distribution<> ud;
    std::normal_distribution<> nd;
    
public:
    void set_seed(int seed);
    
    double unif_rand();
    double norm_rand();
    
};


extern RandomGenerator random_generator;


#endif /* random_generator_hpp */
