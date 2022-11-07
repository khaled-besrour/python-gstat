//
//  random_generator.cpp
//  kriging
//
//  Created by khaled besrour on 19/10/2022.
//

#include "random_generator.hpp"

void RandomGenerator::set_seed(int seed){
    this->gen32.seed(seed);
}


double RandomGenerator::unif_rand(){
    return this->ud(this->gen32);
}

double RandomGenerator::norm_rand(){
    return this->nd(this->gen32);
}

RandomGenerator random_generator;
