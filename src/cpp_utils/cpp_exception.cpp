//
//  cpp_exception.cpp
//  kriging
//
//  Created by khaled besrour on 19/10/2022.
//
#include <stdio.h>
#include <stdexcept>
#include "cpp_exception.hpp"

extern "C" void run_time_error(const char * message){
    throw std::runtime_error(message);
}
