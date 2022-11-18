//
//  python_wrapper.cpp
//  kriging
//
//  Created by khaled besrour on 21/10/2022.
//


#include <stdio.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

#include <string>
#include <vector>
#include "../gstat.hpp"

std::vector<double> gstat_predict_python(long sn,  std::vector<double> & slocs,
                                                 std::vector<double> & sX, std::vector<int> & block_cols,
                                                 std::vector<double> & block, std::vector<double> & weights,
                                                 int nsim, long blue){

    std::vector<double> retvector;
    std::vector<int> retvector_dim;

    gstat_predict(retvector, retvector_dim,  sn,  slocs,
                   sX, block_cols,
                   block,  weights,
                  nsim,  blue);


    return retvector;

}
PYBIND11_MODULE(gstat, m) {
    m.doc() = "gstat python plugin"; // optional module docstring

    m.def("gstat_init", &gstat_init, "A function that init gstat and random seed");
    m.def("gstat_new_data", &gstat_new_data, "Load data in gstat");
    m.def("gstat_load_variogram", &gstat_load_variogram, "Load gstat_load_variogram in gstat");
    m.def("gstat_predict", &gstat_predict_python, "Predict kriging");
    m.def("gstat_exit", &gstat_exit, "Exit clear memory");
}

