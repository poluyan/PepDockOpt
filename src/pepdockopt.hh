/**************************************************************************

   Copyright © 2019 Sergey Poluyan <svpoluyan@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

**************************************************************************/
#ifndef INCLUDED_pepdockopt_hh
#define INCLUDED_pepdockopt_hh

#include <core/pose/Pose.hh>

#include <complex.hh>
#include <spheres.hh>

#include <iostream>

namespace pepdockopt
{

class PepDockOpt
{
    int threads_number;
    core::pose::Pose pose;
    core::pose::Pose native;
    pepdockopt::ComplexInfoNseq param_list;
    
    std::vector<pepdockopt::opt_element> opt_vector;
    pepdockopt::ranges peptide_ranges;
    pepdockopt::ranges protein_ranges;
    
    pepdockopt::spheres::box_trans trans_spheres_obj;
public:
    PepDockOpt();
    void init();
    void set_number_of_threads(size_t n);

};

}


#endif
