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
#ifndef INCLUDED_transform_hh
#define INCLUDED_transform_hh

#include <opt.hh>
#include <bbdep_sm.hh>
#include <spheres.hh>

#include <tuple>

namespace pepdockopt
{

namespace transform
{

//struct ranges
//{
//    std::tuple<bool, size_t, size_t> phipsi;
//    std::tuple<bool, size_t, size_t> omega;
//    std::tuple<bool, size_t, size_t> chi;
//
//    bool do_phipsi;
//    bool do_omega;
//    bool do_chi;
//};

std::vector<double> bbdep_experiment_actual_states_peptide(
    std::vector<double> x,
    const std::vector<pepdockopt::opt_element> &opt_vect,
    const pepdockopt::ranges &range,
    const pepdockopt::bbdep::BBDEP_Dunbrack_sm &bbdep_obj_sm,
    size_t peptide_first_index,
    size_t peptide_last_index);

std::vector<double> bbdep_experiment_actual_states_protein(
    std::vector<double> x,
    const std::vector<pepdockopt::opt_element> &opt_vect,
    const pepdockopt::ranges &range,
    const pepdockopt::bbdep::BBDEP_Dunbrack_sm &bbdep_obj_sm,
    const std::map<core::Size, std::pair<double, double>> &cm_fixed_phipsi,
    const std::vector<core::Size> &protein_first_indices,
    const std::vector<core::Size> &protein_last_indices);

std::vector<double> peptide_quaternion(std::vector<double> x, const std::vector<opt_element> &opt_vect, size_t opt_vect_size);

std::vector<double> twospheres(std::vector<double> x, size_t opt_vect_size, const spheres::box_trans &spheres_obj);

}
}

#endif
