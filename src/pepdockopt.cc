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
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <pepdockopt.hh>
#include <complex.hh>

#include <vector>
#include <omp.h>

namespace pepdockopt
{

PepDockOpt::PepDockOpt()
{

}

void PepDockOpt::init()
{
    core::pose::Pose input, native;
    if(basic::options::option[basic::options::OptionKeys::in::file::s].user())
    {
        std::string fname = basic::options::option[basic::options::OptionKeys::in::file::s].value_string();
        std::cout << "input pose " << fname << std::endl;
        core::import_pose::pose_from_file(input, fname);
        input.dump_pdb("output/pdb/input.pdb");
    }
    else
    {
        std::cout << "native pose not found" << std::endl;
    }
    if(basic::options::option[basic::options::OptionKeys::in::file::native].user())
    {
        std::string fname = basic::options::option[basic::options::OptionKeys::in::file::native].value_string();
        std::cout << "native pose " << fname << std::endl;
        core::import_pose::pose_from_file(native, fname);
        native.dump_pdb("output/pdb/native.pdb");
    }
    else
    {
        std::cout << "native pose not found" << std::endl;
    }
    param_list.set_pose(input);

    //

    param_list.extend_peptide();

    std::cout << param_list.get_max_distance_from_center_of_protein() << std::endl;
    std::cout << param_list.get_extended_peptide_length() << std::endl;

    pose = param_list.get_pose_complex(); //param_list.get_peptide();

    std::vector<core::Size> protein_full_index, peptide_full_index;
    //for( size_t i = param_list.get_protein_first_index(); i <= param_list.get_protein_last_index(); i++ )
    //    protein_full_index.push_back( i );
    for(size_t i = param_list.get_peptide_first_index(); i <= param_list.get_peptide_last_index(); i++)
        peptide_full_index.push_back(i);

    std::vector<pepdockopt::opt_element> peptide_mc_p_info = pepdockopt::get_phi_psi_with_info(pose, peptide_full_index, 2, true);
    //std::vector<opt_element> peptide_mc_o_info = get_omega_with_info(pose, peptide_full_index, 2, true);
    //std::vector<opt_element> peptide_all_chi_info = get_peptide_all_chi_dof_with_info(pose, peptide_full_index, 2, true);
    peptide_mc_p_info.pop_back();

    peptide_ranges.phipsi = std::make_tuple(true, opt_vector.size(), opt_vector.size() + peptide_mc_p_info.size());
    opt_vector.insert(opt_vector.end(), peptide_mc_p_info.begin(), peptide_mc_p_info.end());


    std::vector<core::Size> protein_cm_index;

    size_t use_box = 2; // 0 - sphere, 1 - box, 2 - trans;
    switch(use_box)
    {
        case 0:
        {
            protein_cm_index = param_list.get_protein_index_from_peptide_protein_contact_map(1.0); //5.0
            break;
        }
        case 1:
        {
            protein_cm_index = param_list.get_protein_index_from_peptide_protein_contact_map(60.0); //5.0
            break;
        }
        case 2:
        {
            //do_trans = true;
            //trans_spheres_obj.load_data("input_pdb/DR1-RA_0001.dat", 100);
            trans_spheres_obj.load_data("input/pdb/1JWG_0001.dat", 100); //1JWG_0001_3
            //trans_spheres_obj.load_data("input_pdb/1HQW_0003.dat", 100);

            // from file .dat
            protein_cm_index = param_list.get_protein_index_between_two_spheres(
                                   core::Vector(trans_spheres_obj.spheres.front().x, trans_spheres_obj.spheres.front().y, trans_spheres_obj.spheres.front().z),
                                   core::Vector(trans_spheres_obj.spheres.back().x, trans_spheres_obj.spheres.back().y, trans_spheres_obj.spheres.back().z),
                                   trans_spheres_obj.max_r,
                                   1000);

            // from pose first
            //protein_cm_index = param_list.get_protein_index_between_two_spheres(
            //                       core::Vector(pose.residue(param_list.get_peptide_first_index()).xyz("CA").x(), pose.residue(param_list.get_peptide_first_index()).xyz("CA").y(), pose.residue(param_list.get_peptide_first_index()).xyz("CA").z()),
            //                       core::Vector(pose.residue(param_list.get_peptide_last_index()).xyz("CA").x(), pose.residue(param_list.get_peptide_last_index()).xyz("CA").y(), pose.residue(param_list.get_peptide_last_index()).xyz("CA").z()),
            //                       trans_spheres_obj.max_r,
            //                       1000);


            //std::cout << "center: " << trans_spheres_obj.spheres.back().x << '\t' << trans_spheres_obj.spheres.back().y << '\t' << trans_spheres_obj.spheres.back().z << std::endl;

            break;
        }
    }
}

void PepDockOpt::set_number_of_threads(size_t n)
{
    threads_number = n;
    if(omp_get_num_procs() < threads_number)
    {
        threads_number = omp_get_num_procs();
    }
    std::cout << "number_of_threads: " << threads_number << std::endl;
}


}
