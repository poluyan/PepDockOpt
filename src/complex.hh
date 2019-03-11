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
#ifndef INCLUDED_complex_hh
#define INCLUDED_complex_hh

#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <protocols/simple_filters/AtomicDistanceFilter.hh>
#include <core/kinematics/FoldTree.hh>

#include <set>

#include <get_dof.hh>

namespace pepdockopt
{

class ComplexInfoNseq
{
public:
    ComplexInfoNseq();
    void set_pose(core::pose::Pose const &_pose);
    void show_peptide_info();

    void set_ACE_to_peptide_N_term();

    numeric::Size get_protein_first_index() const;
    numeric::Size get_protein_last_index() const;

    std::vector<numeric::Size> get_protein_first_indices() const;
    std::vector<numeric::Size> get_protein_last_indices() const;

    numeric::Size get_peptide_first_index() const;
    numeric::Size get_peptide_last_index() const;
    core::pose::Pose get_pose_complex() const;
    core::pose::Pose get_ideal_peptide() const;
    core::pose::Pose get_peptide() const;
    core::pose::Pose get_protein() const;

    void extend_peptide();

    void make_ideal_peptide();
    void make_ideal_protein();

    void set_ideal_peptide_values_to_pose_complex(bool include_proline_chi);

    void move_peptide_to_protein_center(bool CA_center);
    double get_max_distance_from_center_of_protein() const;
    double get_extended_peptide_length() const;
    std::set<core::Size> get_protein_core_resiudes(core::Real pore_radius, core::Real core_delta, bool all_atoms) const;
    void remove_core_resiudes(std::vector<core::Size> &cm_res, std::set<core::Size> &core_res) const;

    std::vector<core::Size> get_protein_index_from_peptide_protein_contact_map(core::Real dist);
    std::vector<core::Size> get_protein_index_between_two_spheres(core::Vector a, core::Vector b, core::Real dist, core::Size steps);

    std::vector<opt_element> get_phi_psi_with_info(core::pose::Pose &pose, std::vector<numeric::Size> &ind, size_t chain_id, bool remove_last);
    std::vector<opt_element> get_omega_with_info(core::pose::Pose &pose, std::vector<numeric::Size> &ind, size_t chain_id, bool remove_last);
    std::vector<opt_element> get_peptide_all_chi_dof_with_info(core::pose::Pose &pose, std::vector<numeric::Size> &ind, size_t chain_id, bool include_proline);

protected:
    core::pose::PoseOP pose_complex;

    core::pose::PoseOP protein;
    core::pose::PoseOP peptide;

    core::pose::Pose ideal_peptide;
    core::pose::Pose ideal_protein;

    std::vector<numeric::Size> protein_first_indices;
    std::vector<numeric::Size> protein_last_indices;
    numeric::Size protein_first_index;
    numeric::Size protein_last_index;

    numeric::Size peptide_first_index;
    numeric::Size peptide_last_index;
};

}

#endif
