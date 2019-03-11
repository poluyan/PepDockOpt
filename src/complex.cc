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
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <protocols/simple_moves/DeleteChainMover.hh>

#include <core/pose/annotated_sequence.hh>

#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <core/scoring/sasa.hh>

#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/pose/variant_util.hh>

#include <complex.hh>

namespace pepdockopt
{

ComplexInfoNseq::ComplexInfoNseq()
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
    set_pose(input);
}

void ComplexInfoNseq::set_pose(core::pose::Pose const &_pose)
{
    pose_complex = core::pose::PoseOP(new core::pose::Pose(_pose));

    peptide = pose_complex->split_by_chain(pose_complex->num_chains());

    protein = pose_complex->clone();

    protocols::simple_moves::DeleteChainMover DeleteChainMover;
    DeleteChainMover.chain_num(pose_complex->num_chains());
    DeleteChainMover.apply(*protein);

    if(protein->total_residue() < peptide->total_residue() || peptide->total_residue() > 15)
    {
        std::cout << "error" << std::endl;
    }
    else
    {
        core::kinematics::FoldTree ft = pose_complex->fold_tree();
        protein_first_index = ft.upstream_jump_residue(pose_complex->num_jump());
        protein_last_index = ft.downstream_jump_residue(pose_complex->num_jump()) - 1;

        for(size_t i = 0; i < pose_complex->num_jump(); i++)
        {
            protein_last_indices.push_back(ft.downstream_jump_residue(i + 1) - 1);
        }
        protein_first_indices.push_back(1);
        for(size_t i = 1; i < protein_last_indices.size(); i++)
        {
            protein_first_indices.push_back(protein_last_indices[i - 1] + 1);
        }

        peptide_first_index = protein->total_residue() + 1;
        peptide_last_index = pose_complex->total_residue();
        std::cout << "\nComplex" << std::endl;
        //std::cout << "protein: " << protein_first_index << '-' << protein_last_index << ", length " << protein->total_residue() << std::endl;
        //std::cout << "peptide: " << peptide_first_index << '-' << peptide_last_index << ", length " << peptide->total_residue() << std::endl;
        std::cout << "total atoms: " << pose_complex->total_atoms() << std::endl;
        std::cout << "protein atoms: " << protein->total_atoms() << std::endl;
        std::cout << "peptide atoms: " << peptide->total_atoms() << '\n'
                  << std::endl;
        make_ideal_peptide();
        make_ideal_protein();
    }
}

void ComplexInfoNseq::make_ideal_peptide()
{
    ideal_peptide = *peptide;
    core::conformation::Conformation conf = ideal_peptide.conformation();

    for(numeric::Size seqpos = 1, seqpos_end = ideal_peptide.total_residue(); seqpos <= seqpos_end; ++seqpos)
    {
        core::conformation::idealize_position(seqpos, ideal_peptide.conformation());
    }
}
void ComplexInfoNseq::make_ideal_protein()
{
    ideal_protein = *protein;
    core::conformation::Conformation conf = ideal_protein.conformation();

    for(numeric::Size seqpos = 1, seqpos_end = ideal_protein.total_residue(); seqpos <= seqpos_end; ++seqpos)
    {
        core::conformation::idealize_position(seqpos, ideal_protein.conformation());
    }
}

void ComplexInfoNseq::set_ideal_peptide_values_to_pose_complex(bool include_proline_chi)
{
    for(numeric::Size seqpos = peptide_first_index, seqpos_end = peptide_last_index, pep_seqpos = 1; seqpos <= seqpos_end; ++seqpos, ++pep_seqpos)
    {
        pose_complex->set_phi(seqpos, ideal_peptide.phi(pep_seqpos));
        pose_complex->set_psi(seqpos, ideal_peptide.psi(pep_seqpos));
        pose_complex->set_omega(seqpos, ideal_peptide.omega(pep_seqpos));
    }
    for(numeric::Size seqpos = peptide_first_index, seqpos_end = peptide_last_index, pep_seqpos = 1; seqpos <= seqpos_end; ++seqpos, ++pep_seqpos)
    {
        for(numeric::Size atom = 1, atom_end = ideal_peptide.residue(pep_seqpos).natoms(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(atom, seqpos);
            core::id::AtomID aid_pep(atom, pep_seqpos);
            if(pose_complex->conformation().atom_tree().atom(aid).is_jump()) // skip non-bonded atoms
            {
                continue;
            }
            pose_complex->set_dof(core::id::DOF_ID(aid, core::id::THETA), ideal_peptide.dof(core::id::DOF_ID(aid_pep, core::id::THETA)));
            if(atom > 3)
                pose_complex->set_dof(core::id::DOF_ID(aid, core::id::PHI), ideal_peptide.dof(core::id::DOF_ID(aid_pep, core::id::PHI)));
            pose_complex->set_dof(core::id::DOF_ID(aid, core::id::D), ideal_peptide.dof(core::id::DOF_ID(aid_pep, core::id::D)));

            if((pose_complex->residue(seqpos).type().aa() == core::chemical::aa_pro || pose_complex->residue(seqpos).type().aa() == core::chemical::aa_dpr) && include_proline_chi)
            {
                core::pose::Pose proline;
                core::pose::make_pose_from_sequence(proline, "P", *core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));

                for(core::Size i = 1; i <= pose_complex->residue(seqpos).nchi(); i++)
                {
                    pose_complex->set_chi(i, seqpos, proline.chi(i, 1));
                }
            }
        }
    }
}

numeric::Size ComplexInfoNseq::get_peptide_first_index() const
{
    return peptide_first_index;
}

numeric::Size ComplexInfoNseq::get_peptide_last_index() const
{
    return peptide_last_index;
}

numeric::Size ComplexInfoNseq::get_protein_first_index() const
{
    return protein_first_index;
}
numeric::Size ComplexInfoNseq::get_protein_last_index() const
{
    return protein_last_index;
}

core::pose::Pose ComplexInfoNseq::get_pose_complex() const
{
    return *pose_complex;
}

std::vector<core::Size> ComplexInfoNseq::get_protein_index_from_peptide_protein_contact_map(core::Real dist)
{
    std::vector<core::Size> selected_protein_residues;
    for(core::Size peptide = peptide_first_index, peptide_end = peptide_last_index; peptide <= peptide_end; ++peptide)
    {
        for(core::Size seqpos = 1, seqpos_end = protein->total_residue(); seqpos <= seqpos_end; ++seqpos)
        {
            protocols::simple_filters::AtomicDistanceFilter ADF(peptide, seqpos, "CA", "CA", 0, 0, 4.0);
            if(ADF.compute(*pose_complex) < dist)
                selected_protein_residues.push_back(seqpos);
            //std::cout << seqpos << '\t' << std::fixed << ADF.compute(pose) << std::endl;
            //ADF.report(std::cout, pose);
        }
    }
    std::sort(selected_protein_residues.begin(), selected_protein_residues.end());
    selected_protein_residues.erase(std::unique(selected_protein_residues.begin(), selected_protein_residues.end()), selected_protein_residues.end());
    return selected_protein_residues;
}
std::vector<core::Size> ComplexInfoNseq::get_protein_index_between_two_spheres(core::Vector a, core::Vector b, core::Real radius, core::Size steps)
{
    std::vector<core::Size> selected_protein_residues;
    core::Vector direction((b - a).normalize());
    core::Real dist = a.distance(b);
    for(size_t i = 0; i != steps + 1; i++)
    {
        core::Vector newSpot = a + (direction * (dist * i / steps));
        for(core::Size seqpos = protein_first_index; seqpos <= protein_last_index; ++seqpos)
        {
            for(numeric::Size atom = 1, atom_end = protein->residue(seqpos).natoms(); atom <= atom_end; ++atom)
            {
                core::id::AtomID aid(atom, seqpos);
                if(protein->conformation().atom_tree().atom(aid).is_jump()) // skip non-bonded atoms
                {
                    continue;
                }
                if(newSpot.distance(protein->conformation().xyz(aid)) < radius)
                {
                    selected_protein_residues.push_back(seqpos);
                    break;
                }
            }
        }
    }
    std::sort(selected_protein_residues.begin(), selected_protein_residues.end());
    selected_protein_residues.erase(std::unique(selected_protein_residues.begin(), selected_protein_residues.end()), selected_protein_residues.end());
    return selected_protein_residues;
}
std::vector<numeric::Size> ComplexInfoNseq::get_protein_first_indices() const
{
    return protein_first_indices;
}
std::vector<numeric::Size> ComplexInfoNseq::get_protein_last_indices() const
{
    return protein_last_indices;
}

std::vector<opt_element> ComplexInfoNseq::get_phi_psi_with_info(core::pose::Pose &pose, std::vector<numeric::Size> &ind, size_t chain_id, bool remove_last)
{
    return get_phi_psi_with_info(pose, ind, chain_id, remove_last);
}

std::vector<opt_element> ComplexInfoNseq::get_omega_with_info(core::pose::Pose &pose, std::vector<numeric::Size> &ind, size_t chain_id, bool remove_last)
{
    return get_omega_with_info(pose, ind, chain_id, remove_last);
}

std::vector<opt_element> ComplexInfoNseq::get_peptide_all_chi_dof_with_info(core::pose::Pose &pose, std::vector<numeric::Size> &ind, size_t chain_id, bool include_proline)
{
    return get_peptide_all_chi_dof_with_info(pose, ind, chain_id, include_proline);
}

void ComplexInfoNseq::move_peptide_to_protein_center(bool CA_center)
{
    /// This move depend on conformation!! Different conformation = different centroid!!!
    /// Using core::kinematics::Jump flexible_jump = pose.jump(pose.num_jump()); -- is independent from conformation

    //std::cout << protocols::geometry::center_of_mass(pose, param_list.get_protein_first_index(), param_list.get_protein_last_index()).to_string() << std::endl;

    //protocols::rigid::RigidBodyTransMover RBTM = protocols::rigid::RigidBodyTransMover(*pose_complex, pose_complex->num_jump());
    //core::Vector upstream_dummy, downstream_dummy; //upstream_dummy -- center of protein, downstream_dummy -- center of peptide
    //protocols::geometry::centroids_by_jump(*pose_complex, pose_complex->num_jump(), upstream_dummy, downstream_dummy);
    //numeric::Real current_distance(upstream_dummy.distance(downstream_dummy));
    //RBTM.trans_axis(upstream_dummy - downstream_dummy);
    //RBTM.step_size(current_distance);
    //RBTM.apply(*pose_complex);

    utility::vector1<core::Size> protein_chains;
    for(core::Size i = 1; i < pose_complex->num_chains(); i++)
        protein_chains.push_back(i);

    core::Vector protein_centroid = protocols::geometry::centroid_by_chains(*pose_complex, protein_chains);
    core::Vector peptide_centroid = protocols::geometry::centroid_by_chain(*pose_complex, pose_complex->num_chains());

    protocols::rigid::RigidBodyTransMover RBTM = protocols::rigid::RigidBodyTransMover(*pose_complex, pose_complex->num_jump());

    if(CA_center)
    {
        numeric::Real current_distance(protein_centroid.distance(pose_complex->residue(peptide_first_index).xyz("CA")));
        RBTM.trans_axis(protein_centroid - pose_complex->residue(peptide_first_index).xyz("CA"));
        RBTM.step_size(current_distance);
        RBTM.apply(*pose_complex);
    }
    else
    {
        numeric::Real current_distance(protein_centroid.distance(peptide_centroid));
        RBTM.trans_axis(protein_centroid - peptide_centroid);
        RBTM.step_size(current_distance);
        RBTM.apply(*pose_complex);
    }
}
double ComplexInfoNseq::get_max_distance_from_center_of_protein() const
{
    utility::vector1<core::Size> protein_chains;
    for(core::Size i = 1; i < pose_complex->num_chains(); i++)
        protein_chains.push_back(i);

    core::Vector protein_centroid = protocols::geometry::centroid_by_chains(*pose_complex, protein_chains);

    numeric::Real max_distance(0.0);
    for(numeric::Size seqpos = protein_first_index, seqpos_end = protein_last_index; seqpos <= seqpos_end; ++seqpos)
    {
        for(numeric::Size atom = 1, atom_end = pose_complex->residue(seqpos).natoms(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(atom, seqpos);
            numeric::Real sd(protein_centroid.distance(pose_complex->conformation().atom_tree().atom(aid).xyz()));
            if(max_distance < sd)
            {
                max_distance = sd;
            }
        }
    }
    return max_distance;
}

double ComplexInfoNseq::get_extended_peptide_length() const
{
    core::pose::Pose temp_pose = *peptide;
    for(numeric::Size i = 1; i < temp_pose.total_residue(); i++)
    {
        temp_pose.set_phi(i, -135.0);
        temp_pose.set_psi(i, 135.0);
        temp_pose.set_omega(i, 179.8);
    }
    //    numeric::Real extended_peptide_length(temp_pose.residue(1).atom("CA").xyz().distance(temp_pose.residue(temp_pose.total_residue()).atom("CA").xyz()));
    //    return extended_peptide_length;

    numeric::Real max_distance(0.0);
    for(numeric::Size atom1 = 1, atom_end1 = temp_pose.residue(1).natoms(); atom1 <= atom_end1; ++atom1)
    {
        core::id::AtomID aid1(atom1, 1);
        for(numeric::Size atom2 = 1, atom_end2 = temp_pose.residue(temp_pose.total_residue()).natoms(); atom2 <= atom_end2; ++atom2)
        {
            core::id::AtomID aid2(atom2, temp_pose.total_residue());
            numeric::Real sd(temp_pose.conformation().atom_tree().atom(aid1).xyz().distance(temp_pose.conformation().atom_tree().atom(aid2).xyz()));
            if(max_distance < sd)
            {
                max_distance = sd;
            }
        }
    }
    return max_distance;
}
core::pose::Pose ComplexInfoNseq::get_peptide() const
{
    return *peptide;
}
core::pose::Pose ComplexInfoNseq::get_protein() const
{
    return *protein;
}

std::set<core::Size> ComplexInfoNseq::get_protein_core_resiudes(core::Real pore_radius, core::Real core_delta, bool all_atoms) const
{
    ///core::Real pore_radius = 1.4;
    ///core::Real core_delta = 0.1;

    std::set<core::Size> core_residues;
    core::id::AtomID_Map<bool> atom_map;

    core::pose::Pose pose = *protein;

    if(all_atoms)
    {
        core::pose::initialize_atomid_map(atom_map, pose, true);
    }
    else
    {
        core::pose::initialize_atomid_map(atom_map, pose, false);
        for(core::Size ir = 1; ir <= pose.size(); ++ir)
        {
            for(core::Size j = 1; j <= 5; ++j)
            {
                core::id::AtomID atom(j, ir);
                atom_map.set(atom, true);
            }
        }
    }
    core::id::AtomID_Map<core::Real> atom_sasa;
    utility::vector1<core::Real> sasa_list;
    core::scoring::calc_per_atom_sasa(pose, atom_sasa, sasa_list, pore_radius, false, atom_map);

    for(core::Size i = 1; i <= sasa_list.size(); ++i)
    {
        if(sasa_list[i] < core_delta)
            core_residues.insert(i);
    }
    return core_residues;
}
void ComplexInfoNseq::remove_core_resiudes(std::vector<core::Size> &cm_res, std::set<core::Size> &core_res) const
{
    cm_res.erase(
        std::remove_if(cm_res.begin(), cm_res.end(),
                       [&](const core::Size &t)
    {
        return core_res.find(t) != core_res.end();
    }),
    cm_res.end());
}

void ComplexInfoNseq::set_ACE_to_peptide_N_term()
{
    if(!peptide->residue(peptide_first_index).has_variant_type(core::chemical::ACETYLATED_NTERMINUS_VARIANT))
    {
        core::pose::remove_lower_terminus_type_from_pose_residue(*pose_complex, peptide_first_index); // LOWER_TERMINUS_VARIANT to NO_VARIANT
        core::pose::add_variant_type_to_pose_residue(*pose_complex, core::chemical::ACETYLATED_NTERMINUS_VARIANT, peptide_first_index);

        core::pose::remove_lower_terminus_type_from_pose_residue(*peptide, 1); // LOWER_TERMINUS_VARIANT to NO_VARIANT
        core::pose::add_variant_type_to_pose_residue(*peptide, core::chemical::ACETYLATED_NTERMINUS_VARIANT, 1);

        core::pose::remove_lower_terminus_type_from_pose_residue(ideal_peptide, 1); // LOWER_TERMINUS_VARIANT to NO_VARIANT
        core::pose::add_variant_type_to_pose_residue(ideal_peptide, core::chemical::ACETYLATED_NTERMINUS_VARIANT, 1);
    }
    else
    {
        std::cout << "already has ACE" << std::endl;
    }
}
core::pose::Pose ComplexInfoNseq::get_ideal_peptide() const
{
    return ideal_peptide;
}

void ComplexInfoNseq::extend_peptide()
{
    for(numeric::Size i = peptide_first_index; i < peptide_last_index; i++)
    {
        pose_complex->set_phi(i, -135.0);
        pose_complex->set_psi(i, 135.0);
        pose_complex->set_omega(i, 179.8);
    }
}

}
