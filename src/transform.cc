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
#include <transform.hh>
#include <bbtools.hh>
#include <opt.hh>
#include <bbdep_sm.hh>

#include <numeric/conversions.hh>

namespace pepdockopt
{
namespace transform
{

std::vector<double> bbdep_experiment_actual_states_peptide(
    std::vector<double> x,
    const std::vector<opt_element> &opt_vect,
    const pepdockopt::ranges &range,
    const pepdockopt::bbdep::BBDEP_Dunbrack_sm &bbdep_obj_sm,
    size_t peptide_first_index,
    size_t peptide_last_index)
{
    if(range.do_chi) // bbdep_experiment
    {
        if(std::get<0>(range.chi))
        {
            for(size_t i = std::get<1>(range.chi); i < std::get<2>(range.chi); i++)
            {
                x[i] = -numeric::NumericTraits<core::Real>::pi() + 2.0*x[i]*numeric::NumericTraits<core::Real>::pi();
            }
            double phi = 0, psi = 0;
            for(size_t c = std::get<1>(range.chi); c < std::get<2>(range.chi);)
            {
                pepdockopt::bbdep::Dunbrack_data first_line = bbdep_obj_sm.get_first_line(opt_vect[c].amino_acid);

                if(opt_vect[c].seqpos == peptide_first_index)
                {
                    if(opt_vect[c].nchi == 1)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));
                    }
                    else if(opt_vect[c].nchi == 2)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                        if(first_line.r2 != 0)
                        {
                            double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        1);
                            x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                        }
                    }
                    else if(opt_vect[c].nchi == 3)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                        if(first_line.r2 != 0)
                        {
                            double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        1);
                            x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                        }

                        if(first_line.r3 != 0)
                        {
                            double c3 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 2] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        2);
                            x[c + 2] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
                        }
                    }
                    else if(opt_vect[c].nchi == 4)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                        if(first_line.r2 != 0)
                        {
                            double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        1);
                            x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                        }

                        if(first_line.r3 != 0)
                        {
                            double c3 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 2] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        2);
                            x[c + 2] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
                        }

                        if(first_line.r4 != 0)
                        {
                            double c4 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 3] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        3);
                            x[c + 3] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c4));
                        }
                    }

                    c += opt_vect[c].nchi;
                    continue;
                }
                else if(opt_vect[c].seqpos == peptide_last_index)
                {
                    if(opt_vect[c].nchi == 1)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));
                    }
                    else if(opt_vect[c].nchi == 2)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                        if(first_line.r2 != 0)
                        {
                            double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        1);
                            x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                        }
                    }
                    else if(opt_vect[c].nchi == 3)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                        if(first_line.r2 != 0)
                        {
                            double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        1);
                            x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                        }

                        if(first_line.r3 != 0)
                        {
                            double c3 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 2] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        2);
                            x[c + 2] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
                        }
                    }
                    else if(opt_vect[c].nchi == 4)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                        if(first_line.r2 != 0)
                        {
                            double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        1);
                            x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                        }

                        if(first_line.r3 != 0)
                        {
                            double c3 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 2] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        2);
                            x[c + 2] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
                        }

                        if(first_line.r4 != 0)
                        {
                            double c4 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 3] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        3);
                            x[c + 3] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c4));
                        }
                    }

                    c += opt_vect[c].nchi;
                    continue;
                }
                else
                {
                    numeric::Size sp = opt_vect[c].seqpos;
                    size_t index = std::distance(opt_vect.begin(),
                                                 std::find_if(opt_vect.begin(),
                                                         opt_vect.begin() + std::get<2>(range.phipsi),
                                                         [&sp](const opt_element &s)
                    {
                        return (s.seqpos == sp && s.torsion_name == "bbphi");
                    }));
                    phi = x[index];
                    psi = x[index + 1];
                }

                phi = numeric::conversions::degrees(phi);
                psi = numeric::conversions::degrees(psi);

                size_t index = bbdep_obj_sm.find_index_for_cdf_chi234(opt_vect[c].amino_acid, phi, psi);

                if(opt_vect[c].nchi == 1)
                {
                    double c1 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chinumber(index, opt_vect[c].amino_acid,
                                (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                (2 * numeric::NumericTraits<core::Real>::pi()),
                                0);
                    x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));
                }
                else if(opt_vect[c].nchi == 2)
                {
                    double c1 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chinumber(index, opt_vect[c].amino_acid,
                                (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                (2 * numeric::NumericTraits<core::Real>::pi()),
                                0);
                    x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                    if(first_line.r2 != 0)
                    {
                        double c2 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi1_dep_actual_states(index,
                                    opt_vect[c].amino_acid, numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c])),
                                    (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()));
                        x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                    }
                }
                else if(opt_vect[c].nchi == 3)
                {
                    double c1 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chinumber(index, opt_vect[c].amino_acid,
                                (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                (2 * numeric::NumericTraits<core::Real>::pi()),
                                0);
                    x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                    if(first_line.r2 != 0)
                    {
                        double c2 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi1_dep_actual_states(index,
                                    opt_vect[c].amino_acid, numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c])),
                                    (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()));
                        x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                    }

                    if(first_line.r3 != 0)
                    {
                        double c3 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi12_dep_actual_states(index,
                                    opt_vect[c].amino_acid, numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c])),
                                    numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c + 1])),
                                    (x[c + 2] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()));
                        x[c + 2] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
                    }
                }
                else if(opt_vect[c].nchi == 4)
                {
                    double c1 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chinumber(index, opt_vect[c].amino_acid,
                                (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                (2 * numeric::NumericTraits<core::Real>::pi()),
                                0);
                    x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                    if(first_line.r2 != 0)
                    {
                        double c2 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi1_dep_actual_states(index,
                                    opt_vect[c].amino_acid, numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c])),
                                    (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()));
                        x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                    }

                    if(first_line.r3 != 0)
                    {
                        double c3 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi12_dep_actual_states(index,
                                    opt_vect[c].amino_acid, numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c])),
                                    numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c + 1])),
                                    (x[c + 2] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()));
                        x[c + 2] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
                    }

                    if(first_line.r4 != 0)
                    {
                        double c4 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi123_dep_actual_states(index,
                                    opt_vect[c].amino_acid, numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c])),
                                    numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c + 1])),
                                    numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c + 2])),
                                    (x[c + 3] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()));
                        x[c + 3] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c4));
                    }
                }
                c += opt_vect[c].nchi;
            }
        }
    }
    return x;
}

std::vector<double> bbdep_experiment_actual_states_protein(
    std::vector<double> x,
    const std::vector<opt_element> &opt_vect,
    const pepdockopt::ranges &range,
    const pepdockopt::bbdep::BBDEP_Dunbrack_sm &bbdep_obj_sm,
    const std::map<core::Size, std::pair<double, double>> &cm_fixed_phipsi,
    const std::vector<core::Size> &protein_first_indices,
    const std::vector<core::Size> &protein_last_indices)
{
    if(range.do_chi) // bbdep_experiment
    {
        if(std::get<0>(range.chi))
        {
            for(size_t i = std::get<1>(range.chi); i < std::get<2>(range.chi); i++)
            {
                x[i] = -numeric::NumericTraits<core::Real>::pi() + 2.0*x[i]*numeric::NumericTraits<core::Real>::pi();
            }
            double phi = 0, psi = 0;
            for(size_t c = std::get<1>(range.chi); c < std::get<2>(range.chi);)
            {
                pepdockopt::bbdep::Dunbrack_data first_line = bbdep_obj_sm.get_first_line(opt_vect[c].amino_acid);

                if(std::find(protein_first_indices.begin(), protein_first_indices.end(), opt_vect[c].seqpos) != protein_first_indices.end())
                {
                    if(opt_vect[c].nchi == 1)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));
                    }
                    else if(opt_vect[c].nchi == 2)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                        if(first_line.r2 != 0)
                        {
                            double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        1);
                            x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                        }
                    }
                    else if(opt_vect[c].nchi == 3)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                        if(first_line.r2 != 0)
                        {
                            double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        1);
                            x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                        }

                        if(first_line.r3 != 0)
                        {
                            double c3 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 2] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        2);
                            x[c + 2] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
                        }
                    }
                    else if(opt_vect[c].nchi == 4)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                        if(first_line.r2 != 0)
                        {
                            double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        1);
                            x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                        }

                        if(first_line.r3 != 0)
                        {
                            double c3 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 2] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        2);
                            x[c + 2] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
                        }

                        if(first_line.r4 != 0)
                        {
                            double c4 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 3] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        3);
                            x[c + 3] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c4));
                        }
                    }

                    c += opt_vect[c].nchi;
                    continue;
                }
                else if(std::find(protein_last_indices.begin(), protein_last_indices.end(), opt_vect[c].seqpos) != protein_last_indices.end())
                {
                    if(opt_vect[c].nchi == 1)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));
                    }
                    else if(opt_vect[c].nchi == 2)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                        if(first_line.r2 != 0)
                        {
                            double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        1);
                            x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                        }
                    }
                    else if(opt_vect[c].nchi == 3)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                        if(first_line.r2 != 0)
                        {
                            double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        1);
                            x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                        }

                        if(first_line.r3 != 0)
                        {
                            double c3 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 2] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        2);
                            x[c + 2] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
                        }
                    }
                    else if(opt_vect[c].nchi == 4)
                    {
                        double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                    (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()),
                                    0);
                        x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                        if(first_line.r2 != 0)
                        {
                            double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        1);
                            x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                        }

                        if(first_line.r3 != 0)
                        {
                            double c3 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 2] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        2);
                            x[c + 2] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
                        }

                        if(first_line.r4 != 0)
                        {
                            double c4 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid,
                                        (x[c + 3] + numeric::NumericTraits<core::Real>::pi()) /
                                        (2 * numeric::NumericTraits<core::Real>::pi()),
                                        3);
                            x[c + 3] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c4));
                        }
                    }

                    c += opt_vect[c].nchi;
                    continue;
                }
                else
                {
                    auto it = cm_fixed_phipsi.find(opt_vect[c].seqpos);
                    phi = it->second.first;
                    psi = it->second.second;
                }

//                phi = numeric::conversions::degrees(phi);
//                psi = numeric::conversions::degrees(psi);

                size_t index = bbdep_obj_sm.find_index_for_cdf_chi234(opt_vect[c].amino_acid, phi, psi);

                if(opt_vect[c].nchi == 1)
                {
                    double c1 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chinumber(index, opt_vect[c].amino_acid,
                                (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                (2 * numeric::NumericTraits<core::Real>::pi()),
                                0);
                    x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));
                }
                else if(opt_vect[c].nchi == 2)
                {
                    double c1 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chinumber(index, opt_vect[c].amino_acid,
                                (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                (2 * numeric::NumericTraits<core::Real>::pi()),
                                0);
                    x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                    if(first_line.r2 != 0)
                    {
                        double c2 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi1_dep_actual_states(index,
                                    opt_vect[c].amino_acid, numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c])),
                                    (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()));
                        x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                    }
                }
                else if(opt_vect[c].nchi == 3)
                {
                    double c1 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chinumber(index, opt_vect[c].amino_acid,
                                (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                (2 * numeric::NumericTraits<core::Real>::pi()),
                                0);
                    x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                    if(first_line.r2 != 0)
                    {
                        double c2 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi1_dep_actual_states(index,
                                    opt_vect[c].amino_acid, numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c])),
                                    (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()));
                        x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                    }

                    if(first_line.r3 != 0)
                    {
                        double c3 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi12_dep_actual_states(index,
                                    opt_vect[c].amino_acid, numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c])),
                                    numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c + 1])),
                                    (x[c + 2] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()));
                        x[c + 2] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
                    }
                }
                else if(opt_vect[c].nchi == 4)
                {
                    double c1 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chinumber(index, opt_vect[c].amino_acid,
                                (x[c] + numeric::NumericTraits<core::Real>::pi()) /
                                (2 * numeric::NumericTraits<core::Real>::pi()),
                                0);
                    x[c] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

                    if(first_line.r2 != 0)
                    {
                        double c2 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi1_dep_actual_states(index,
                                    opt_vect[c].amino_acid, numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c])),
                                    (x[c + 1] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()));
                        x[c + 1] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
                    }

                    if(first_line.r3 != 0)
                    {
                        double c3 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi12_dep_actual_states(index,
                                    opt_vect[c].amino_acid, numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c])),
                                    numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c + 1])),
                                    (x[c + 2] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()));
                        x[c + 2] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
                    }

                    if(first_line.r4 != 0)
                    {
                        double c4 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi123_dep_actual_states(index,
                                    opt_vect[c].amino_acid, numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c])),
                                    numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c + 1])),
                                    numeric::conversions::degrees(pepdockopt::bbtools::to_positive_radians(x[c + 2])),
                                    (x[c + 3] + numeric::NumericTraits<core::Real>::pi()) /
                                    (2 * numeric::NumericTraits<core::Real>::pi()));
                        x[c + 3] = pepdockopt::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c4));
                    }
                }
                c += opt_vect[c].nchi;
            }
        }
    }
    return x;
}

std::vector<double> peptide_quaternion(std::vector<double> x, const std::vector<opt_element> &opt_vect, size_t opt_vect_size)
{
    // if search space from [0;1] move it to [-1;1]
    x[opt_vect_size + 3] = 2.0*x[opt_vect_size + 3] - 1.0;
    x[opt_vect_size + 4] = 2.0*x[opt_vect_size + 4] - 1.0;
    x[opt_vect_size + 5] = 2.0*x[opt_vect_size + 5] - 1.0;
    // for angle from [0;1] to [0;180]
    x[opt_vect_size + 6] = 180.0*x[opt_vect_size + 6];
    return x;
}

std::vector<double> twospheres(std::vector<double> x, size_t opt_vect_size, const spheres::box_trans &spheres_obj)
{
    std::vector<double> xyz01 = { x[opt_vect_size], x[opt_vect_size + 1], x[opt_vect_size + 2] };
    std::vector<double> rez = bbutils::get_3d_from_dst(spheres_obj.dist, xyz01);
    x[opt_vect_size] = rez[0];
    x[opt_vect_size + 1] = rez[1];
    x[opt_vect_size + 2] = rez[2];
    return x;
}

}
}
