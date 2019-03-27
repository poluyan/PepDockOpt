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

#include <numeric/NumericTraits.hh>

#include <pepdockopt.hh>
#include <complex.hh>
#include <transform.hh>
#include <cso.hh>
#include <data_io.hh>

#include <vector>
#include <random>
#include <omp.h>

namespace pepdockopt
{

PepDockOpt::PepDockOpt() {}
template <typename T>
bool increase(const std::vector<std::vector<T>>& v, std::vector<std::size_t>& it)
{
    for(std::size_t i = 0, size = it.size(); i != size; ++i)
    {
        const std::size_t index = size - 1 - i;
        ++it[index];
        if(it[index] == v[index].size())
        {
            it[index] = 0;
        }
        else
        {
            return true;
        }
    }
    return false;
}

template <typename T>
std::vector<T> do_job(const std::vector<std::vector<T>>& v, std::vector<std::size_t>& it)
{
    std::vector<T> rez(v.size());
    for(std::size_t i = 0, size = v.size(); i != size; ++i)
    {
        rez[i] = v[i][it[i]];
    }
    return rez;
}
template <typename T>
std::vector<std::vector<T>> iterate(const std::vector<std::vector<T>>& v)
{
    std::vector<std::size_t> it(v.size(), 0);
    std::vector<std::vector<T>> values;
    do
    {
        values.push_back(do_job(v, it));
    }
    while(increase(v, it));
    std::cout << "Permutations " << values.size() << std::endl;
    return values;
}
std::vector<char> findTheDot_MultipleGrids_Moore(std::vector<std::vector<double>>& grids,
        std::vector<std::vector<char>> &points,
        std::set<std::vector<char>> &visited,
        std::vector<double> &dx,
        std::function<double(std::vector<double>)> f)
{
    /// generate all permutations
    std::vector<std::vector<char>> variable_values(points.back().size() - 7, std::vector<char>(3));
    for(size_t i = 0; i != variable_values.size(); i++)
    {
        variable_values[i][0] = -1;
        variable_values[i][1] = 0;
        variable_values[i][2] = 1;
    }
    std::vector<std::vector<char>> permut = iterate(variable_values);


    std::vector<char> t = points.back(), point;
    points.pop_back();

    for(size_t i = 0; i != t.size(); i++)
    {
        point = t;
        point[i] = point[i] + 1;

        if(point[i] < 0 || point[i] > grids[i].size() - 1)
            continue;

        points.push_back(point);
    }

    for(size_t i = 0; i != t.size(); i++)
    {
        point = t;
        point[i] = point[i] - 1;

        if(point[i] < 0 || point[i] > grids[i].size() - 1)
            continue;

        points.push_back(point);
    }

    for(size_t i = 0; i != t.size(); i++)
    {
        point = t;
        point[i] = point[i] - 1;
        for(size_t j = 0; j != t.size(); j++)
        {
            point[i] = point[i] - 1;

            if(point[i] < 0 || point[i] > grids[i].size() - 1)
                continue;

            points.push_back(point);
        }
    }

    for(size_t i = 0; i != permut.size(); i++)
    {
        point = t;
        bool flag = true;
        for(size_t j = 0; j != variable_values.size(); j++)
        {
            point[j] = point[j] + permut[i][j];
            if(point[j] < 0 || point[j] > grids[j].size() - 1)
            {
                flag = false;
                break;
            }
        }
        if(!flag)
            continue;

        points.push_back(point);
    }

    while(!points.empty())
    {
        t = points.back();
        points.pop_back();

        std::vector<double> dot(t.size());
        for(size_t i = 0; i != dot.size(); i++)
        {
            dot[i] = grids[i][t[i]] + dx[i];
        }

        double val = f(dot);

        if(val > 1.0)
        {
            std::cout << "found!" << std::endl;
            break;
        }
    }
    return t;
}

std::vector<char> findTheDot_MultipleGrids(std::vector<std::vector<double>>& grids,
        std::vector<std::vector<char>> &points,
        std::set<std::vector<char>> &visited,
        std::vector<double> &dx,
        std::function<double(std::vector<double>)> f)
{
    std::vector<char> t = points.back(), point;
    points.pop_back();

    for(size_t i = 0; i != t.size(); i++)
    {
        point = t;
        point[i] = point[i] + 1;

        if(point[i] < 0 || point[i] > grids[i].size() - 1)
            continue;

        points.push_back(point);
    }

    for(size_t i = 0; i != t.size(); i++)
    {
        point = t;
        point[i] = point[i] - 1;

        //if(point[i] < 0 || point[i] > grids[i].size() - 1)
        //    continue;

        if(point[i] < 0)
        {
            point[i] = grids[i].size() - 1;
        }
        if(point[i] > grids[i].size() - 1)
        {
            point[i] = 0;
        }

        points.push_back(point);
    }

    for(size_t i = 0; i != t.size(); i++)
    {
        point = t;
        point[i] = point[i] - 1;
        for(size_t j = 0; j != t.size(); j++)
        {
            point[i] = point[i] - 1;

            //if(point[i] < 0 || point[i] > grids[i].size() - 1)
            //    continue;
            if(point[i] < 0)
            {
                point[i] = grids[i].size() - 1;
            }
            if(point[i] > grids[i].size() - 1)
            {
                point[i] = 0;
            }

            points.push_back(point);
        }
    }

    while(!points.empty())
    {
        t = points.back();
        points.pop_back();

        std::vector<double> dot(t.size());
        for(size_t i = 0; i != dot.size(); i++)
        {
            dot[i] = grids[i][t[i]] + dx[i];
        }

        double val = f(dot);

        if(val > 1.0)
        {
            std::cout << "found!" << std::endl;
            return t;
        }
    }
    std::cout << "still miss the point" << std::endl;
    return t;
}

double Quadratic_ValueSimple3(std::vector<double> x,
                              core::pose::Pose& pose,
                              PoseShift &shift_params,
                              const std::vector<opt_element> &opt_vector,
                              const ComplexInfoNseq &complex_stuff,
                              const pepdockopt::ranges &peptide_ranges,
                              const spheres::box_trans &trans_spheres_obj,
                              const std::pair<size_t, size_t> sphere_number)
{
    double pi = numeric::NumericTraits<core::Real>::pi();
    for(size_t i = std::get<1>(peptide_ranges.phipsi); i != std::get<2>(peptide_ranges.phipsi); i++)
    {
        x[i] = pi*(2.0*x[i] - 1.0);
    }
    for(size_t i = 0, end = opt_vector.size(); i < end; i++)
    {
        pose.set_dof(opt_vector[i].dofid, x[i]);
    }

    x = transform::peptide_quaternion(x, opt_vector, opt_vector.size());
    x = transform::twospheres(x, opt_vector.size(), trans_spheres_obj);

    ///

    core::Vector new_position = shift_params.InitPeptidePosition;
    new_position.at(0) = x[opt_vector.size()];
    new_position.at(1) = x[opt_vector.size() + 1];
    new_position.at(2) = x[opt_vector.size() + 2];

    shift_params.FlexibleJump.reset();

    shift_params.FlexibleJump.set_rotation(shift_params.InitRm);
    pose.set_jump(pose.num_jump(), shift_params.FlexibleJump);

    numeric::Real current_distance(pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA").distance(new_position));
    shift_params.FlexibleJump.translation_along_axis(shift_params.UpstreamStub, new_position - pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA"), current_distance);
    pose.set_jump(pose.num_jump(), shift_params.FlexibleJump);

    core::Vector peptide_c_alpha_centroid(pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA"));
    core::Vector axis(x[opt_vector.size() + 3], x[opt_vector.size() + 4], x[opt_vector.size() + 5]);

    shift_params.SpinMover.rot_center(peptide_c_alpha_centroid);
    shift_params.SpinMover.spin_axis(axis);
    shift_params.SpinMover.angle_magnitude(x.back());
    shift_params.SpinMover.apply(pose);

    ///

    std::vector<double> exponential_center =
    {
        trans_spheres_obj.spheres[sphere_number.second].x,
        trans_spheres_obj.spheres[sphere_number.second].y,
        trans_spheres_obj.spheres[sphere_number.second].z
    };

    core::Vector lCA = pose.residue(complex_stuff.get_peptide_last_index()).xyz("CA");

    double rez = 0;
    double temp = 0;
    for(size_t i = 0; i != exponential_center.size(); i++)
    {
        rez += std::pow(lCA.at(i) - exponential_center[i], 2.0);
    }
    temp = std::sqrt(rez);
    if(temp < trans_spheres_obj.spheres[sphere_number.second].r)
    {
        rez = 1.1;
    }
    else
        rez = 0;
    return rez;
}
double Quadratic_ValueSimple3_omp(std::vector<double> x,
                                  int th_id,
                                  std::vector<core::pose::Pose> &pose,
                                  std::vector<PoseShift> &shift_params,
                                  const std::vector<opt_element> &opt_vector,
                                  const ComplexInfoNseq &complex_stuff,
                                  const pepdockopt::ranges &peptide_ranges,
                                  const spheres::box_trans &trans_spheres_obj,
                                  const std::pair<size_t, size_t> sphere_number)
{
    double pi = numeric::NumericTraits<core::Real>::pi();
    for(size_t i = std::get<1>(peptide_ranges.phipsi); i != std::get<2>(peptide_ranges.phipsi); i++)
    {
        x[i] = pi*(2.0*x[i] - 1.0);
    }
    for(size_t i = 0, end = opt_vector.size(); i < end; i++)
    {
        pose[th_id].set_dof(opt_vector[i].dofid, x[i]);
    }

    x = transform::peptide_quaternion(x, opt_vector, opt_vector.size());
    x = transform::twospheres(x, opt_vector.size(), trans_spheres_obj);

    ///

    core::Vector new_position = shift_params[th_id].InitPeptidePosition;
    new_position.at(0) = x[opt_vector.size()];
    new_position.at(1) = x[opt_vector.size() + 1];
    new_position.at(2) = x[opt_vector.size() + 2];

    shift_params[th_id].FlexibleJump.reset();

    shift_params[th_id].FlexibleJump.set_rotation(shift_params[th_id].InitRm);
    pose[th_id].set_jump(pose[th_id].num_jump(), shift_params[th_id].FlexibleJump);

    numeric::Real current_distance(pose[th_id].residue(complex_stuff.get_peptide_first_index()).xyz("CA").distance(new_position));
    shift_params[th_id].FlexibleJump.translation_along_axis(shift_params[th_id].UpstreamStub, new_position - pose[th_id].residue(complex_stuff.get_peptide_first_index()).xyz("CA"), current_distance);
    pose[th_id].set_jump(pose[th_id].num_jump(), shift_params[th_id].FlexibleJump);

    core::Vector peptide_c_alpha_centroid(pose[th_id].residue(complex_stuff.get_peptide_first_index()).xyz("CA"));
    core::Vector axis(x[opt_vector.size() + 3], x[opt_vector.size() + 4], x[opt_vector.size() + 5]);

    shift_params[th_id].SpinMover.rot_center(peptide_c_alpha_centroid);
    shift_params[th_id].SpinMover.spin_axis(axis);
    shift_params[th_id].SpinMover.angle_magnitude(x.back());
    shift_params[th_id].SpinMover.apply(pose[th_id]);

    ///

    std::vector<double> exponential_center =
    {
        trans_spheres_obj.spheres[sphere_number.second].x,
        trans_spheres_obj.spheres[sphere_number.second].y,
        trans_spheres_obj.spheres[sphere_number.second].z
    };

    core::Vector lCA = pose[th_id].residue(complex_stuff.get_peptide_last_index()).xyz("CA");

    double rez = 0;
    double temp = 0;
    for(size_t i = 0; i != exponential_center.size(); i++)
    {
        rez += std::pow(lCA.at(i) - exponential_center[i], 2.0);
    }
    temp = std::sqrt(rez);
    if(temp < trans_spheres_obj.spheres[sphere_number.second].r)
    {
        rez = 1.1;
    }
    else
        rez = 0;
    return rez;
}
double Quadratic_ValueSimple2(std::vector<double> x,
                              core::pose::Pose& pose,
                              PoseShift &shift_params,
                              const std::vector<opt_element> &opt_vector,
                              const ComplexInfoNseq &complex_stuff,
                              double width,
                              const pepdockopt::ranges &peptide_ranges,
                              const spheres::box_trans &trans_spheres_obj,
                              const std::pair<size_t, size_t> sphere_number)
{
    double pi = numeric::NumericTraits<core::Real>::pi();
    for(size_t i = std::get<1>(peptide_ranges.phipsi); i != std::get<2>(peptide_ranges.phipsi); i++)
    {
        x[i] = pi*(2.0*x[i] - 1.0);
    }
    for(size_t i = 0, end = opt_vector.size(); i < end; i++)
    {
        pose.set_dof(opt_vector[i].dofid, x[i]);
    }

    x = transform::peptide_quaternion(x, opt_vector, opt_vector.size());
    x = transform::twospheres(x, opt_vector.size(), trans_spheres_obj);

    ///

    core::Vector new_position = shift_params.InitPeptidePosition;
    new_position.at(0) = x[opt_vector.size()];
    new_position.at(1) = x[opt_vector.size() + 1];
    new_position.at(2) = x[opt_vector.size() + 2];

    shift_params.FlexibleJump.reset();

    shift_params.FlexibleJump.set_rotation(shift_params.InitRm);
    pose.set_jump(pose.num_jump(), shift_params.FlexibleJump);

    numeric::Real current_distance(pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA").distance(new_position));
    shift_params.FlexibleJump.translation_along_axis(shift_params.UpstreamStub, new_position - pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA"), current_distance);
    pose.set_jump(pose.num_jump(), shift_params.FlexibleJump);

    core::Vector peptide_c_alpha_centroid(pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA"));
    core::Vector axis(x[opt_vector.size() + 3], x[opt_vector.size() + 4], x[opt_vector.size() + 5]);

    shift_params.SpinMover.rot_center(peptide_c_alpha_centroid);
    shift_params.SpinMover.spin_axis(axis);
    shift_params.SpinMover.angle_magnitude(x.back());
    shift_params.SpinMover.apply(pose);

    ///

    std::vector<double> exponential_center =
    {
        trans_spheres_obj.spheres[sphere_number.second].x,
        trans_spheres_obj.spheres[sphere_number.second].y,
        trans_spheres_obj.spheres[sphere_number.second].z
    };

    core::Vector lCA = pose.residue(complex_stuff.get_peptide_last_index()).xyz("CA");

    double rez = 0;
    double temp = 0;
    for(size_t i = 0; i != exponential_center.size(); i++)
    {
        rez += std::pow(lCA.at(i) - exponential_center[i], 2.0);
    }
    temp = std::sqrt(rez);
    if(temp < trans_spheres_obj.spheres[sphere_number.second].r)
    {
        rez = 1.1;
    }
    else
        rez = std::exp(-temp*width);
    return -rez;
}


double find_sphere_quant(std::vector<double> x,
                         std::shared_ptr<empirical_quantile::ImplicitQuantile<size_t, double>> quant,
                         core::pose::Pose& pose,
                         PoseShift &shift_params,
                         const std::vector<opt_element> &opt_vector,
                         const ComplexInfoNseq &complex_stuff,
                         double width,
                         const pepdockopt::ranges &peptide_ranges,
                         const spheres::box_trans &trans_spheres_obj,
                         size_t sp1)
{
    double pi = numeric::NumericTraits<core::Real>::pi();
    for(size_t i = std::get<1>(peptide_ranges.phipsi); i != std::get<2>(peptide_ranges.phipsi); i++)
    {
        x[i] = pi*(2.0*x[i] - 1.0);
    }
    for(size_t i = 0, end = opt_vector.size(); i < end; i++)
    {
        pose.set_dof(opt_vector[i].dofid, x[i]);
    }

    x = transform::peptide_quaternion(x, opt_vector, opt_vector.size());
    //x = transform::twospheres(x, opt_vector.size(), trans_spheres_obj);

    std::vector<double> t1 = {x[opt_vector.size()], x[opt_vector.size() + 1], x[opt_vector.size() + 2]};
    std::vector<double> s1(t1.size());
    quant->transform(t1, s1);
    x[opt_vector.size()] = s1[0];
    x[opt_vector.size() + 1] = s1[1];
    x[opt_vector.size() + 2] = s1[2];

    ///

    core::Vector new_position = shift_params.InitPeptidePosition;
    new_position.at(0) = x[opt_vector.size()];
    new_position.at(1) = x[opt_vector.size() + 1];
    new_position.at(2) = x[opt_vector.size() + 2];

    shift_params.FlexibleJump.reset();

    shift_params.FlexibleJump.set_rotation(shift_params.InitRm);
    pose.set_jump(pose.num_jump(), shift_params.FlexibleJump);

    numeric::Real current_distance(pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA").distance(new_position));
    shift_params.FlexibleJump.translation_along_axis(shift_params.UpstreamStub, new_position - pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA"), current_distance);
    pose.set_jump(pose.num_jump(), shift_params.FlexibleJump);

    core::Vector peptide_c_alpha_centroid(pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA"));
    core::Vector axis(x[opt_vector.size() + 3], x[opt_vector.size() + 4], x[opt_vector.size() + 5]);

    shift_params.SpinMover.rot_center(peptide_c_alpha_centroid);
    shift_params.SpinMover.spin_axis(axis);
    shift_params.SpinMover.angle_magnitude(x.back());
    shift_params.SpinMover.apply(pose);

    ///

    size_t sp2 = sp1 ? 0 : 1;
    std::vector<double> center1 =
    {
        trans_spheres_obj.spheres[sp1].x,
        trans_spheres_obj.spheres[sp1].y,
        trans_spheres_obj.spheres[sp1].z
    };
    std::vector<double> center2 =
    {
        trans_spheres_obj.spheres[sp2].x,
        trans_spheres_obj.spheres[sp2].y,
        trans_spheres_obj.spheres[sp2].z
    };

    core::Vector fCA = pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA");
    core::Vector lCA = pose.residue(complex_stuff.get_peptide_last_index()).xyz("CA");

    double rl1 = 0, rl2 = 0, rf1 = 0, rf2 = 0;
    double dl1 = 0, dl2 = 0, df1 = 0, df2 = 0;
    for(size_t i = 0; i != center1.size(); i++)
    {
        rl1 += std::pow(lCA.at(i) - center1[i], 2.0);
        rl2 += std::pow(lCA.at(i) - center2[i], 2.0);
        rf1 += std::pow(fCA.at(i) - center1[i], 2.0);
        rf2 += std::pow(fCA.at(i) - center2[i], 2.0);
    }
    dl1 = std::sqrt(rl1);
    dl2 = std::sqrt(rl2);
    df1 = std::sqrt(rf1);
    df2 = std::sqrt(rf2);
    double rez = 0;
    if(dl1 < trans_spheres_obj.spheres[sp1].r/2.0 && df2 < trans_spheres_obj.spheres[sp2].r/2.0)
        //if(dl1 < trans_spheres_obj.spheres[sp1].r)
    {
        rez = 1.1;
    }
    else
        rez = 0.5*std::exp(-dl1*width) + 0.5*std::exp(-df2*width);
    if(df2 > df1 || dl1 > dl2)
        rez = -1;
    return -rez;
}

double find_sphere_quant_check(std::vector<double> x,
                               std::shared_ptr<empirical_quantile::ImplicitQuantile<size_t, double>> quant,
                               core::pose::Pose& pose,
                               PoseShift &shift_params,
                               const std::vector<opt_element> &opt_vector,
                               const ComplexInfoNseq &complex_stuff,
                               const pepdockopt::ranges &peptide_ranges,
                               const spheres::box_trans &trans_spheres_obj,
                               size_t sp1)
{
    double pi = numeric::NumericTraits<core::Real>::pi();
    for(size_t i = std::get<1>(peptide_ranges.phipsi); i != std::get<2>(peptide_ranges.phipsi); i++)
    {
        x[i] = pi*(2.0*x[i] - 1.0);
    }
    for(size_t i = 0, end = opt_vector.size(); i < end; i++)
    {
        pose.set_dof(opt_vector[i].dofid, x[i]);
    }

    x = transform::peptide_quaternion(x, opt_vector, opt_vector.size());
    //x = transform::twospheres(x, opt_vector.size(), trans_spheres_obj);

    std::vector<double> t1 = {x[opt_vector.size()], x[opt_vector.size() + 1], x[opt_vector.size() + 2]};
    std::vector<double> s1(t1.size());
    quant->transform(t1, s1);
    x[opt_vector.size()] = s1[0];
    x[opt_vector.size() + 1] = s1[1];
    x[opt_vector.size() + 2] = s1[2];

    ///

    core::Vector new_position = shift_params.InitPeptidePosition;
    new_position.at(0) = x[opt_vector.size()];
    new_position.at(1) = x[opt_vector.size() + 1];
    new_position.at(2) = x[opt_vector.size() + 2];

    shift_params.FlexibleJump.reset();

    shift_params.FlexibleJump.set_rotation(shift_params.InitRm);
    pose.set_jump(pose.num_jump(), shift_params.FlexibleJump);

    numeric::Real current_distance(pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA").distance(new_position));
    shift_params.FlexibleJump.translation_along_axis(shift_params.UpstreamStub, new_position - pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA"), current_distance);
    pose.set_jump(pose.num_jump(), shift_params.FlexibleJump);

    core::Vector peptide_c_alpha_centroid(pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA"));
    core::Vector axis(x[opt_vector.size() + 3], x[opt_vector.size() + 4], x[opt_vector.size() + 5]);

    shift_params.SpinMover.rot_center(peptide_c_alpha_centroid);
    shift_params.SpinMover.spin_axis(axis);
    shift_params.SpinMover.angle_magnitude(x.back());
    shift_params.SpinMover.apply(pose);

    ///

    size_t sp2 = sp1 ? 0 : 1;
    std::vector<double> center1 =
    {
        trans_spheres_obj.spheres[sp1].x,
        trans_spheres_obj.spheres[sp1].y,
        trans_spheres_obj.spheres[sp1].z
    };
    std::vector<double> center2 =
    {
        trans_spheres_obj.spheres[sp2].x,
        trans_spheres_obj.spheres[sp2].y,
        trans_spheres_obj.spheres[sp2].z
    };

    core::Vector fCA = pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA");
    core::Vector lCA = pose.residue(complex_stuff.get_peptide_last_index()).xyz("CA");

    double rl1 = 0, rl2 = 0, rf1 = 0, rf2 = 0;
    double dl1 = 0, dl2 = 0, df1 = 0, df2 = 0;
    for(size_t i = 0; i != center1.size(); i++)
    {
        rl1 += std::pow(lCA.at(i) - center1[i], 2.0);
        rl2 += std::pow(lCA.at(i) - center2[i], 2.0);
        rf1 += std::pow(fCA.at(i) - center1[i], 2.0);
        rf2 += std::pow(fCA.at(i) - center2[i], 2.0);
    }
    dl1 = std::sqrt(rl1);
    dl2 = std::sqrt(rl2);
    df1 = std::sqrt(rf1);
    df2 = std::sqrt(rf2);
    double rez = 0;
    if(dl1 < trans_spheres_obj.spheres[sp1].r && df2 < trans_spheres_obj.spheres[sp2].r)
        //if(dl1 < trans_spheres_obj.spheres[sp1].r)
    {
        rez = 1.1;
    }
    else
        rez = 0;
    if(df2 > df1 || dl1 > dl2)
        rez = 0;
    return rez;
}

double find_sphere_quant_check_omp(std::vector<double> x,
                                   int th_id,
                                   std::shared_ptr<empirical_quantile::ImplicitQuantile<size_t, double>> quant,
                                   std::vector<core::pose::Pose>& pose,
                                   std::vector<PoseShift> &shift_params,
                                   const std::vector<opt_element> &opt_vector,
                                   const ComplexInfoNseq &complex_stuff,
                                   const pepdockopt::ranges &peptide_ranges,
                                   const spheres::box_trans &trans_spheres_obj,
                                   size_t sp1)
{
    double pi = numeric::NumericTraits<core::Real>::pi();
    for(size_t i = std::get<1>(peptide_ranges.phipsi); i != std::get<2>(peptide_ranges.phipsi); i++)
    {
        x[i] = pi*(2.0*x[i] - 1.0);
    }
    for(size_t i = 0, end = opt_vector.size(); i < end; i++)
    {
        pose[th_id].set_dof(opt_vector[i].dofid, x[i]);
    }

    x = transform::peptide_quaternion(x, opt_vector, opt_vector.size());
    //x = transform::twospheres(x, opt_vector.size(), trans_spheres_obj);

    std::vector<double> t1 = {x[opt_vector.size()], x[opt_vector.size() + 1], x[opt_vector.size() + 2]};
    std::vector<double> s1(t1.size());
    quant->transform(t1, s1);
    x[opt_vector.size()] = s1[0];
    x[opt_vector.size() + 1] = s1[1];
    x[opt_vector.size() + 2] = s1[2];

    ///

    core::Vector new_position = shift_params[th_id].InitPeptidePosition;
    new_position.at(0) = x[opt_vector.size()];
    new_position.at(1) = x[opt_vector.size() + 1];
    new_position.at(2) = x[opt_vector.size() + 2];

    shift_params[th_id].FlexibleJump.reset();

    shift_params[th_id].FlexibleJump.set_rotation(shift_params[th_id].InitRm);
    pose[th_id].set_jump(pose[th_id].num_jump(), shift_params[th_id].FlexibleJump);

    numeric::Real current_distance(pose[th_id].residue(complex_stuff.get_peptide_first_index()).xyz("CA").distance(new_position));
    shift_params[th_id].FlexibleJump.translation_along_axis(shift_params[th_id].UpstreamStub, new_position - pose[th_id].residue(complex_stuff.get_peptide_first_index()).xyz("CA"), current_distance);
    pose[th_id].set_jump(pose[th_id].num_jump(), shift_params[th_id].FlexibleJump);

    core::Vector peptide_c_alpha_centroid(pose[th_id].residue(complex_stuff.get_peptide_first_index()).xyz("CA"));
    core::Vector axis(x[opt_vector.size() + 3], x[opt_vector.size() + 4], x[opt_vector.size() + 5]);

    shift_params[th_id].SpinMover.rot_center(peptide_c_alpha_centroid);
    shift_params[th_id].SpinMover.spin_axis(axis);
    shift_params[th_id].SpinMover.angle_magnitude(x.back());
    shift_params[th_id].SpinMover.apply(pose[th_id]);

    ///

    size_t sp2 = sp1 ? 0 : 1;
    std::vector<double> center1 =
    {
        trans_spheres_obj.spheres[sp1].x,
        trans_spheres_obj.spheres[sp1].y,
        trans_spheres_obj.spheres[sp1].z
    };
    std::vector<double> center2 =
    {
        trans_spheres_obj.spheres[sp2].x,
        trans_spheres_obj.spheres[sp2].y,
        trans_spheres_obj.spheres[sp2].z
    };

    core::Vector fCA = pose[th_id].residue(complex_stuff.get_peptide_first_index()).xyz("CA");
    core::Vector lCA = pose[th_id].residue(complex_stuff.get_peptide_last_index()).xyz("CA");

    double rl1 = 0, rl2 = 0, rf1 = 0, rf2 = 0;
    double dl1 = 0, dl2 = 0, df1 = 0, df2 = 0;
    for(size_t i = 0; i != center1.size(); i++)
    {
        rl1 += std::pow(lCA.at(i) - center1[i], 2.0);
        rl2 += std::pow(lCA.at(i) - center2[i], 2.0);
        rf1 += std::pow(fCA.at(i) - center1[i], 2.0);
        rf2 += std::pow(fCA.at(i) - center2[i], 2.0);
    }
    dl1 = std::sqrt(rl1);
    dl2 = std::sqrt(rl2);
    df1 = std::sqrt(rf1);
    df2 = std::sqrt(rf2);
    double rez = 0;
    if(dl1 < trans_spheres_obj.spheres[sp1].r && df2 < trans_spheres_obj.spheres[sp2].r)
        //if(dl1 < trans_spheres_obj.spheres[sp1].r)
    {
        rez = 1.1;
    }
    else
        rez = 0;
    if(df2 > df1 || dl1 > dl2)
        rez = 0;
    return rez;
}


template <typename T, typename V>
struct cell
{
    std::vector<T> dot;
    V value;

    cell() {}
    cell(std::vector<T> _dot, V _value):dot(_dot),value(_value) {}
    cell(const cell &a):dot(a.dot), value(a.value) {}
    bool operator==(const cell& a) const
    {
        return dot == a.dot;
    }
    bool operator<(const cell& a) const
    {
        return dot < a.dot;
    }
};

void boundaryFill4MultipleGrids_omp_trie_only(std::vector<std::vector<double>>& grids,
        std::vector<std::uint8_t> &start,
        std::shared_ptr<trie_based::TrieBased<trie_based::NodeCount<std::uint8_t>,std::uint8_t>> samples,
        std::vector<double> &dx,
        int total_threads,
        std::function<double(std::vector<double>, int)> f)
{
    /// generate all permutations
//    std::vector<int> ind;
//    ind.push_back(0);
//    ind.push_back(grids.size() - 4);
//    ind.push_back(grids.size() - 3);
//    ind.push_back(grids.size() - 2);
//    ind.push_back(grids.size() - 1);
//    std::vector<std::vector<char>> variable_values(ind.size(), std::vector<char>(3));
//    for(size_t i = 0; i != variable_values.size(); i++)
//    {
//        variable_values[i][0] = -1;
//        variable_values[i][1] = 0;
//        variable_values[i][2] = 1;
//    }
//    std::vector<std::vector<char>> permut = iterate(variable_values);

    int dim = grids.size();
    std::vector<double> dot(grids.size());

    trie_based::TrieBased<trie_based::Node<std::uint8_t>,std::uint8_t> visited;
    visited.set_dimension(dim);
    trie_based::TrieBased<trie_based::Node<std::uint8_t>,std::uint8_t> points;
    points.set_dimension(dim);
    trie_based::TrieBased<trie_based::Node<std::uint8_t>,std::uint8_t> not_coumputed;
    not_coumputed.set_dimension(dim);

    for(size_t i = 0; i != dot.size(); i++)
    {
        dot[i] = grids[i][start[i]] + dx[i];
    }
    if(f(dot, 0) > 1.05)
        points.insert(start);
    else
        return;
    
    while(!points.empty() || !not_coumputed.empty())
    {
        while(!points.empty())
        {
            auto point = points.get_and_remove_last();
            if(visited.search(point) || samples->search(point))
            {
                //points.pop_front();
                continue;
            }
            //visited.insert(point);
            samples->insert(point);

            auto init_point = point;
            for(size_t i = 0; i != point.size(); i++)
            {
                point = init_point;
                //point[i] = point[i] + 1;
                int new_val = point[i] + 1;

                if(new_val < 0)
                {
                    new_val = grids[i].size() - 1;
                }
                if(new_val > grids[i].size() - 1)
                {
                    new_val = 0;
                }
                point[i] = new_val;
//                if(point[i] < 0 || point[i] > grids[i].size() - 1)
//                    continue;
                if(!visited.search(point) && !samples->search(point))
                {
                    not_coumputed.insert(point);
                }
            }
            for(size_t i = 0; i != point.size(); i++)
            {
                point = init_point;
                //point[i] = point[i] - 1;
                int new_val = point[i] - 1;

                if(new_val < 0)
                {
                    new_val = grids[i].size() - 1;
                }
                if(new_val > grids[i].size() - 1)
                {
                    new_val = 0;
                }
                point[i] = new_val;

//                if(point[i] < 0 || point[i] > grids[i].size() - 1)
//                    continue;

                if(!visited.search(point) && !samples->search(point))
                {
                    not_coumputed.insert(point);
                }
            }

            ///Moore for 4 dots only
//            for(size_t i = 0; i != permut.size(); i++)
//            {
//                point = init_point;
//                bool flag = true;
//                for(size_t j = 0; j != variable_values.size(); j++)
//                {
//                    point[ind[j]] = point[ind[j]] + permut[i][j];
////                    if(point[ind[j]] < 0 || point[ind[j]] > grids[ind[j]].size() - 1)
////                    {
////                        flag = false;
////                        break;
////                    }
//
//                    if(point[ind[j]] < 0)
//                    {
//                        point[ind[j]] = grids[ind[j]].size() - 1;
//                    }
//                    if(point[ind[j]] > grids[ind[j]].size() - 1)
//                    {
//                        point[ind[j]] = 0;
//                    }
//                }
//                if(!flag)
//                    continue;
//
//                if(!visited.search(point) && !samples->search(point) && !not_coumputed.search(point))
//                {
//                    not_coumputed.insert(point);
//                }
//            }
        }
        size_t number_to_points = 0;
        while(!not_coumputed.empty())
        {
            std::vector<cell<std::uint8_t,char>> to_compute;
            for(size_t i = 0; i != total_threads*100000; i++)
            {
                if(not_coumputed.empty())
                    break;

                auto point = not_coumputed.get_and_remove_last();
                to_compute.push_back(cell<std::uint8_t,char>(point, -1));
            }
            int omp_size = to_compute.size();
            while(omp_size%total_threads)
            {
                omp_size--;
            }

            //omp_size = 0;

            int th_id;
            #pragma omp parallel private(th_id)
            {
                th_id = omp_get_thread_num();
                if(th_id < total_threads)
                {
                    for(int i = th_id * omp_size / total_threads; i < (th_id + 1) * omp_size / total_threads; ++i)
                    {
                        std::vector<double> values(dim);
                        for(size_t j = 0; j != values.size(); j++)
                        {
                            values[j] = grids[j][to_compute[i].dot[j]] + dx[j];
                        }
                        to_compute[i].value = f(values, th_id);
                    }
                }
            }

            for(size_t i = omp_size; i != to_compute.size(); i++)
            {
                std::vector<double> values(dim);
                for(size_t j = 0; j != values.size(); j++)
                {
                    values[j] = grids[j][to_compute[i].dot[j]] + dx[j];
                }
                to_compute[i].value = f(values, 0);
            }

            for(size_t i = 0; i != to_compute.size(); i++)
            {
                if(to_compute[i].value == 1)
                {
                    points.insert(to_compute[i].dot);
                    ++number_to_points;
                }
                else if(to_compute[i].value == 0)
                {
                    visited.insert(to_compute[i].dot);
                }
            }
        }
        std::cout << number_to_points << std::endl;
    }
    visited.remove_tree();
    points.remove_tree();
    not_coumputed.remove_tree();
}

void boundaryFill4MultipleGrids_omp(std::vector<std::vector<double>>& grids,
                                    std::vector<std::vector<char>> &pp,
                                    std::set<std::vector<char>> &visited,
                                    std::vector<std::vector<char>> &samples,
                                    std::vector<double> &dx,
                                    std::function<double(std::vector<double>, int)> f,
                                    int total_threads)
{
    /// generate all permutations
    std::vector<std::vector<char>> variable_values(pp.back().size() - 7, std::vector<char>(3));
    for(size_t i = 0; i != variable_values.size(); i++)
    {
        variable_values[i][0] = -1;
        variable_values[i][1] = 0;
        variable_values[i][2] = 1;
    }
    std::vector<std::vector<char>> permut = iterate(variable_values);


    std::set<cell<char,double>> computed;
    std::vector<double> dot(pp.back().size());
    std::vector<std::vector<double>> dots;
    std::deque<cell<char,double>> points;
    for(size_t i = 0; i != dot.size(); i++)
    {
        dot[i] = grids[i][pp.back()[i]] + dx[i];
    }
    points.push_back(cell<char,double>(pp.back(),-1.0));

    std::vector<size_t> index;

    while(!points.empty())
    {
        //auto t = points.front();
        //points.pop_front();

        auto visit_it = visited.find(points.front().dot);
        if(!(visit_it == visited.end()))
        {
            points.pop_front();
            continue;
        }
        visited.insert(points.front().dot);

        if(points.front().value < 0)
        {
            index.clear();
            dots.clear();

            for(size_t i = 0; i < points.size(); i++)
            {
                auto comp_it = computed.find(points[i]);
                if(comp_it == computed.end())
                {
                    index.push_back(i);
                }
                else
                {
                    points[i].value = comp_it->value;
                }
            }

            std::cout << index.size() << std::endl;
            dots.resize(index.size());

            for(size_t i = 0; i < dots.size(); i++)
            {
                for(size_t j = 0; j != dot.size(); j++)
                {
                    dot[j] = grids[j][points[index[i]].dot[j]] + dx[j];
                }
                dots[i] = dot;
                //points[index[i]].value = pdf(dot);
            }
            /*#pragma omp parallel for
            for(size_t i = 0; i < index.size(); i++)
            {
                points[index[i]].value = Quadratic_ValueSimple3(dots[i], );
            }*/

            int omp_size = index.size();
            while(omp_size%total_threads)
            {
                omp_size--;
            }

            int th_id = 0;
            #pragma omp parallel private(th_id)
            {
                th_id = omp_get_thread_num();
                if(th_id < total_threads)
                {
                    for(int i = th_id * omp_size / total_threads; i < (th_id + 1) * omp_size / total_threads; ++i)
                    {
                        points[index[i]].value = f(dots[i], th_id);
                    }
                }
            }

            for(size_t i = omp_size; i != index.size(); i++)
            {
                points[index[i]].value = f(dots[i], th_id);
            }

            for(size_t i = 0; i < index.size(); i++)
            {
                computed.insert(points[index[i]]);
                //fe_count++;
            }
        }

        if(points.front().value > 1.0)
        {

            if(!(samples.size()%100))
            {
                std::cout << "Samples   " << samples.size() << std::endl;
                std::cout << "Points    " << points.size() << std::endl;
                std::cout << "Computed  " << computed.size() << std::endl;
                std::cout << "Visited   " << visited.size() << std::endl;
                std::cout << std::endl;
            }


            std::vector<char> point = points.front().dot;
            samples.push_back(point);

            for(size_t i = 0; i != point.size(); i++)
            {
                point = points.front().dot;
                point[i] = point[i] + 1;

                if(point[i] < 0 || point[i] > grids[i].size() - 1)
                    continue;

                visit_it = visited.find(point);
                if(visit_it == visited.end())
                {
                    points.push_back(cell<char,double>(point, -1.0));
                }

                //points.push_back(cell<char,float>(point, -1.0));
            }
            for(size_t i = 0; i != point.size(); i++)
            {
                point = points.front().dot;
                point[i] = point[i] - 1;

                if(point[i] < 0 || point[i] > grids[i].size() - 1)
                    continue;

                visit_it = visited.find(point);
                if(visit_it == visited.end())
                {
                    points.push_back(cell<char,double>(point, -1.0));
                }

                //points.push_back(cell<char,float>(point, -1.0));
            }

            /*for(size_t i = 0; i != permut.size(); i++)
            {
                point = points.front().dot;
                bool flag = true;
                for(size_t j = 0; j != variable_values.size(); j++)
                {
                    point[j] = point[j] + permut[i][j];
                    if(point[j] < 0 || point[j] > grids[j].size() - 1)
                    {
                        flag = false;
                        break;
                    }
                }
                if(!flag)
                    continue;

                visit_it = visited.find(point);
                if(visit_it == visited.end())
                {
                    points.push_back(cell<char,double>(point, -1.0));
                }
            }*/
        }
        points.pop_front();

        //if(samples.size() > 100000)
        //{
        //    break;
        //}
    }
}

/*
class mc
{
public:
    mc(std::vector<double> _lb, std::vector<double> _ub, size_t _fe_max)
    {
        lb = _lb;
        ub = _ub;
        dim = lb.size();
        fe_max = _fe_max;
    }
    void set_generator_seed(int _seed)
    {
        generator.seed(_seed);
    }
    void run(core::pose::Pose& pose,
             protocols::rigid::RigidBodyDeterministicSpinMover& SpinMover,
             std::vector<opt_element> &opt_vector,
             core::kinematics::Jump &flexible_jump,
             core::kinematics::Stub &upstream_stub,
             core::Vector &init_peptide_position,
             core::kinematics::RT::Matrix &init_rm,
             ComplexInfoNseq &complex_stuff,
             double width,
             pepdockopt::ranges peptide_ranges,
             spheres::box_trans trans_spheres_obj,
             std::pair<size_t, size_t> sphere_number)
    {
        best = std::numeric_limits<double>::max();
        std::uniform_real_distribution<double> uniform_real01(0.0, 1.0);
        std::vector<double> x(dim);
        std::vector<double> xx(dim);
        for(size_t i = 0; i != dim; i++)
        {
            std::uniform_real_distribution<double> uniform_real_dist(lb[i], ub[i]);
            x[i] = uniform_real_dist(generator);
        }
        double F = Quadratic_ValueSimple2(x,
                                          pose,
                                          SpinMover,
                                          opt_vector,
                                          flexible_jump,
                                          upstream_stub,
                                          init_peptide_position,
                                          init_rm,
                                          complex_stuff,
                                          width,
                                          peptide_ranges,
                                          trans_spheres_obj,
                                          sphere_number);

        double theta = 1.0;
        double dtheta = 0.999;
        std::vector<double> dx(lb.size());
        for(size_t i = 0; i != dx.size(); i++)
        {
            dx[i] = (ub[i] - lb[i])/1000.0;
        }
        for(size_t fe_count = 0; fe_count != fe_max; fe_count++)
        {
            for(size_t i = 0; i != dim; i++)
            {
                std::uniform_real_distribution<double> uniform_real_dist(-dx[i], dx[i]);
                xx[i] = x[i] + uniform_real_dist(generator);
                if(xx[i] < lb[i])
                    xx[i] = lb[i];
                if(xx[i] > ub[i])
                    xx[i] = ub[i];
            }
            double FF = Quadratic_ValueSimple2(xx,
                                               pose,
                                               SpinMover,
                                               opt_vector,
                                               flexible_jump,
                                               upstream_stub,
                                               init_peptide_position,
                                               init_rm,
                                               complex_stuff,
                                               width,
                                               peptide_ranges,
                                               trans_spheres_obj,
                                               sphere_number);
            if(best > FF)
            {
                best = FF;
                std::cout << best << std::endl;
                best_vector = xx;
            }
            if(FF < F || uniform_real01(generator) < std::exp((F - FF) / theta))
            {
                std::swap(x, xx);
                F = FF;
            }
            theta *= dtheta;
            for(size_t i = 0; i != dim; i++)
                dx[i] -= ((ub[i] - lb[i])/1000.0) / fe_max;
        }
    }
    void set_generator(std::mt19937_64& _generator)
    {
        generator = _generator;
    }
    std::vector<double> get_best_vector()
    {
        return best_vector;
    }
    double get_best()
    {
        return best;
    }
    double (*FitnessFunction)(const std::vector<double> &x) = nullptr;
private:
    size_t dim;

    size_t fe_max;
    double best;
    std::vector<double> best_vector;
    std::vector<double> lb;
    std::vector<double> ub;

    std::mt19937_64 generator;
};*/


std::vector<double> PepDockOpt::get_position(std::vector<double> _lb,
        std::vector<double> _ub,
        double width,
        spheres::box_trans trans_sp,
        std::pair<size_t, size_t> spheres_number)
{
    std::vector<double> start(_lb.size());
    std::vector<double> start_temp(start.size());

    double best = 0;
    bool is_in_sphere = false;
    bool reduce_step = true;
    std::mt19937_64 generator_;
    generator_.seed(1);
    for(size_t j = 0; j != start.size(); j++)
    {
        std::cout << _lb[j] << '\t' << _ub[j] << std::endl;
        std::uniform_real_distribution<double> distribution(_lb[j], _ub[j]);
        start[j] = distribution(generator_);
    }

    /*mc obj1(_lb, _ub, 1e5);
    obj1.set_generator(generator_);
    do
    {
        obj1.run(pose.front(),
                 pose_shift.front().SpinMover,
                 opt_vector,
                 pose_shift.front().FlexibleJump,
                 pose_shift.front().UpstreamStub,
                 pose_shift.front().InitPeptidePosition,
                 pose_shift.front().InitRm,
                 param_list,
                 width,
                 peptide_ranges,
                 trans_spheres_obj,
                 spheres_number);
        //std::cout << obj1.get_best() << std::endl;
        start = obj1.get_best_vector();
    }
    while(obj1.get_best() > -1.05);*/

    pepdockopt::cso::CSO obj1(1000, start.size(), 0.01, 1e6, generator_);
    obj1.FitnessFunction = std::bind(&Quadratic_ValueSimple2,
                                     std::placeholders::_1,
                                     pose.front(),
                                     pose_shift.front(),
                                     opt_vector,
                                     param_list,
                                     width,
                                     peptide_ranges,
                                     trans_sp,
                                     spheres_number);
    do
    {
        obj1.set_bounds(lb, ub);
        obj1.init();
        obj1.optimization();
        std::cout << obj1.get_best() << std::endl;
        start = obj1.get_best_vector();
        //std::cin.get();
    }
    while(obj1.get_best() > -1.05);

    //std::function<double(std::vector<double>)> FitnessFunction = std::bind(&lol, std::placeholders::_1);
    /*std::function<double(std::vector<double>)> FitnessFunction = std::bind(&lol2,
            std::placeholders::_1,
            pose.front(),
            pose_shift);*/

    return start;
}

void PepDockOpt::set_score_function()
{
    score_func.resize(threads_number);
    core::scoring::ScoreFunctionOP score_fn = core::scoring::get_score_function();
    std::cout << "score function: " << score_fn->get_name() << std::endl;
    for(auto &i : score_func)
        i = score_fn;
}

void PepDockOpt::init(size_t _threads_number)
{
    threads_number = _threads_number;
    if(omp_get_num_procs() < threads_number)
    {
        threads_number = omp_get_num_procs();
    }
    std::cout << "number_of_threads: " << threads_number << std::endl;

    pose.resize(threads_number);
    set_score_function();

    /// load input structures
    core::pose::Pose input;
    if(basic::options::option[basic::options::OptionKeys::in::file::s].user())
    {
        std::string fname = basic::options::option[basic::options::OptionKeys::in::file::s].value_string();
        std::cout << "input pose " << fname << std::endl;
        core::import_pose::pose_from_file(input, fname);
        input.dump_pdb("output/pdb/input.pdb");

        spheres_fname = fname;
        spheres_fname = spheres_fname.substr(0, fname.length()-4);
        spheres_fname += ".dat";
        std::cout << "spheres " << spheres_fname << std::endl;
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

    for(auto &i : pose)
        i = param_list.get_pose_complex(); //param_list.get_peptide();

    pose_shift.resize(threads_number);
    for(int i = 0; i != threads_number; i++)
    {
        pose_shift[i].SpinMover.rb_jump(pose[i].num_jump());
        pose_shift[i].FlexibleJump = pose[i].jump(pose[i].num_jump());
        pose_shift[i].InitPeptidePosition = pose[i].residue(param_list.get_peptide_first_index()).xyz("CA");
        pose_shift[i].InitRm = pose_shift[i].FlexibleJump.get_rotation();
        pose_shift[i].UpstreamStub = pose[i].conformation().upstream_jump_stub(pose[i].num_jump());
    }
}
void PepDockOpt::set_opt()
{
    std::vector<core::Size> protein_full_index, peptide_full_index;
    //for( size_t i = param_list.get_protein_first_index(); i <= param_list.get_protein_last_index(); i++ )
    //    protein_full_index.push_back( i );
    for(size_t i = param_list.get_peptide_first_index(); i <= param_list.get_peptide_last_index(); i++)
        peptide_full_index.push_back(i);

    std::vector<pepdockopt::opt_element> peptide_mc_p_info = pepdockopt::get_phi_psi_with_info(pose.front(), peptide_full_index, 2, true);
    //std::vector<opt_element> peptide_mc_o_info = get_omega_with_info(pose, peptide_full_index, 2, true);
    //std::vector<opt_element> peptide_all_chi_info = get_peptide_all_chi_dof_with_info(pose, peptide_full_index, 2, true);
    peptide_mc_p_info.pop_back();

    peptide_ranges.phipsi = std::make_tuple(true, opt_vector.size(), opt_vector.size() + peptide_mc_p_info.size());
    opt_vector.insert(opt_vector.end(), peptide_mc_p_info.begin(), peptide_mc_p_info.end());

    std::vector<core::Size> protein_cm_index;
    trans_spheres_obj.load_data(spheres_fname, 100);
    //trans_spheres_obj.load_data("input/pdb/DR1-RA_0001.dat", 100);
    //trans_spheres_obj.load_data("input/pdb/1JWG_0001.dat", 100);
//    trans_spheres_obj1.load_data("input/pdb/1JWG_0001.dat", 100, false); //1JWG_0001_3
//    trans_spheres_obj2.load_data("input/pdb/1JWG_0001.dat", 100, true); //1JWG_0001_3

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


    //
    lb.clear();
    ub.clear();
    size_t k = 0;
    for(size_t i = 0; i < peptide_mc_p_info.size(); i++, k++)
    {
        lb.push_back(/*peptide_mc_p_info[i].bounds.first*/0);
        ub.push_back(/*peptide_mc_p_info[i].bounds.second*/1);
    }
    //for(size_t i = 0; i < peptide_mc_o_info.size(); i++, k++)
    //{
    //    lb.push_back(3.0 * numeric::NumericTraits<core::Real>::pi() - omega_delta);
    //    ub.push_back(3.0 * numeric::NumericTraits<core::Real>::pi() + omega_delta);
    //}

//    for(size_t i = 0; i < peptide_all_chi_info.size(); i++, k++)
//    {
//        lb.push_back(peptide_all_chi_info[i].bounds.first);
//        ub.push_back(peptide_all_chi_info[i].bounds.second);
//    }

//    for(size_t i = 0; i < protein_all_chi_info.size(); i++, k++)
//    {
//        lb.push_back(protein_all_chi_info[i].bounds.first);
//        ub.push_back(protein_all_chi_info[i].bounds.second);
//    }


    for(size_t i = 0; i != 3; i++, k++)
    {
        lb.push_back(0.0);
        ub.push_back(1.0);
    }
    for(size_t i = 0; i != 3; i++, k++)
    {
        lb.push_back(0.0);
        ub.push_back(1.0);
    }
    lb.push_back(0.0);
    ub.push_back(1.0);
    k++;

    if(k != opt_vector.size() + 7 || lb.size() - 7 != opt_vector.size() || ub.size() - 7 != opt_vector.size())
    {
        std::cout << "fail\n";
        return;
    }

    std::cout << opt_vector.size() << std::endl;

    int dim = opt_vector.size() + 7;
    std::cout << "dim = " << dim << std::endl;
}

void PepDockOpt::set_grid()
{
    std::vector<std::vector<double>> best_samples;

    gridN.resize(lb.size());
    /*gridN[0] = 36;//psi
    gridN[1] = 36;//2
    gridN[2] = 36;//2
    gridN[3] = 36;//3
    gridN[4] = 36;//3
    gridN[5] = 36;//4
    gridN[6] = 36;//4

    gridN[7] = 100;
    gridN[8] = 100;
    gridN[9] = 100;
    gridN[10] = 40;
    gridN[11] = 40;
    gridN[12] = 40;
    gridN[13] = 100;*/

    for(size_t i = 0; i != lb.size(); i++)
    {
        gridN[i] = 18;
    }
    gridN[gridN.size() - 7] = 10;
    gridN[gridN.size() - 6] = 10;
    gridN[gridN.size() - 5] = 10;
    gridN[gridN.size() - 4] = 60;
    gridN[gridN.size() - 3] = 60;
    gridN[gridN.size() - 2] = 60;
    gridN[gridN.size() - 1] = 180;

    grids.resize(gridN.size());
    dx.resize(gridN.size());

    for(size_t i = 0; i != grids.size(); i++)
    {
        double grid_N = gridN[i];
        std::vector<double> grid(grid_N);
        double startp = 0;
        double endp = 1;
        double es = endp - startp;
        for(size_t j = 0; j != grid.size(); j++)
        {
            grid[j] = startp + j*es/grid_N;
        }
        grids[i] = grid;
        dx[i] = es/(grid_N*2);
    }

//    std::cout << "grids:" << std::endl;
//    for(size_t i = 0; i != grids.size(); i++)
//    {
//        for(size_t j = 0; j != grids[i].size(); j++)
//        {
//            std::cout << grids[i][j] << ' ';
//        }
//        std::cout << std::endl;
//        for(size_t j = 0; j != grids[i].size(); j++)
//        {
//            std::cout << grids[i][j] + dx[i] << ' ';
//        }
//        std::cout << std::endl;
//        std::cout << std::endl;
//    }
//    std::cout << std::endl;


    structures_triebased = std::make_shared<frag_type>();
    structures_quant = std::make_shared<empirical_quantile::ImplicitQuantile<std::uint8_t, double>>(
                           std::vector<double>(lb.size(), 0.0),
                           std::vector<double>(ub.size(), 1.0),
                           gridN);
    structures_quant->set_sample_shared(structures_triebased);
}
void PepDockOpt::set_quantile1()
{
    std::vector<char> startdot(start1.size());
    for(size_t i = 0; i != startdot.size(); i++)
    {
        auto pos1 = std::lower_bound(grids[i].begin(), grids[i].end(), start1[i]);
        startdot[i] = std::distance(grids[i].begin(), pos1) - 1;
    }
    ///
    std::vector<double> tt(startdot.size());
    for(size_t i = 0; i != tt.size(); i++)
    {
        tt[i] = grids[i][startdot[i]] + dx[i];
    }
    if(find_sphere_quant_check(tt,
                               two_spheres_quant,
                               pose.front(),
                               pose_shift.front(),
                               opt_vector,
                               param_list,
                               peptide_ranges,
                               trans_spheres_obj,
                               1) < 1.0)
    {
        core::pose::PoseOP pep = pose[0].split_by_chain(pose[0].num_chains());
        //to_allatom(*pep);
        pep->dump_pdb("output/pdb/missed.pdb");
        std::cout << "miss the point!" << std::endl;

        /// finding the dot
        std::vector<std::vector<char>> pp;
        pp.push_back(startdot);
        std::set<std::vector<char>> vv;

        std::function<double(std::vector<double>)> f = std::bind(&find_sphere_quant_check,
                std::placeholders::_1,
                two_spheres_quant,
                std::ref(pose.front()),
                std::ref(pose_shift.front()),
                std::ref(opt_vector),
                std::ref(param_list),
                std::ref(peptide_ranges),
                std::ref(trans_spheres_obj),
                1);

        if(startdot.size() <= 15)
        {
            startdot = findTheDot_MultipleGrids_Moore(grids, pp, vv, dx, f);
        }
        else
        {
            startdot = findTheDot_MultipleGrids(grids, pp, vv, dx, f);
        }
    }

    for(size_t j = 0; j != startdot.size(); j++)
    {
        start1[j] = grids[j][startdot[j]] + dx[j];
    }
//    best_samples.push_back(start);
//
//    std::vector<std::vector<char>> points;
//    points.push_back(startdot);
//
//    std::set<std::vector<char>> visited;
//    std::vector<std::vector<char> > samples;

    //boundaryFill4MultipleGrids(grids, points, visited, samples, dx);

    std::function<double(std::vector<double>, int)> ff = std::bind(&find_sphere_quant_check_omp,
            std::placeholders::_1,
            std::placeholders::_2,
            two_spheres_quant,
            std::ref(pose),
            std::ref(pose_shift),
            std::ref(opt_vector),
            std::ref(param_list),
            std::ref(peptide_ranges),
            std::ref(trans_spheres_obj),
            1);

    //boundaryFill4MultipleGrids_omp(grids, points, visited, samples, dx, ff, threads_number);

    std::vector<std::uint8_t> st(startdot.size());
    for(size_t i = 0; i != st.size(); i++)
    {
        st[i] = startdot[i];
    }

    boundaryFill4MultipleGrids_omp_trie_only(grids, st, structures_triebased, dx, threads_number, ff);

    std::cout << "samples size " << structures_triebased->get_total_count() << std::endl;


    // quant
//    std::vector<std::uint8_t> st(startdot.size());
//    for(size_t i = 0; i != st.size(); i++)
//    {
//        st[i] = startdot[i];
//    }
//    structures_triebased->insert(st);
}

void PepDockOpt::set_quantile2()
{
    std::vector<char> startdot(start2.size());
    for(size_t i = 0; i != startdot.size(); i++)
    {
        auto pos1 = std::lower_bound(grids[i].begin(), grids[i].end(), start2[i]);
        startdot[i] = std::distance(grids[i].begin(), pos1) - 1;
    }
    ///
    std::vector<double> tt(startdot.size());
    for(size_t i = 0; i != tt.size(); i++)
    {
        tt[i] = grids[i][startdot[i]] + dx[i];
    }

    if(find_sphere_quant_check(tt,
                               two_spheres_quant,
                               pose.front(),
                               pose_shift.front(),
                               opt_vector,
                               param_list,
                               peptide_ranges,
                               trans_spheres_obj,
                               0) < 1.0)
    {
        core::pose::PoseOP pep = pose[0].split_by_chain(pose[0].num_chains());
        //to_allatom(*pep);
        pep->dump_pdb("output/pdb/missed.pdb");
        std::cout << "miss the point!" << std::endl;

        /// finding the dot
        std::vector<std::vector<char>> pp;
        pp.push_back(startdot);
        std::set<std::vector<char>> vv;

        std::function<double(std::vector<double>)> f = std::bind(&find_sphere_quant_check,
                std::placeholders::_1,
                two_spheres_quant,
                std::ref(pose.front()),
                std::ref(pose_shift.front()),
                std::ref(opt_vector),
                std::ref(param_list),
                std::ref(peptide_ranges),
                std::ref(trans_spheres_obj),
                0);

        if(startdot.size() <= 15)
        {
            startdot = findTheDot_MultipleGrids_Moore(grids, pp, vv, dx, f);
        }
        else
        {
            startdot = findTheDot_MultipleGrids(grids, pp, vv, dx, f);
        }
        //std::cout << "the point!" << std::endl;
    }


    for(size_t j = 0; j != startdot.size(); j++)
    {
        start2[j] = grids[j][startdot[j]] + dx[j];
    }
//    best_samples.push_back(start);
//
//    std::vector<std::vector<char>> points;
//    points.push_back(startdot);
//
//    std::set<std::vector<char>> visited;
//    std::vector<std::vector<char> > samples;

    //boundaryFill4MultipleGrids(grids, points, visited, samples, dx);

    /*  std::function<double(std::vector<double>, int)> ff = std::bind(&Quadratic_ValueSimple3_omp,
              std::placeholders::_1,
              std::placeholders::_2,
              std::ref(pose),
              std::ref(pose_shift),
              std::ref(opt_vector),
              std::ref(param_list),
              std::ref(peptide_ranges),
              std::ref(trans_spheres_obj),
              spheres_number);

      //boundaryFill4MultipleGrids_omp(grids, points, visited, samples, dx, ff, threads_number);


      std::vector<std::uint8_t> st(startdot.size());
      for(size_t i = 0; i != st.size(); i++)
      {
          st[i] = startdot[i];
      }

      structures_triebased = std::make_shared<frag_type>();
      //structures_triebased->insert(st);
      boundaryFill4MultipleGrids_omp_trie_only(grids, st, structures_triebased, dx, threads_number, ff);

      std::cout << "samples size " << structures_triebased->get_total_count() << std::endl;*/

    // quant
    std::vector<std::uint8_t> st(startdot.size());
    for(size_t i = 0; i != st.size(); i++)
    {
        st[i] = startdot[i];
    }
    structures_triebased->insert(st);
}
void PepDockOpt::check()
{
    structures_triebased->fill_tree_count();
    // check
    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<double> ureal01(0.0,1.0);

    std::vector<std::vector<double> > values01;
    std::vector<std::vector<double> > sampled;
    std::vector<double> temp1(gridN.size());
    std::vector<double> temp2(temp1.size());
    for(size_t i = 0; i != 100; ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }
    for(size_t j = 0; j != values01.size(); j++)
        structures_quant->transform(values01[j], sampled[j]);

    std::function<double(std::vector<double>)> fff = std::bind(&find_sphere_quant_check,
            std::placeholders::_1,
            two_spheres_quant,
            std::ref(pose.front()),
            std::ref(pose_shift.front()),
            std::ref(opt_vector),
            std::ref(param_list),
            std::ref(peptide_ranges),
            std::ref(trans_spheres_obj),
            1);

    for(size_t j = 0; j != values01.size(); j++)
    {
        fff(sampled[j]);
//        Quadratic_ValueSimple3(sampled[j], pose.front(),
//                               pose_shift.front(),
//                               opt_vector,
//                               param_list,
//                               peptide_ranges,
//                               trans_spheres_obj,
//                               spheres_number);
        core::pose::PoseOP pep = pose[0].split_by_chain(pose[0].num_chains());
        //to_allatom(*pep);
        pep->dump_pdb("output/pdb/" + std::to_string(int(1e4) + j) + ".pdb");
    }
}

void PepDockOpt::sphere_quant()
{
    //std::vector<size_t> grid_n = {50, 5, 5};
    std::vector<size_t> grid_n = {500, 50, 50};
    std::vector<double> llb = {trans_spheres_obj.space.x.first, trans_spheres_obj.space.y.first, trans_spheres_obj.space.z.first};
    std::vector<double> uub = {trans_spheres_obj.space.x.second, trans_spheres_obj.space.y.second, trans_spheres_obj.space.z.second};

    for(size_t i = 0; i != llb.size(); i++)
    {
        std::cout << llb[i] << '\t' << uub[i] << std::endl;
    }

    two_spheres_sample = std::make_shared<trie_based::TrieBased<trie_based::NodeCount<size_t>, size_t>>();
    two_spheres_quant = std::make_shared<empirical_quantile::ImplicitQuantile<size_t, double>>(llb, uub, grid_n);

    std::vector<std::vector<double>> gridss(grid_n.size());
    std::vector<double> dxs(grid_n.size());

    for(size_t i = 0; i != gridss.size(); i++)
    {
        double grid_N = grid_n[i];
        std::vector<double> grid(grid_N);
        double startp = llb[i];
        double endp = uub[i];
        double es = endp - startp;
        for(size_t j = 0; j != grid.size(); j++)
        {
            grid[j] = startp + j*es/grid_N;
        }
        gridss[i] = grid;
        dxs[i] = es/(grid_N*2);
    }
//    std::vector<std::vector<double>> sampled;
    for(size_t i = 0; i != grid_n[0]; i++)
    {
        for(size_t j = 0; j != grid_n[1]; j++)
        {
            for(size_t k = 0; k != grid_n[2]; k++)
            {
                //std::cout << gridss[0][i] + dxs[0] << '\t' << gridss[1][j] + dxs[1] << '\t' << gridss[2][k] + dxs[2] << std::endl;
                bool in = false;
                for(size_t m = 0; m != trans_spheres_obj.spheres.size(); m++) // if 1 trans only to one first sphere
                {
                    double myDistance = std::sqrt(
                                            std::pow(gridss[0][i] + dxs[0] - trans_spheres_obj.spheres[m].x, 2.0) +
                                            std::pow(gridss[1][j] + dxs[1] - trans_spheres_obj.spheres[m].y, 2.0) +
                                            std::pow(gridss[2][k] + dxs[2] - trans_spheres_obj.spheres[m].z, 2.0));
                    if(myDistance < trans_spheres_obj.spheres[m].r)
                    {
                        in = true;
                        break;
                    }
                }
                if(in)
                {
                    std::vector<size_t> t = {i, j, k};
                    if(!two_spheres_sample->search(t))
                        two_spheres_sample->insert(t);
                }
            }
        }
    }

    two_spheres_quant->set_sample_shared(two_spheres_sample);

//    std::mt19937_64 generator;
//    generator.seed(1);
//    std::uniform_real_distribution<double> ureal01(0.0,1.0);
//
//    std::vector<std::vector<double>> values01;
//    std::vector<std::vector<double>> sampled;
//    std::vector<double> temp1(grid_n.size());
//    std::vector<double> temp2(temp1.size());
//    for(size_t i = 0; i != 10000; ++i)
//    {
//        for(size_t j = 0; j != temp1.size(); j++)
//        {
//            temp1[j] = ureal01(generator);
////            if(temp1[j] > 0.5)
////                temp1[j] = temp1[j]/2.0;
//        }
//        values01.push_back(temp1);
//        sampled.push_back(temp2);
//    }
//    for(size_t j = 0; j != values01.size(); j++)
//        two_spheres_quant->transform(values01[j], sampled[j]);
//    std::cout << two_spheres_sample->get_total_count() << std::endl;
//    write_default2d("spheres.data", sampled, 5);


    start2.resize(opt_vector.size() + 7);
    std::vector<double> sphere_center1 = {trans_spheres_obj.spheres[0].x, trans_spheres_obj.spheres[0].y, trans_spheres_obj.spheres[0].z};
    std::vector<double> sphere_center2 = {trans_spheres_obj.spheres[1].x, trans_spheres_obj.spheres[1].y, trans_spheres_obj.spheres[1].z};

    double rez = 0;
    for(size_t i = 0; i != start2.size(); i++)
        rez += std::pow(sphere_center1[i] - sphere_center2[i], 2.0);
    double width = 20.0/std::sqrt(rez);//20.0/rez;

    std::mt19937_64 generator_;
    generator_.seed(1);
    std::vector<double> start(lb.size());
    for(size_t j = 0; j != start.size(); j++)
    {
        std::cout << lb[j] << '\t' << ub[j] << std::endl;
        std::uniform_real_distribution<double> distribution(lb[j], ub[j]);
        start[j] = distribution(generator_);
    }

    cso::CSO obj1(1000, lb.size(), 0.01, 1e6, generator_);
    obj1.FitnessFunction = std::bind(&find_sphere_quant,
                                     std::placeholders::_1,
                                     two_spheres_quant,
                                     std::ref(pose.front()),
                                     std::ref(pose_shift.front()),
                                     std::ref(opt_vector),
                                     std::ref(param_list),
                                     width,
                                     std::ref(peptide_ranges),
                                     std::ref(trans_spheres_obj),
                                     1);
    do
    {
        obj1.set_bounds(lb, ub);
        obj1.init();
        obj1.optimization();
        std::cout << obj1.get_best() << std::endl;
        start = obj1.get_best_vector();
        //std::cin.get();
    }
    while(obj1.get_best() > -1.05);
    start1 = start;

    cso::CSO obj2(1000, lb.size(), 0.01, 1e6, generator_);
    obj2.FitnessFunction = std::bind(&find_sphere_quant,
                                     std::placeholders::_1,
                                     two_spheres_quant,
                                     std::ref(pose.front()),
                                     std::ref(pose_shift.front()),
                                     std::ref(opt_vector),
                                     std::ref(param_list),
                                     width,
                                     std::ref(peptide_ranges),
                                     std::ref(trans_spheres_obj),
                                     0);
    do
    {
        obj2.set_bounds(lb, ub);
        obj2.init();
        obj2.optimization();
        std::cout << obj2.get_best() << std::endl;
        start = obj2.get_best_vector();
        //std::cin.get();
    }
    while(obj2.get_best() > -1.05);
    start2 = start;
}

}
