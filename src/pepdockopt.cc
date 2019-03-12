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

#include <vector>
#include <random>
#include <omp.h>

namespace pepdockopt
{

PepDockOpt::PepDockOpt() {}


double Quadratic_ValueSimple2(std::vector<double> x,
                            core::pose::Pose& pose,
        protocols::rigid::RigidBodyDeterministicSpinMover& SpinMover,
        std::vector<opt_element> &opt_vector,
        core::kinematics::Jump &flexible_jump,
        core::kinematics::Stub &upstream_stub,
        core::Vector &init_peptide_position,
        core::kinematics::RT::Matrix &init_rm,
        ComplexInfoNseq &complex_stuff,
        spheres::box_trans &trans_obj,
        double width,
        pepdockopt::ranges peptide_ranges,
        spheres::box_trans trans_spheres_obj,
        std::pair<size_t, size_t> sphere_number)
{
    bool do_trans = true;
    bool must_do_experiment = true;
    if(must_do_experiment)
    {

//        x = transformation::bbdep_experiment_peptide_actual_states(x, OPT_VECTOR_ALL_INFO, peptide_ranges, bbdep_obj_em, peptide_phipsi_2d,
//                param_list.get_peptide_first_index(), param_list.get_peptide_last_index());

        double pi = numeric::NumericTraits<core::Real>::pi();
        x[std::get<1>(peptide_ranges.phipsi)] = pi*(2.0*x[std::get<1>(peptide_ranges.phipsi)] - 1.0);
        //x[7] = pi*(2.0*x[7] - 1.0);
        x[std::get<2>(peptide_ranges.phipsi) - 1] = pi*(2.0*x[std::get<2>(peptide_ranges.phipsi) - 1] - 1.0);

        //            pokt = transformation::bbdep_experiment_protein_actual_states(pokt, OPT_VECTOR_ALL_INFO, protein_ranges, bbdep_obj_em, cm_fixed_phipsi, param_list.get_protein_first_indices(),
        //                    param_list.get_protein_last_indices());

        //std::cout << "ok2" << std::endl;

        x = transform::peptide_quaternion(x, opt_vector, opt_vector.size());

        //std::cout << "ok3" << std::endl;

        if(do_trans)
            x = transform::twospheres(x, opt_vector.size(), trans_spheres_obj);

        //std::cout << "ok4" << std::endl;
    }

    for(size_t i = 0, end = opt_vector.size(); i < end; i++)
    {
        pose.set_dof(opt_vector[i].dofid, x[i]);
        //std::cout << i << '\t' << x[i] << std::endl;
    }

    ///

    core::Vector new_position = init_peptide_position;
    new_position.at(0) = x[opt_vector.size()];
    new_position.at(1) = x[opt_vector.size() + 1];
    new_position.at(2) = x[opt_vector.size() + 2];

    flexible_jump.reset();

    flexible_jump.set_rotation(init_rm);
    pose.set_jump(pose.num_jump(), flexible_jump);
    
    numeric::Real current_distance(pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA").distance(new_position));
    flexible_jump.translation_along_axis(upstream_stub, new_position - pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA"), current_distance);
    pose.set_jump(pose.num_jump(), flexible_jump);

    core::Vector peptide_c_alpha_centroid(pose.residue(complex_stuff.get_peptide_first_index()).xyz("CA"));
    core::Vector axis(x[opt_vector.size() + 3], x[opt_vector.size() + 4], x[opt_vector.size() + 5]);

    SpinMover.rot_center(peptide_c_alpha_centroid);
    SpinMover.spin_axis(axis);
    SpinMover.angle_magnitude(x.back());
    SpinMover.apply(pose);

    ///

    std::vector<double> exponential_center = {trans_obj.spheres[sphere_number.second].x,
                                              trans_obj.spheres[sphere_number.second].y,
                                              trans_obj.spheres[sphere_number.second].z
                                             };

    core::Vector lCA = pose.residue(complex_stuff.get_peptide_last_index()).xyz("CA");

    double rez = 0;
    double temp = 0;
    for(size_t i = 0; i != exponential_center.size(); i++)
    {
        rez += std::pow(lCA.at(i) - exponential_center[i], 2.0);
    }
    temp = std::sqrt(rez);
    if(temp < trans_obj.spheres[sphere_number.second].r)
    {
        rez = 1.1;

        //std::cout << "center: " << exponential_center[0] << '\t' << exponential_center[1] << '\t' << exponential_center[2] << std::endl;
        //std::cout << "last r: " << lCA.at(0) << '\t' << lCA.at(1) << '\t' << lCA.at(2) << std::endl;
        //std::cout << std::endl;
    }
    else
        rez = std::exp(-temp*width);

    //std::cout << bfgs_fe_count << '\t' << -rez << std::endl;

    //std::cout << "stop" << std::endl;
    //std::cin.get();
    return -rez;
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
    void run()
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
        double F = FitnessFunction(x);

        double theta = 1.0;
        double dtheta = 0.999;
        std::vector<double> dx(lb.size());
        for(size_t i = 0; i != dx.size(); i++)
        {
            dx[i] = (ub[i] - lb[i])/100.0;
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
            double FF = FitnessFunction(xx);
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
                dx[i] -= ((ub[i] - lb[i])/100.0) / fe_max;
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
    double (*FitnessFunction)(const std::vector<double> &x) = nullptr;
private:
    size_t dim;

    size_t fe_max;
    double best;
    std::vector<double> best_vector;
    std::vector<double> lb;
    std::vector<double> ub;

    std::mt19937_64 generator;
};


std::vector<double> PepDockOpt::get_position()
{
    std::vector<double> start_temp(start.size());

    double best = 0;
    bool is_in_sphere = false;
    bool reduce_step = true;

    for(size_t j = 0; j != start.size(); j++)
    {
        std::cout << lb[j] << '\t' << ub[j] << std::endl;
        std::uniform_real_distribution<double> distribution(lb[j], ub[j]);
        start[j] = distribution(generator_);
    }

    _GlobalObjs = GlobalObj;
    _SpinMovers = SpinMover;
    _opt_vectors = OPT_VECTORS;
    _flexible_jumps = FlexibleJumps;
    _upstream_stubs = UpstreamStubs;
    _init_peptide_positions = InitPeptidePositions;
    _init_rms = InitRms;
    _complex_stuff = param_list;
    _trans_obj = trans_spheres_obj;
    _width = width;
    _sphere_numbers = spheres_number;
    _is_in_sphere = is_in_sphere;

    CSO::CSO obj1(1000, start.size(), 0.01, 1e6, generator_);
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
    return
}
*/
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
    SpinMover.resize(threads_number);
    for(int j = 0; j != threads_number; j++)
    {
        SpinMover[j].rb_jump(pose[j].num_jump());
    }
    set_score_function();

    /// load input structures
    core::pose::Pose input;
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

    for(auto &i : pose)
        i = param_list.get_pose_complex(); //param_list.get_peptide();

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


    //
    std::vector<double> lb, ub;
    numeric::Real omega_delta = 0.2;
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

    std::vector<std::vector<core::Real>> rb_peptide_space;
    switch(use_box)
    {
        case 0:
        {
            /// r >= 0
            lb.push_back(0);
            //ub.push_back(param_list.get_max_distance_from_center_of_protein() + param_list.get_extended_peptide_length() + 7.0);
            ub.push_back(0.1);

            std::cout << "Sphere radius = " << ub.back() << std::endl;

            /// tau [0, pi]
            lb.push_back(0);
            ub.push_back(numeric::NumericTraits<core::Real>::pi());

            /// phi [-pi, pi]
            lb.push_back(-numeric::NumericTraits<core::Real>::pi());
            ub.push_back(numeric::NumericTraits<core::Real>::pi());

            k += 3;
            break;
        }
        case 1:
        {
            numeric::Real box_size(60.0); //60
            core::kinematics::Jump flexible_jump = pose.front().jump(pose.front().num_jump());
            core::Vector current_peptide_position = flexible_jump.get_translation();

            for(size_t i = 0; i != 3; i++, k++)
            {
                lb.push_back(current_peptide_position.at(i) - box_size);
                ub.push_back(current_peptide_position.at(i) + box_size);

                std::vector<numeric::Real> temp(2);
                temp[0] = lb.back();
                temp[1] = ub.back();
                rb_peptide_space.push_back(temp);
            }
            break;
        }
        case 2:
        {
            for(size_t i = 0; i != 3; i++, k++)
            {
                lb.push_back(0.0);
                ub.push_back(1.0);
            }
            break;
        }
    }

    if(peptide_ranges.do_quaternion)
    {
        for(size_t i = 0; i != 3; i++, k++)
        {
            lb.push_back(0.0);
            ub.push_back(1.0);
        }

        lb.push_back(0.0);
        ub.push_back(1.0);
        k++;
    }
    else
    {
        for(size_t i = 0; i != 3; i++, k++)
        {
            lb.push_back(-1.0);
            ub.push_back(1.0);
        }

        lb.push_back(0.0);
        ub.push_back(180.0);
        k++;
    }

    if(k != opt_vector.size() + 7 || lb.size() - 7 != opt_vector.size() || ub.size() - 7 != opt_vector.size())
    {
        std::cout << "fail\n";
        return;
    }

    std::cout << opt_vector.size() << std::endl;

    int dim = opt_vector.size() + 7;
    std::cout << "dim = " << dim << std::endl;

    bool ca_trans = true;
    for(int j = 0; j != threads_number; j++)
    {
        FlexibleJumps.push_back(pose[j].jump(pose[j].num_jump()));

        if(ca_trans)
        {
            InitPeptidePositions.push_back(pose[j].residue(param_list.get_peptide_first_index()).xyz("CA"));
        }
        else
        {
            InitPeptidePositions.push_back(FlexibleJumps.back().get_translation());
        }

        InitRms.push_back(FlexibleJumps.back().get_rotation());

        UpstreamStubs.push_back(pose[j].conformation().upstream_jump_stub(pose[j].num_jump()));
    }


    std::pair<size_t, size_t> spheres_number = std::make_pair(0, 1);

    std::mt19937_64 generator_;
    generator_.seed(1);

    std::vector<double> start(lb.size());

    std::vector<double> sphere_center1 = {trans_spheres_obj.spheres[0].x, trans_spheres_obj.spheres[0].y, trans_spheres_obj.spheres[0].z};
    std::vector<double> sphere_center2 = {trans_spheres_obj.spheres[1].x, trans_spheres_obj.spheres[1].y, trans_spheres_obj.spheres[1].z};

    double rez = 0;
    for(size_t i = 0; i != start.size(); i++)
        rez += std::pow(sphere_center1[i] - sphere_center2[i], 2.0);
    double width = 100.0/rez;//20.0/rez;

    size_t N = 1e7;

    std::vector<double> step_sizes(lb.size());
    for(size_t i = 0; i != step_sizes.size() - 1; i++)
    {
        step_sizes[i] = std::abs(ub[i] - lb[i])/100.0;
    }
    step_sizes[step_sizes.size() - 1] = std::abs(ub[step_sizes.size() - 1] - lb[step_sizes.size() - 1])/1000.0;

    step_sizes[step_sizes.size() - 5] = std::abs(ub[step_sizes.size() - 5] - lb[step_sizes.size() - 5])/1000.0;
    step_sizes[step_sizes.size() - 6] = std::abs(ub[step_sizes.size() - 6] - lb[step_sizes.size() - 6])/1000.0;
    step_sizes[step_sizes.size() - 7] = std::abs(ub[step_sizes.size() - 7] - lb[step_sizes.size() - 7])/1000.0;

    //mh(generator_, N, start, shift, width, true, sphere_radius, 0.1);

    std::normal_distribution<double> norm01(0, 1);
    std::uniform_real_distribution<double> ureal01(0, 1);

    peptide_ranges.do_chi = false;

    //std::vector<std::vector<double>> samples(N, std::vector<double>(start.size()));
    std::vector<std::vector<double>> best_samples;
    //samples[0] = start;

//
//    get_position();
//
//    std::vector<double> gridN(start.size());
//    gridN[0] = 8;//psi
//    gridN[1] = 8;//2
//    gridN[2] = 8;//2
//    gridN[3] = 8;//3
//    gridN[4] = 8;//3
//    gridN[5] = 8;//4
//    gridN[6] = 8;//4
//    gridN[7] = 8;//phi
//
//    gridN[8] = 8;
//    gridN[9] = 8;
//    gridN[10] = 8;
//    gridN[11] = 8;
//    gridN[12] = 8;
//    gridN[13] = 8;
//    //gridN[14] = 8;
//
//    /*gridN[0] = 50;//psi
//    gridN[1] = 50;//2
//    gridN[2] = 50;//2
//    gridN[3] = 50;//phi
//
//    gridN[4] = 1;
//    gridN[5] = 1;
//    gridN[6] = 1;
//    gridN[7] = 1;
//    gridN[8] = 1;
//    gridN[9] = 1;
//    gridN[10] = 1;*/
//
//    std::vector<std::vector<double>> grids(gridN.size());
//    std::vector<double> dx(gridN.size());
//
//    for(size_t i = 0; i != grids.size(); i++)
//    {
//        size_t grid_N = gridN[i];
//        std::vector<double> grid(grid_N);
//        double startp = 0;
//        double endp = 1;
//        double es = endp - startp;
//        for(size_t j = 0; j != grid.size(); j++)
//        {
//            grid[j] = startp + j*es/(grid_N);
//        }
//        grids[i] = grid;
//        dx[i] = es/(grid_N*2);
//    }
//
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
//
//    std::vector<char> startdot(start.size());
//    for(size_t i = 0; i != startdot.size(); i++)
//    {
//        auto pos1 = std::lower_bound(grids[i].begin(), grids[i].end(), start[i]);
//        startdot[i] = std::distance(grids[i].begin(), pos1) - 1;
//    }

}


}
