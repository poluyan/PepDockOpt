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

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/rigid/RigidBodyMover.hh>

#include <complex.hh>
#include <spheres.hh>

#include <quantile.hh>

#include <iostream>

namespace pepdockopt
{
    
struct PoseShift
{
    protocols::rigid::RigidBodyDeterministicSpinMover SpinMover;
    core::kinematics::Jump FlexibleJump;
    core::Vector InitPeptidePosition;
    core::kinematics::RT::Matrix InitRm;
    core::kinematics::Stub UpstreamStub;
};

class PepDockOpt
{
protected:
    int threads_number;
    
    std::vector<core::pose::Pose> pose;
    core::pose::Pose native;
    
    pepdockopt::ComplexInfoNseq param_list;
    pepdockopt::ranges peptide_ranges;
    pepdockopt::ranges protein_ranges;
    std::vector<PoseShift> pose_shift;
    
    std::vector<core::scoring::ScoreFunctionOP> score_func;
        
    std::vector<pepdockopt::opt_element> opt_vector;
    std::vector<double> lb;
    std::vector<double> ub;
    std::vector<double> start1;
    std::vector<double> start2;
    
    std::vector<std::shared_ptr<trie_based::TrieBased<trie_based::NodeCount<int>,int>>> phipsi_rama2_sample;
    std::vector<std::shared_ptr<empirical_quantile::ImplicitQuantile<int, double>>> phipsi_rama2_quantile;
    std::vector<std::shared_ptr<empirical_quantile::ImplicitQuantile<int, double>>> omega_quantile;
    
    pepdockopt::spheres::box_trans trans_spheres_obj1;
    pepdockopt::spheres::box_trans trans_spheres_obj2;
    
    std::shared_ptr<trie_based::TrieBased<trie_based::NodeCount<size_t>, size_t>> two_spheres_sample;
    std::shared_ptr<empirical_quantile::ImplicitQuantile<size_t, double>> two_spheres_quant;
    
    std::vector<size_t> gridN;
    std::vector<std::vector<double>> grids;
    std::vector<double> dx;
    
    typedef trie_based::TrieBased<trie_based::NodeCount<std::uint8_t>,std::uint8_t> frag_type;  
    std::shared_ptr<frag_type> structures_triebased;
    std::shared_ptr<empirical_quantile::ImplicitQuantile<std::uint8_t, double>> structures_quant;
    
    void set_score_function();
public:
    PepDockOpt();
    void init(size_t _threads_number);
    void set_opt();
    void start_position1();
    void start_position2();
    void set_grid();
    void set_quantile1();
    void set_quantile2();
    void check();
    
    void sphere_quant();
    
    std::vector<double> get_position(std::vector<double> _lb, std::vector<double> _ub, double width, spheres::box_trans trans_sp, std::pair<size_t, size_t> spheres_number); 
};

}


#endif
