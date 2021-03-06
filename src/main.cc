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
#include <devel/init.hh>

#include <iostream>
#include <vector>
#include <random>
#include <limits>

#include <pepdockopt.hh>

int main(int argc, char *argv[])
{    
    devel::init(argc, argv);
    std::cout << "Start..." << std::endl;

    size_t thread_num = 4;

    pepdockopt::PepDockOpt obj;
    obj.init(thread_num);
    obj.set_opt();
    obj.set_grid();
    obj.sphere_quant();
    obj.set_quantile1();
    obj.set_quantile2();
    obj.check();
    
    obj.set_objective();
    obj.set_bbdep(72);
    obj.set_omega_quantile(360);
        
    std::vector<double> x(obj.get_objective_dimension(),0.6);
    std::cout << obj.objective(x,0) << std::endl;
}
