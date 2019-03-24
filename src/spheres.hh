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
#ifndef INCLUDED_spheres_hh
#define INCLUDED_spheres_hh

#include <bbutils.hh>
#include <string>

namespace pepdockopt
{
namespace spheres
{

struct sphere_info
{
    double x;
    double y;
    double z;
    double r;
};

struct box
{
    std::pair<double, double> x;
    std::pair<double, double> y;
    std::pair<double, double> z;
};

class box_trans
{
public:
    box space;

    std::vector<sphere_info> spheres;
    bbutils::distribution_3d dist;
    double max_r;

    box_trans();
    void load_data(std::string fname, size_t steps, bool change_spheres);
    bbutils::distribution_3d make_cdf(size_t steps);
};

}
}

#endif
