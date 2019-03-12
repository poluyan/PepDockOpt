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
#ifndef INCLUDED_bbtools_hh
#define INCLUDED_bbtools_hh

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>

namespace pepdockopt
{
namespace bbtools
{

core::Real to_positive_radians(core::Real radian);
core::Real to_negative_radians(core::Real radian);
core::Real normalize_to_mpi_to_ppi(core::Real radian);
void make_ideal_peptide(core::pose::Pose& ideal_peptide, const core::pose::Pose& peptide);
void to_centroid(core::pose::Pose &pose);

}
}

#endif
