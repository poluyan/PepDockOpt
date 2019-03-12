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
#include <bbutils.hh>

namespace pepdockopt
{
namespace bbutils
{

double get_1d_from_dst(const distribution_1d &dst, double value)
{
    size_t index = std::distance(dst.cdf.begin(), std::upper_bound(dst.cdf.begin(), dst.cdf.end(), value)) - 1;
    double x0 = dst.grid[index], y0 = dst.cdf[index], x1 = dst.grid[index + 1], y1 = dst.cdf[index + 1];
    return x0 + (value - y0) * (x1 - x0) / (y1 - y0);
}
double get_inverse_1d_from_dst(const distribution_1d &dst, double value) // here value is the degree, must be from -180 to 180
{
    size_t index = std::distance(dst.grid.begin(), std::upper_bound(dst.grid.begin(), dst.grid.end(), value)) - 1;
    double x0 = dst.cdf[index], y0 = dst.grid[index], x1 = dst.cdf[index + 1], y1 = dst.grid[index + 1];
    return x0 + (value - y0) * (x1 - x0) / (y1 - y0);
}
std::vector<double> get_3d_from_dst(const distribution_3d &dst, const std::vector<double> &values)
{
    size_t index1 = std::distance(dst.col_dist.begin(), std::upper_bound(dst.col_dist.begin(), dst.col_dist.end(), values[0])) - 1;

    double x0 = dst.grid[0][index1], y0 = dst.col_dist[index1], x1 = dst.grid[0][index1 + 1], y1 = dst.col_dist[index1 + 1];
    double fx = x0 + (values[0] - y0) * (x1 - x0) / (y1 - y0);

    size_t index2 = std::distance(dst.row_dist1[index1].begin(), std::upper_bound(dst.row_dist1[index1].begin(), dst.row_dist1[index1].end(), values[1])) - 1;

    x0 = dst.grid[1][index2];
    y0 = dst.row_dist1[index1][index2];
    x1 = dst.grid[1][index2 + 1];
    y1 = dst.row_dist1[index1][index2 + 1];
    double fy = x0 + (values[1] - y0) * (x1 - x0) / (y1 - y0);

    size_t index3 = std::distance(dst.row_dist2[index1][index2].begin(), std::upper_bound(dst.row_dist2[index1][index2].begin(), dst.row_dist2[index1][index2].end(), values[2]))
                    - 1;

    x0 = dst.grid[2][index3];
    y0 = dst.row_dist2[index1][index2][index3];
    x1 = dst.grid[2][index3 + 1];
    y1 = dst.row_dist2[index1][index2][index3 + 1];
    double fz = x0 + (values[2] - y0) * (x1 - x0) / (y1 - y0);

    std::vector<double> result = { fx, fy, fz };
    return result;
}

}
}
