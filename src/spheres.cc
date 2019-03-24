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
#include <spheres.hh>

#include <vector>
#include <deque>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>
#include <functional>   // std::multiplies
#include <numeric>      // std::partial_sum

namespace pepdockopt
{
namespace spheres
{

box_trans::box_trans() {}

void box_trans::load_data(std::string fname, size_t steps, bool change_spheres)
{
    std::ifstream fIn;
    fIn.open(fname.c_str());
    if(!fIn.is_open())
    {
        std::cout << "Error opening file." << std::endl;
    }
    std::string temp;

    while(!fIn.eof())
    {
        fIn >> temp;
        if(temp.find("Sphere") != std::string::npos)
        {
            sphere_info sphere;
            fIn >> sphere.x;
            fIn >> sphere.y;
            fIn >> sphere.z;
            fIn >> sphere.r;
            spheres.push_back(sphere);
        }
    }
    fIn.close();
    
    if(change_spheres)
    {
        std::reverse(spheres.begin(), spheres.end());
    }

    auto minmax_x = std::minmax_element(spheres.begin(), spheres.end(),
        [](sphere_info const &lhs, sphere_info const &rhs) { return lhs.x < rhs.x; });
    auto minmax_y = std::minmax_element(spheres.begin(), spheres.end(),
        [](sphere_info const &lhs, sphere_info const &rhs) { return lhs.y < rhs.y; });
    auto minmax_z = std::minmax_element(spheres.begin(), spheres.end(),
        [](sphere_info const &lhs, sphere_info const &rhs) { return lhs.z < rhs.z; });
    auto minmax_r = std::minmax_element(spheres.begin(), spheres.end(),
        [](sphere_info const &lhs, sphere_info const &rhs) { return lhs.r < rhs.r; });

    space.x = std::make_pair(minmax_x.first->x - minmax_r.second->r, minmax_x.second->x + minmax_r.second->r);
    space.y = std::make_pair(minmax_y.first->y - minmax_r.second->r, minmax_y.second->y + minmax_r.second->r);
    space.z = std::make_pair(minmax_z.first->z - minmax_r.second->r, minmax_z.second->z + minmax_r.second->r);
    
    max_r = minmax_r.second->r;
    
    dist = make_cdf(steps);
}

bbutils::distribution_3d box_trans::make_cdf(size_t steps)
{
    std::vector<double> es;
    std::vector<size_t> step(3, steps);

    es.push_back(space.x.second - space.x.first);
    es.push_back(space.y.second - space.y.first);
    es.push_back(space.z.second - space.z.first);

    std::deque<std::deque<std::deque<double>>> pdf;
    std::deque<double> x, y, z;
    for(size_t i = 0; i != step[0]; i++)
        x.push_back(space.x.first + i * es[0] / step[0]);
    for(size_t i = 0; i != step[1]; i++)
        y.push_back(space.y.first + i * es[1] / step[1]);
    for(size_t i = 0; i != step[2]; i++)
        z.push_back(space.z.first + i * es[2] / step[2]);

    x.push_back(space.x.second);
    y.push_back(space.y.second);
    z.push_back(space.z.second);

    for(size_t i = 0; i != step[0]; i++)
    {
        std::deque<std::deque<double>> temp1;
        for(size_t j = 0; j != step[1]; j++)
        {
            std::deque<double> temp2;
            for(size_t k = 0; k != step[2]; k++)
            {
                bool in = false;
                double value = 0.0;
                for(size_t m = 0; m != /*spheres.size()*/1; m++) // if 1 trans only to one first sphere
                {
                    double myDistance = std::sqrt(
                        std::pow(space.x.first + i * (space.x.second - space.x.first) / step[0] - spheres[m].x, 2.0) +
                        std::pow(space.y.first + j * (space.y.second - space.y.first) / step[1] - spheres[m].y, 2.0) +
                        std::pow(space.z.first + k * (space.z.second - space.z.first) / step[2] - spheres[m].z, 2.0));

                    if(myDistance < spheres[m].r)
                    {
                        //value = spheres[m].r - myDistance;
                        value = 1.0;
                        in = true;
                        break;
                    }
                }
                if(in)
                    temp2.push_back(value);
                else
                    temp2.push_back(0);
            }
            temp1.push_back(temp2);
        }
        pdf.push_back(temp1);
    }
    //print2file_3d_d("/work/ProjectsCPP/two_spheres/maps/spheres.dat", pdf);

    // col cdf

    std::deque<std::deque<double>> first_2d;
    for(size_t i = 0; i != pdf.size(); i++)
    {
        std::deque<double> temp1;
        for(size_t j = 0; j != pdf[i].size(); j++)
        {
            double t = 0.0;
            for(size_t k = 0; k != pdf[i][j].size(); k++)
            {
                t += pdf[i][j][k];
            }
            temp1.push_back(t);
        }
        first_2d.push_back(temp1);
    }

    std::deque<std::deque<double>> first_2d_temp = first_2d;
    for(size_t i = 0; i != pdf.size(); i++)
    {
        std::partial_sum(
            first_2d_temp[i].begin(), first_2d_temp[i].end(), first_2d_temp[i].begin(), std::plus<double>());
    }
    std::deque<double> col_pdf;
    for(size_t i = 0; i != pdf.size(); i++)
    {
        col_pdf.push_back(first_2d_temp[i].back());
    }

    double S = std::accumulate(col_pdf.begin(), col_pdf.end(), 0.0);
    for(auto &a : col_pdf)
    {
        a /= S;
    }

    std::deque<double> col_cdf(col_pdf.size());
    std::partial_sum(col_pdf.begin(), col_pdf.end(), col_cdf.begin(), std::plus<double>());
    col_cdf.push_front(0);

    std::deque<std::deque<double>> row_dist_bank1(pdf.size());
    for(size_t i = 0; i != row_dist_bank1.size(); i++)
    {
        row_dist_bank1[i] = first_2d[i];

        double S = std::accumulate(row_dist_bank1[i].begin(), row_dist_bank1[i].end(), 0.0);
        for(auto &a : row_dist_bank1[i])
        {
            a /= S;
        }

        std::partial_sum(
            row_dist_bank1[i].begin(), row_dist_bank1[i].end(), row_dist_bank1[i].begin(), std::plus<double>());
        row_dist_bank1[i].push_front(0);
    }

    std::deque<std::deque<std::deque<double>>> row_dist_bank2(
        pdf.size(), std::deque<std::deque<double>>(pdf[0].size()));

    for(size_t i = 0; i != row_dist_bank2.size(); i++)
    {
        for(size_t j = 0; j != row_dist_bank2[i].size(); j++)
        {
            row_dist_bank2[i][j] = pdf[i][j];

            double S = std::accumulate(row_dist_bank2[i][j].begin(), row_dist_bank2[i][j].end(), 0.0);
            for(auto &a : row_dist_bank2[i][j])
            {
                a /= S;
            }

            std::partial_sum(row_dist_bank2[i][j].begin(), row_dist_bank2[i][j].end(), row_dist_bank2[i][j].begin(),
                std::plus<double>());
            row_dist_bank2[i][j].push_front(0);
        }
    }

    bbutils::distribution_3d result;
    // result.pdf = pdf;

    result.row_dist1 = row_dist_bank1;
    result.col_dist = col_cdf;

    result.row_dist2 = row_dist_bank2;

    result.grid.push_back(x);
    result.grid.push_back(y);
    result.grid.push_back(z);

    return result;
}
}
}