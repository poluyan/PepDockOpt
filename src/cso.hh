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
#ifndef INCLUDED_cso_hh
#define INCLUDED_cso_hh

#include <random>
#include <iostream>
#include <vector>

namespace pepdockopt
{
namespace cso
{
class CSO
{
public:
    void set_generator_seed(int _seed)
    {
        generator.seed(_seed);
    }
    void init()
    {
        for(size_t i = 0; i < p.size(); i++)
        {
            for(size_t j = 0; j < p[0].size(); j++)
            {
                std::uniform_real_distribution<double> uniform_real_dist(lb[j], ub[j]);
                p[i][j] = uniform_real_dist(generator);
            }
        }

        for(size_t i = 0; i < fitness.size(); i++)
        {
            fitness[i] = FitnessFunction(p[i]);
            if(best > fitness[i] && fe_count < fe_max)
            {
                best = fitness[i];
                bestvector = p[i];
                std::cout << fe_count << '\t' << best << std::endl;
            }
            fe_count++;
        }
    }
    void optimization()
    {
        if(best < -1.05)
            return;
        for(size_t g = 0; fe_count < fe_max; g++)
        {
            std::vector<size_t> rlist(pop_size);
            for(size_t i = 0; i < rlist.size(); i++)
            {
                rlist[i] = i;
            }
            std::shuffle(rlist.begin(), rlist.end(), generator);

            std::vector<std::pair<size_t, size_t> > rpairs(pop_size/2);
            for(size_t i = 0, j = pop_size/2; i < rpairs.size(); i++, j++)
            {
                rpairs[i].first = rlist[i];
                rpairs[i].second = rlist[j];
            }

            std::vector<double> sum(dim, 0.0);
            std::vector<double> mean_p(dim, 0.0);
            for(size_t i = 0; i < p.size(); i++)
            {
                for(size_t j = 0; j < p[0].size(); j++)
                {
                    sum[j] += p[i][j];
                }
            }
            for(size_t i = 0; i < mean_p.size(); i++)
            {
                mean_p[i] = sum[i]/p.size();
            }

            std::vector<std::vector<double> > center(p.size()/2, std::vector<double>(dim));
            for(size_t i = 0; i < center.size(); i++)
            {
                for(size_t j = 0; j < center[0].size(); j++)
                {
                    center[i][j] = mean_p[j];
                }
            }

            std::vector<int> mask(rpairs.size(), 0);
            std::vector<size_t> losers(rpairs.size());
            std::vector<size_t> winners(rpairs.size());

            for(size_t i = 0; i < mask.size(); i++)
            {
                if(fitness[rpairs[i].first] > fitness[rpairs[i].second])
                {
                    mask[i] = 1;
                    losers[i] = rpairs[i].first;
                    winners[i] = rpairs[i].second;
                }
                else
                {
                    mask[i] = 0;
                    losers[i] = rpairs[i].second;
                    winners[i] = rpairs[i].first;
                }
            }

            std::vector<std::vector<double> > randco1(p.size()/2, std::vector<double>(dim));
            std::vector<std::vector<double> > randco2(p.size()/2, std::vector<double>(dim));
            std::vector<std::vector<double> > randco3(p.size()/2, std::vector<double>(dim));

            for(size_t i = 0; i < randco1.size(); i++)
            {
                for(size_t j = 0; j < randco1[0].size(); j++)
                {
                    std::uniform_real_distribution<double> uniform_real_01(0, 1);
                    randco1[i][j] = uniform_real_01(generator);
                    randco2[i][j] = uniform_real_01(generator);
                    randco3[i][j] = uniform_real_01(generator);
                }
            }

            for(size_t i = 0; i < center.size(); i++)
            {
                for(size_t j = 0; j < center[0].size(); j++)
                {
                    center[i][j] = phi*randco3[i][j]*(center[i][j] - p[losers[i]][j]);
                }
            }

            for(size_t i = 0; i < randco1.size(); i++)
            {
                for(size_t j = 0; j < randco1[0].size(); j++)
                {
                    v[losers[i]][j] = randco1[i][j]*v[losers[i]][j] +
                                      randco2[i][j]*(p[winners[i]][j] - p[losers[i]][j]) +
                                      center[i][j];
                }
            }

            for(size_t i = 0; i < losers.size(); i++)
            {
                for(size_t j = 0; j < p[0].size(); j++)
                {
                    p[losers[i]][j] = p[losers[i]][j] + v[losers[i]][j];
                }
            }

            for(size_t i = 0; i < p.size(); i++)
            {
                for(size_t j = 0; j < p[0].size(); j++)
                {
                    if(p[i][j] > ub[j])
                        p[i][j] = ub[j];
                    if(p[i][j] < lb[j])
                        p[i][j] = lb[j];
                }
            }

            for(size_t i = 0; i < losers.size(); i++)
            {
                fitness[losers[i]] = FitnessFunction(p[losers[i]]);
            }

            for(size_t i = 0; i < losers.size(); i++)
            {
                if(best > fitness[losers[i]] && fe_count < fe_max)
                {
                    best = fitness[losers[i]];
                    std::cout << fe_count << '\t' << best << std::endl;
                    bestvector = p[losers[i]];
                }
                fe_count++;
            }

            if(best < -1.05)
            {
                break;
            }
        }
    }
    CSO(size_t _pop_size, size_t _dim, double _phi, size_t _fe_max, std::mt19937_64 &generator_)
    {
        pop_size = _pop_size;
        dim = _dim;
        phi = _phi;

        fe_max = _fe_max;
        fe_count = 0;
        best = std::numeric_limits<double>::max();

        p.resize(pop_size);
        for(size_t i = 0; i < p.size(); i++)
        {
            p[i].resize(dim);
        }

        v.resize(pop_size);
        for(size_t i = 0; i < v.size(); i++)
        {
            v[i].resize(dim);
        }

        fitness.resize(pop_size);

        generator = generator_;
    }
    void set_bounds(std::vector<double> _lb, std::vector<double> _ub)
    {
        lb = _lb;
        ub = _ub;
    }
    void set_population(std::vector<std::vector<double> >& _p)
    {
        if(_p.size() != pop_size || _p[0].size() != dim)
        {
            std::cout << "error!\n";
        }
        p = _p;
    }
    std::vector<double> get_best_vector()
    {
        return bestvector;
    }
    std::vector<std::vector<double> > get_population()
    {
        return p;
    }
    std::vector<double> get_fitness()
    {
        return fitness;
    }
    void set_fe_max(size_t _fe_max)
    {
        fe_max = _fe_max;
    }
    void reset_fe_count()
    {
        fe_count = 0;
    }
    double get_best()
    {
        return best;
    }
    std::function<double(std::vector<double>)> FitnessFunction;
protected:
    size_t pop_size;
    size_t dim;

    double phi;

    size_t fe_max;
    size_t fe_count;
    double best;

    std::vector<std::vector<double> > p;
    std::vector<std::vector<double> > v;
    std::vector<double> fitness;

    std::vector<double> lb;
    std::vector<double> ub;

    std::vector<double> bestvector;

    std::mt19937_64 generator;
};

}
}


#endif