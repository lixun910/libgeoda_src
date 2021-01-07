//
// Created by Xun Li on 9/26/19.
//

#include <boost/algorithm/string.hpp>
#include "../weights/GalWeight.h"
#include "../GenUtils.h"
#include "cluster.h"
#include "azp.h"
#include "maxp_wrapper.h"

maxp_wrapper::~maxp_wrapper() {

}

maxp_wrapper::maxp_wrapper(GeoDaWeight *w,
                           const std::vector<std::vector<double> >& data,
                           int iterations,
                           int inits,
                           const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                           const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                           const std::vector<int>& _init_regions,
                           const std::string &distance_method,
                           int _rnd_seed)
{
    num_obs = 0;
    init_regions = _init_regions;
    rnd_seed = _rnd_seed;

    if (w) {
        num_obs = w->num_obs;
        gal = Gda::GetGalElement(w);
        if (gal) {
            // get distance matrix
            n_cols = data.size();
            input_data = new double*[num_obs];
            int** mask = new int*[num_obs];
            for (size_t i=0; i<num_obs; ++i) {
                input_data[i] = new double[n_cols];
                mask[i] = new int[n_cols];
                for (size_t j=0; j<n_cols; ++j) mask[i][j] = 1.0;
            }
            for (size_t i=0; i<n_cols; ++i) {
                std::vector<double> vals = data[i];
                GenUtils::StandardizeData(vals);
                for (size_t r=0; r<num_obs; ++r) {
                    input_data[r][i] = vals[r];
                }
            }
            char dist = 'e';
            if (boost::iequals(distance_method, "manhattan")) dist = 'b';
            int transpose = 0; // row wise
            double* weight = new double[n_cols];
            for (size_t i=0; i<n_cols; ++i) weight[i] = 1.0;

            double** ragged_distances = distancematrix(num_obs, n_cols, input_data,  mask, weight, dist, transpose);
            dm = new RawDistMatrix(ragged_distances);

            // create bounds
            CreateController(min_bounds, max_bounds);

            MaxpRegion* maxp = RunMaxp();

            std::vector<int> final_solution = maxp->GetResults();

            delete maxp;

            std::map<int, std::vector<int> > solution;
            for (int i=0; i<final_solution.size(); ++i) {
                solution[final_solution[i]].push_back(i);
            }
            std::map<int, std::vector<int> >::iterator it;
            for (it = solution.begin(); it != solution.end(); ++it) {
                cluster_ids.push_back(it->second);
            }

            for (int i = 1; i < num_obs; i++) free(ragged_distances[i]);
            free(ragged_distances);

            delete dm;
        }
    }
}

void maxp_wrapper::CreateController(const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                    const std::vector<std::pair<double, std::vector<double> > >& max_bounds)
{
    // min bounds
    for (int i=0; i<min_bounds.size(); ++i) {
        const std::pair<double, std::vector<double> >& bound = min_bounds[i];
        double min_bound = bound.first;
        std::vector<double> bound_vals = bound.second;
        ZoneControl zc(bound_vals);
        zc.AddControl(ZoneControl::SUM,
                      ZoneControl::MORE_THAN, min_bound);
        controllers.push_back(zc);
    }
    // max bounds
    for (int i=0; i<max_bounds.size(); ++i) {
        const std::pair<double, std::vector<double> >& bound = max_bounds[i];
        double max_bound = bound.first;
        std::vector<double> bound_vals = bound.second;
        ZoneControl zc(bound_vals);
        zc.AddControl(ZoneControl::SUM,
                      ZoneControl::LESS_THAN, max_bound);
        controllers.push_back(zc);
    }
}

const std::vector<std::vector<int> > maxp_wrapper::GetClusters() {
    return cluster_ids;
}

MaxpRegion* maxp_wrapper::RunMaxp() {
    MaxpRegion* maxp = new MaxpRegion(iterations, gal, input_data, dm, num_obs, n_cols, controllers, inits, init_regions,
            rnd_seed);

    return maxp;
}