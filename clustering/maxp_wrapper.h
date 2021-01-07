//
// Created by Xun Li on 9/26/19.
//

#ifndef GEODA_MAXP_WRAPPER_H
#define GEODA_MAXP_WRAPPER_H


#include <vector>

class GeoDa;
class GalElement;
class GeoDaWeight;
class ZoneControl;
class MaxpRegion;
class RawDistMatrix;

class maxp_wrapper {
public:
    maxp_wrapper(GeoDaWeight *w,
                 const std::vector<std::vector<double> >& data,
                 int iterations,
                 int inits,
                 const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                 const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                 const std::vector<int>& init_regions,
                 const std::string &distance_method,
                 int rnd_seed);

    virtual ~maxp_wrapper();

    virtual const std::vector<std::vector<int> > GetClusters();

    virtual MaxpRegion* RunMaxp();

    virtual void CreateController(const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                  const std::vector<std::pair<double, std::vector<double> > >& max_bounds);

protected:
    int num_obs;

    int n_cols;

    int iterations;

    int inits;

    GalElement *gal;

    double **input_data;

    RawDistMatrix *dm;

    std::vector<ZoneControl> controllers;

    std::vector<int> init_regions;

    int rnd_seed;

    std::vector<std::vector<int> > cluster_ids;
};


#endif //GEODA_MAXP_WRAPPER_H
