#ifndef __JSGEODSA_GDA_CLUSTERING__
#define __JSGEODSA_GDA_CLUSTERING__

#include <vector>
#include <string>

class GeoDaWeight;

// APIs of clustering

/**
 *
 * @param w
 * @param data
 * @param bound_vals
 * @param min_bound
 * @param local_search_method
 * @param initial
 * @param tabu_length
 * @param cool_rate
 * @param seeds
 * @param distance_method
 * @param rand_seed
 * @return
 */
const std::vector<std::vector<int> > gda_maxp(GeoDaWeight *w,
                                              const std::vector<std::vector<double> > &data,
                                              int iterations,
                                              int inits,
                                              const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                              const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                              const std::vector<int>& init_regions,
                                              const std::string &distance_method,
                                              int rnd_seed);

/**
 *
 * @param k
 * @param w
 * @param data
 * @param redcap_method
 * @param distance_method
 * @param bound_vals
 * @param min_bound
 * @param rand_seed
 * @return
 */
const std::vector<std::vector<int> > gda_redcap(unsigned int k,
                                                GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &data,
                                                const std::string &redcap_method,
                                                const std::string &distance_method,
                                                const std::vector<double>& bound_vals,
                                                double min_bound,
                                                int rand_seed);

/**
 *
 * @param k
 * @param w
 * @param data
 * @param distance_method
 * @param bound_vals
 * @param min_bound
 * @param rand_seed
 * @return
 */
const std::vector<std::vector<int> > gda_skater(unsigned int k,
                                                GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &data,
                                                const std::string &distance_method,
                                                const std::vector<double>& bound_vals,
                                                double min_bound,
                                                int rand_seed);


/**
 *
 * @param vals
 * @return
 */
double gda_sumofsquares(const std::vector<double>& vals);

/**
 *
 * @param vals
 * @return
 */
double gda_totalsumofsquare(const std::vector<std::vector<double> >& vals);

/**
 *
 * @param solution
 * @param vals
 * @return
 */
double gda_withinsumofsquare(const std::vector<std::vector<int> >& solution,
                             const std::vector<std::vector<double> >& vals);

/**
 *
 * @param solution
 * @param data
 * @return
 */
double gda_betweensumofsquare(const std::vector<std::vector<int> >& solution,
                              const std::vector<std::vector<double> >& data);

#endif

