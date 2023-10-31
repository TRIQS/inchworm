#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>
#include <inchworm/diagram/diagram.hpp>
#include <inchworm/diagram/inclusion_exclusion.hpp>
#include <inchworm/diagram/print.hpp>
#include <inchworm/atom_diag.hpp>
#include <inchworm/u_frame.hpp>
#include <inchworm/impurity_product.hpp>
#include <inchworm/util.hpp>
#include <inchworm/interpolator.hpp>
#include "./hubbard.hpp"
#include "./integral_util_naive.hpp"

// Short alias for this namespace
namespace pt = boost::property_tree;
int main() {
    pt::ptree root;
    pt::read_json("/Users/yangyu/src/inchworm/apps/parameters.json", root);
    int n_site = root.get<int>("n_site");
    double beta = root.get<double>("cp.beta");
    gf_struct_t gf_struct;
    for (pt::ptree::value_type &g_s: root.get_child("cp.gf_struct"))
    {
        std::string name = g_s.first;
        int size = g_s.second.get_value<int>();
        gf_struct.push_back(std::make_pair(name, size));
    }
    std::vector<double> epsilon;
    for (pt::ptree::value_type &ep: root.get_child("epsilon"))
    {
        epsilon.push_back(ep.second.get_value<double>());
    }
    std::vector<std::vector<double>> theta;
    for (pt::ptree::value_type &th: root.get_child("theta"))
    {
        std::vector<double> theta_i;
        for (pt::ptree::value_type &th_i: th.second)
        {
            theta_i.push_back(th_i.second.get_value<double>());
        }
        theta.push_back(theta_i);
    }
    std::cout << beta << std::endl;
    std::cout << n_site << std::endl;
    for(auto const& g_s: gf_struct)
    {
        std::cout << g_s.first << std::endl;
        std::cout << g_s.second << std::endl;
    }
    for(auto const& t: epsilon)
    {
        std::cout << t << std::endl;
    }
    for(auto const& t: theta)
    {
        for(auto const& t_i: t)
        {
            std::cout << t_i << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}
