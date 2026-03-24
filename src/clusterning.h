#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <vector>
#include "constant.h"

struct ClusteringResult {
    std::vector<std::vector<cell_t>> clusters_of_new_net; // net -> cluster indices connected to this net
    std::vector<std::vector<net_t>> new_nets_of_cluster; // cluster -> net indices connected to this cluster
    std::vector<cell_t> cluster_of_cell; // cell -> cluster index
    std::vector<size_t> size_of_cluster; // cluster -> cell count of this cluster
};

ClusteringResult first_choice_clustering(
    const std::vector<std::vector<cell_t>>& cells_of_net, 
    const std::vector<std::vector<net_t>>& nets_of_cell,
    const std::vector<std::size_t>& size_of_cell
);

#endif