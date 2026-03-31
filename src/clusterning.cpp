#include "clusterning.h"
#include <numeric>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <iostream>

using std::vector;
using std::iota;
using std::unordered_map;
using std::unordered_set;
using std::cout;
using std::endl;

enum LogLevel {
    QUIET = 0,
    FETAL = 1,
    ERROR = 2,
    WARNING = 3,
    INFO = 4,
    DEBUG = 5,
    TRACE = 6
};
constexpr LogLevel LOG_LEVEL = INFO;
#define RUN(level) if (level > LOG_LEVEL) {} else
#define LOG(level) RUN(level) std::cout

ClusteringResult first_choice_clustering(
    const vector<vector<cell_t>>& cells_of_net, 
    const vector<vector<net_t>>& nets_of_cell,
    const vector<size_t>& size_of_cell
) {
    auto num_of_cells = nets_of_cell.size();
    auto num_of_nets = cells_of_net.size();
    vector<cell_t> root_of_cell(num_of_cells); // cell -> root cell represent a unique cluster
    iota(BEGIN_END(root_of_cell), 0);
    vector<size_t> size_of_root(BEGIN_END(size_of_cell)); // root cell -> total cell count of this cluster
    auto is_matched = [&](cell_t cell) { return size_of_root[root_of_cell[cell]] > size_of_cell[cell]; };

    for (cell_t cell = 0; cell < num_of_cells; cell++) {
        if (is_matched(cell)) {
            LOG(TRACE) << "Cell c" << cell << " is already matched to cluster c" << root_of_cell[cell] << ", skipping" << endl;
            continue;
        }
        LOG(TRACE) << "Evaluating cell c" << cell << " for clustering:" << endl;

        double max_weight = -1.0;
        cell_t best_root = -1;
        unordered_map<cell_t, double> weight_of_root;

        // calculate the weight of each possible cluster
        for (net_t net : nets_of_cell[cell]) {
            // prevent net with too many cells to save time
            // their contribution to weight is small enough to ignore
            if (cells_of_net[net].size() > 64) continue; 
            
            double edge_weight = 1.0 / (cells_of_net[net].size() - 1);
            for (cell_t neighbor : cells_of_net[net]) {
                if (neighbor == cell) continue;
                weight_of_root[root_of_cell[neighbor]] += edge_weight;
            }
        }

        // find the best cluster for this cell
        for (const auto& pair : weight_of_root) {
            cell_t root = pair.first;
            double weight = pair.second;
            LOG(TRACE) << "  - Cluster c" << root << " weight: " << weight << ", size: " << size_of_root[root] << endl;
            if (weight > max_weight ||
            (weight == max_weight && size_of_root[root] < size_of_root[best_root])) {
                max_weight = weight;
                best_root = root;
            }
        }

        // add cell to cluster
        if (best_root != -1) {
            root_of_cell[cell] = best_root;
            size_of_root[best_root] += size_of_cell[cell];
            LOG(TRACE) << "c" << cell << " clustered with c" << best_root 
                << " (weight: " << max_weight << ", size: " << size_of_root[best_root] << ")" << endl;
        }
    }

    unordered_map<cell_t, cell_t> cluster_of_root; // root cell -> cluster index
    vector<cell_t> root_of_cluster; // cluster -> root cell represent this cluster
    vector<unordered_set<cell_t>> clusters_of_new_net; // net -> cluster indices connected to this net
    for (net_t net = 0; net < num_of_nets; net++) {
        unordered_set<cell_t> clusters;
        for (cell_t cell : cells_of_net[net]) {
            cell_t root = root_of_cell[cell];
            if (!cluster_of_root.count(root)) {
                cluster_of_root[root] = root_of_cluster.size();
                root_of_cluster.emplace_back(root);
            }
            cell_t cluster = cluster_of_root[root];
            clusters.insert(cluster);
        }
        if (clusters.size() > 1) {
            clusters_of_new_net.emplace_back(std::move(clusters));
        } else if (clusters.size() == 1) {
            LOG(TRACE) << "Net n" << net << " is fully absorbed into cluster " << *clusters.begin() << endl;
        }
    }

    vector<vector<net_t>> new_nets_of_cluster(root_of_cluster.size()); // cluster -> net indices connected to this cluster
    for (net_t new_net = 0; new_net < clusters_of_new_net.size(); new_net++) {
        for (cell_t cluster : clusters_of_new_net[new_net]) {
            new_nets_of_cluster[cluster].push_back(new_net);
        }
    }

    vector<cell_t> cluster_of_cell(num_of_cells); // cell -> cluster index
    for (cell_t cell = 0; cell < num_of_cells; cell++) {
        cluster_of_cell[cell] = cluster_of_root[root_of_cell[cell]];
    }

    vector<size_t> size_of_cluster(root_of_cluster.size()); // cluster -> cell count of this cluster
    for (cell_t cluster = 0; cluster < root_of_cluster.size(); cluster++) {
        size_of_cluster[cluster] = size_of_root[root_of_cluster[cluster]];
    }

    RUN(TRACE) {
        for (cell_t cell = 0; cell < num_of_cells; cell++) {
            LOG(TRACE) << "Cell c" << cell << " -> Cluster " << cluster_of_cell[cell] << endl;
        }
    }

    vector<vector<cell_t>> _clusters_of_new_net;
    for (const auto& clusters : clusters_of_new_net) {
        _clusters_of_new_net.emplace_back(BEGIN_END(clusters));
    }

    return ClusteringResult({
        .clusters_of_new_net = std::move(_clusters_of_new_net),
        .new_nets_of_cluster = std::move(new_nets_of_cluster),
        .cluster_of_cell = std::move(cluster_of_cell),
        .size_of_cluster = std::move(size_of_cluster)
    });
}