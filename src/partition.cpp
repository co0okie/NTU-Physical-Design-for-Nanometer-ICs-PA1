#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <cmath>
#include <list>
#include <algorithm>
#include <cassert>
#include <array>
#include <unordered_set>
#include <chrono>
#include <numeric>
#include <exception>
#include <random>
#include <tuple>
#include "partition.h"
#include "clusterning.h"

constexpr int PARTITION_REPEAT = 16;

using std::vector; 
using std::string; 
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::unordered_map;
using std::unordered_set;
using std::list;
using std::tuple;
using std::max;
using std::min;
using std::array;
using std::chrono::high_resolution_clock;

typedef std::chrono::high_resolution_clock::time_point time_point;

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

size_t get_cut_size(const vector<array<size_t, 2>>& cell_count_of_net_by_group) {
    size_t cut_size = 0;
    for (net_t net = 0; net < cell_count_of_net_by_group.size(); net++) {
        if (cell_count_of_net_by_group[net][0] != 0 && cell_count_of_net_by_group[net][1] != 0) 
            cut_size++;
    }
    return cut_size;
}

class FiducciaMattheysesPartitioner {
public:
    FiducciaMattheysesPartitioner(
        double _balance_degree, 
        const vector<vector<cell_t>>& _cells_of_net, const vector<vector<net_t>>& _nets_of_cell,
        const vector<bool>& _group_of_cell, const vector<size_t>& _size_of_cell,
        time_point _start_time, double _time_limit_seconds
    );
    // return 0: success
    // return 1: incorrect cut size after trace back, which should not happen
    // return 2: no valid partitioning found
    int solve();
    uint64_t best_cut_size;
    vector<bool> group_of_cell; // cell -> group index (0 or 1)

private:
    // return 0: success
    // return 1: cut size does not decrease
    // return 2: no valid partitioning found in this iteration
    int do_one_iteration(size_t& cut_size, size_t& iteration_min_cut_size);

    // return 0: success
    // return 1: no cell can be moved, ending the iteration
    int move_one_cell(size_t& cut_size, size_t& min_cut_size, size_t& picked_No_after_min_cut, 
        array<size_t, 2>& unlocked_size_of_group);
    
    const time_point start_time;
    const double time_limit_seconds;

    const double balance_degree;
    const double min_group_size;
    const double max_group_size;
    bool is_balanced() const {
        return min_group_size <= size_of_group[0] && size_of_group[0] <= max_group_size;
    }
    const vector<size_t>& size_of_cell; // cell -> size of this cell, used for clustering

    const gain_t gain_offset; // max possible gain = max number of nets connected to a cell

    size_t num_of_nets() { return cells_of_net.size(); }
    size_t num_of_cells() { return nets_of_cell.size(); }

    const vector<vector<cell_t>>& cells_of_net; // net -> cell indices
    const vector<vector<net_t>>& nets_of_cell; // cell -> net indices

    vector<cell_t> cell_of_picked_No; // cell picked order sequence -> cell index
    uint32_t gain_update_ts;

    vector<gain_t> gain_of_cell; // cell -> gain of this cell
    vector<uint32_t> gain_update_ts_of_cell; // cell -> timestamp of the last gain update of this cell, to choose the latest updated cell
    vector<bool> locked_of_cell; // cell -> {0: unlocked, 1: locked}
    void increase_gain_of_cell(cell_t cell, gain_t gain_change);

    vector<array<size_t, 2>> cell_count_of_net_by_group; // net -> group -> number of cells in this group connected to this net

    array<vector<cell_t>, 2> head_cell_of_group_by_offset_gain; // group -> gain -> gain list head cell with this gain in this group
    cell_t& head_cell_of_group_by_gain(bool group, gain_t gain) {
        return head_cell_of_group_by_offset_gain[group][gain + gain_offset];
    }
    vector<cell_t> prev_of_cell; // cell -> previous cell in the same gain list
    vector<cell_t> next_of_cell; // cell -> next cell in the same gain list
    void push_cell_to_front_of_gain(cell_t cell, gain_t gain, bool group);
    void erase_cell_of_gain(cell_t cell, gain_t gain, bool group);

    array<size_t, 2> size_of_group; // group -> cell size sum in this group
    array<gain_t, 2> max_gain_of_group;
};

tuple<size_t, vector<bool>> hMetis(double balance_degree, const vector<vector<cell_t>>& cells_of_net, 
    const vector<vector<net_t>>& nets_of_cell, const vector<size_t>& size_of_cell,
    time_point start_time, double time_limit_seconds, int depth = 1);

PartitionOutput partition(
    const PartitionInput& input, 
    uint32_t repeat,
    double time_limit_seconds
) {
    auto start_time = high_resolution_clock::now();

    vector<vector<cell_t>> cells_of_net;
    vector<vector<net_t>> nets_of_cell;
    cells_of_net.reserve(input.cells_of_net.size());
    for (net_t net = 0; net < input.cells_of_net.size(); net++) {
        cells_of_net.emplace_back(BEGIN_END(input.cells_of_net[net]));
        for (cell_t cell : cells_of_net[net]) {
            nets_of_cell.resize(max(nets_of_cell.size(), (size_t) cell + 1));
            nets_of_cell[cell].push_back(net);
        }
    }

    vector<size_t> size_of_cell(nets_of_cell.size(), 1); // cell -> size of this cell, used for clustering

    LOG(INFO) << "start hMetis with " << cells_of_net.size() << " nets " 
        << nets_of_cell.size() << " cells." << endl;
    
    auto result = hMetis(input.balance_degree, cells_of_net, nets_of_cell, size_of_cell, 
        start_time, time_limit_seconds);
    auto& cut_size = std::get<0>(result);
    auto& group_of_cell = std::get<1>(result);

    if (group_of_cell.empty()) {
        LOG(WARNING) << "hMetis fails to find a valid partitioning, randomly partitioning" << endl;
        group_of_cell.resize(nets_of_cell.size());
        for (cell_t cell = 0; cell < nets_of_cell.size(); cell++) {
            group_of_cell[cell] = cell % 2;
        }
        vector<array<size_t, 2>> cell_count_of_net_by_group(cells_of_net.size(), {0, 0});
        for (cell_t cell = 0; cell < nets_of_cell.size(); cell++) {
            for (net_t net : nets_of_cell[cell]) {
                cell_count_of_net_by_group[net][group_of_cell[cell]]++;
            }
        }
        cut_size = get_cut_size(cell_count_of_net_by_group);
    }

    PartitionOutput output;
    output.cut_size = cut_size;
    for (cell_t cell = 0; cell < nets_of_cell.size(); cell++) {
        output.cells_of_group[group_of_cell[cell]].push_back(cell);
    }

    LOG(INFO) << "final cut size = " << cut_size << endl;
    LOG(INFO) << "group0_size " << output.cells_of_group[0].size() << " group1_size " << output.cells_of_group[1].size() << endl;

    return output;
}

tuple<size_t, vector<bool>> hMetis(
    double balance_degree,
    const vector<vector<cell_t>>& cells_of_net, 
    const vector<vector<net_t>>& nets_of_cell,
    const vector<size_t>& size_of_cell,
    time_point start_time, double time_limit_seconds, int depth
) {
    LOG(INFO) << "*================================== depth " << depth << 
        " hMetis ==================================* " << endl;
    tuple<size_t, vector<bool>> best_result = {UINT64_MAX, {}};

    LOG(INFO) << "trying pre-partitioning..." << endl;

    FiducciaMattheysesPartitioner pre_partitioner(
        balance_degree, cells_of_net, nets_of_cell, {},
        size_of_cell, start_time, time_limit_seconds
    );
    int result = pre_partitioner.solve();

    if (result != 0) {
        LOG(INFO) << "pre-partitioning fails, maybe it's too deep" << endl;
        LOG(INFO) << "*================================== end of depth " << depth << 
            " hMetis ==================================* " << endl;
        return {UINT64_MAX, {}};
    }

    best_result = {pre_partitioner.best_cut_size, pre_partitioner.group_of_cell};

    LOG(INFO) << "try clustering..." << endl;

    auto cluster_result = first_choice_clustering(cells_of_net, nets_of_cell, size_of_cell);

    auto min_max = std::minmax_element(BEGIN_END(cluster_result.size_of_cluster));
    LOG(INFO) << "after clustering: " << cluster_result.clusters_of_new_net.size() << " nets "
        << cluster_result.new_nets_of_cluster.size() << " clusters" << endl;
    LOG(INFO) << "min cluster size " << *min_max.first << " max cluster size " << *min_max.second << endl;
    
    auto hMetis_result = hMetis(
        balance_degree, 
        cluster_result.clusters_of_new_net, cluster_result.new_nets_of_cluster, 
        cluster_result.size_of_cluster, start_time, time_limit_seconds, depth + 1
    );

    if (std::get<1>(hMetis_result).empty()) {
        LOG(INFO) << "problem size small enough to directly partitioning" << endl;
        for (uint32_t r = 1; r <= PARTITION_REPEAT; r++) {
            LOG(INFO) << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX repeat " 
                << r << "/" << PARTITION_REPEAT << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
            FiducciaMattheysesPartitioner partitioner(
                balance_degree, cells_of_net, nets_of_cell, 
                {}, // random initial partition
                size_of_cell, start_time, time_limit_seconds
            );
            partitioner.solve();

            if (partitioner.best_cut_size < std::get<0>(best_result)) {
                best_result = {partitioner.best_cut_size, partitioner.group_of_cell};
            }
        }

        LOG(INFO) << "best cut size = " << std::get<0>(best_result) 
            << " after " << PARTITION_REPEAT << " repeat" << endl;
    } else {
        auto& group_of_cluster = std::get<1>(hMetis_result);
        vector<bool> group_of_cell(nets_of_cell.size());
        for (cell_t cell = 0; cell < nets_of_cell.size(); cell++) {
            group_of_cell[cell] = group_of_cluster[cluster_result.cluster_of_cell[cell]];
        }

        LOG(INFO) << "start refinement for depth " << depth << " hMetis" << endl;

        FiducciaMattheysesPartitioner partitioner(
            balance_degree, cells_of_net, nets_of_cell, group_of_cell,
            size_of_cell, start_time, time_limit_seconds
        );
        int result = partitioner.solve();

        best_result = {partitioner.best_cut_size, partitioner.group_of_cell};
    }

    
    LOG(INFO) << "*================================== end of depth " << depth << 
        " hMetis ==================================* " << endl;
    return best_result;
}

FiducciaMattheysesPartitioner::FiducciaMattheysesPartitioner(
    double _balance_degree, 
    const vector<vector<cell_t>>& _cells_of_net, const vector<vector<net_t>>& _nets_of_cell,
    const vector<bool>& _group_of_cell, const vector<size_t>& _size_of_cell,
    time_point _start_time, double _time_limit_seconds
) : balance_degree(_balance_degree), size_of_cell(_size_of_cell),
    min_group_size((1 - _balance_degree) / 2 * std::accumulate(BEGIN_END(_size_of_cell), 0)), 
    max_group_size((1 + _balance_degree) / 2 * std::accumulate(BEGIN_END(_size_of_cell), 0)),
    cells_of_net(_cells_of_net), nets_of_cell(_nets_of_cell), 
    start_time(_start_time), time_limit_seconds(_time_limit_seconds),
    gain_offset(std::max_element(BEGIN_END(_nets_of_cell), [](const vector<net_t>& a, const vector<net_t>& b) {
        return a.size() < b.size();
    })->size()) 
{
    cell_of_picked_No.reserve(num_of_cells());
    gain_of_cell.reserve(num_of_cells());
    gain_update_ts_of_cell.reserve(num_of_cells());
    locked_of_cell.reserve(num_of_cells());
    cell_count_of_net_by_group.reserve(num_of_nets());
    head_cell_of_group_by_offset_gain[0].reserve(2 * gain_offset + 1);
    head_cell_of_group_by_offset_gain[1].reserve(2 * gain_offset + 1);
    prev_of_cell.reserve(num_of_cells());
    next_of_cell.reserve(num_of_cells());

    LOG(INFO) << "balance_degree: " << balance_degree << endl;
    LOG(INFO) << "\nmax possible gain: " << gain_offset << endl;
    LOG(INFO) << "\nnum_of_nets: " << num_of_nets() << endl;
    for (net_t net = 0; net < num_of_nets(); ++net) {
        LOG(DEBUG) << "n" << net << "{";
        for (size_t j = 0; j < cells_of_net[net].size(); ++j) {
            if (j) { LOG(DEBUG) << " "; }
            LOG(DEBUG) << "c" << cells_of_net[net][j];
        }
        LOG(DEBUG) << "}\n";
    }
    LOG(DEBUG) << endl;
    LOG(INFO) << "\nnum_of_cells: " << num_of_cells() << endl;
    for (cell_t cell = 0; cell < num_of_cells(); ++cell) {
        LOG(DEBUG) << "c" << cell << "{";
        for (size_t i = 0; i < nets_of_cell[cell].size(); ++i) {
            if (i) { LOG(DEBUG) << " "; }
            LOG(DEBUG) << "n" << nets_of_cell[cell][i];
        }
        LOG(DEBUG) << "}\n";
    }
    LOG(DEBUG) << endl;
    
    if (_group_of_cell.size()) {
        LOG(INFO) << "Use the given initial partition." << endl;
        group_of_cell = _group_of_cell;
    } else {
        LOG(INFO) << "Use random initial partition." << endl;
        group_of_cell.resize(num_of_cells());

        std::vector<size_t> random_cells(num_of_cells());
        std::iota(BEGIN_END(random_cells), 0);
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(BEGIN_END(random_cells), g);
        
        for (cell_t cell = 0; cell < num_of_cells(); cell++) {
            group_of_cell[cell] = random_cells[cell] < num_of_cells() / 2;
        }
    }
}

int FiducciaMattheysesPartitioner::solve() {
    cell_count_of_net_by_group.assign(num_of_nets(), {0, 0});
    size_of_group = {0, 0};
    for (cell_t cell = 0; cell < num_of_cells(); cell++) {
        bool group = group_of_cell[cell];
        size_of_group[group] += size_of_cell[cell];
        for (net_t net : nets_of_cell[cell]) {
            cell_count_of_net_by_group[net][group]++;
        }
    }
    
    LOG(INFO) << "\ngroup0_size " << size_of_group[0] << " group1_size " << size_of_group[1] << endl;
    LOG(DEBUG) << "\n<cell_index>:<group_index> " << endl;
    for (cell_t cell = 0; cell < num_of_cells(); cell++) {
        LOG(DEBUG) << "c" << cell << ":" << group_of_cell[cell] << " ";
    }
    LOG(DEBUG) << endl;
    LOG(DEBUG) << "\n<net_index>:<group 0 cell count>/<group 1 cell count> " << endl;
    for (net_t net = 0; net < num_of_nets(); net++) {
        LOG(DEBUG) << "n" << net << ":" << cell_count_of_net_by_group[net][0] 
            << "/" << cell_count_of_net_by_group[net][1] << " ";
    }
    LOG(DEBUG) << endl;

    uint64_t cut_size = get_cut_size(cell_count_of_net_by_group); // current cut size, not always valid
    // record only when satisfying the balance constraint
    uint64_t iteration_min_cut_size = is_balanced() ? cut_size : UINT64_MAX;
    LOG(INFO) << "\ninitial cut size: " << cut_size << endl;
    if (!is_balanced()) LOG(INFO) << "initial partition does not satisfy the balance constraint, initial cut size is invalid." << endl;
    for (size_t iter = 0;; iter++) {
        LOG(INFO) << "================================ iteration " << iter + 1 
            << " ================================" << endl;

        int result = do_one_iteration(cut_size, iteration_min_cut_size);
        if (result == 1) break;
        else if (result == 2) return 2;
        
        if (time_limit_seconds > 0) {
            std::chrono::duration<double> elapsed = high_resolution_clock::now() - start_time;
            if (elapsed.count() > time_limit_seconds) break;
        }
    }

    cut_size = get_cut_size(cell_count_of_net_by_group);
    if (cut_size != iteration_min_cut_size) {
        LOG(FETAL) << "Error: cut size " << iteration_min_cut_size 
            << " after trace back does not equal to the minimum cut size " << cut_size 
            << " found in this iteration." << endl;
        return 1;
    }

    best_cut_size = cut_size;
    
    return 0;
}

int FiducciaMattheysesPartitioner::do_one_iteration(size_t& cut_size, size_t& iteration_min_cut_size) {
    max_gain_of_group = {MIN_GAIN, MIN_GAIN};
    head_cell_of_group_by_offset_gain[0].assign(2 * gain_offset + 1, -1);
    head_cell_of_group_by_offset_gain[1].assign(2 * gain_offset + 1, -1);
    prev_of_cell.assign(num_of_cells(), -1);
    next_of_cell.assign(num_of_cells(), -1);
    gain_of_cell.assign(num_of_cells(), 0);
    gain_update_ts = 0;
    gain_update_ts_of_cell.assign(num_of_cells(), gain_update_ts);
    
    for (cell_t cell = 0; cell < num_of_cells(); cell++) {
        gain_t gain = 0;
        bool group = group_of_cell[cell];
        for (size_t j = 0; j < nets_of_cell[cell].size(); j++) {
            const net_t net = nets_of_cell[cell][j];
            size_t cell_count = cell_count_of_net_by_group[net][group];
            if (cell_count == 1) gain++;
            if (cell_count == cells_of_net[net].size()) gain--;
        }
        push_cell_to_front_of_gain(cell, gain, group);
        max_gain_of_group[group] = max(max_gain_of_group[group], gain);
        gain_of_cell[cell] = gain;
    }
    
    RUN(DEBUG) {
        for (int group = 0; group < 2; group++) {
            LOG(DEBUG) << "group" << group << " gain table: " << endl;
            size_t cell_count = size_of_group[group];
            for (gain_t gain = max_gain_of_group[group]; cell_count; gain--) {
                LOG(DEBUG) << "gain " << gain << ": ";
                for (cell_t cell = head_cell_of_group_by_gain(group, gain); cell != -1; cell = next_of_cell[cell]) { 
                    LOG(DEBUG) << "c" << cell << " ";
                    cell_count -= size_of_cell[cell];
                }
                LOG(DEBUG) << endl;
            }
        }
    }
    
    size_t min_cut_size = iteration_min_cut_size;
    cell_of_picked_No.clear();
    cell_of_picked_No.reserve(num_of_cells());
    size_t picked_No_after_min_cut = 0; // {cut 5} {pick cell 1} {cut 2} {pick cell 0} {cut 6}, in this case picked_No_after_min_cut = 1
    array<size_t, 2> unlocked_size_of_group = {size_of_group[0], size_of_group[1]}; // group -> number of unlocked cells in this group
    auto unlocked_total_size = [&]() { return unlocked_size_of_group[0] + unlocked_size_of_group[1]; };
    locked_of_cell.assign(num_of_cells(), false);

    for (size_t i = 1 ; unlocked_total_size(); i++) {
        LOG(DEBUG) << "---------------------------------- move " << i << "th cell ----------------------------------" << endl;

        int result = move_one_cell(cut_size, min_cut_size, picked_No_after_min_cut, unlocked_size_of_group);
        if (result == 1) break;
    }
    
    if (min_cut_size == UINT64_MAX) {
        LOG(INFO) << "no valid partitioning found in this iteration, ending the loop" << endl;
        return 2;
    }
    LOG(INFO) << "min_cut_size = " << min_cut_size << " after " 
        << picked_No_after_min_cut << " cells moved" << endl;
    LOG(DEBUG) << "order of picked cells: ";
    RUN(DEBUG) for (size_t i = 0; i < cell_of_picked_No.size(); i++) {
        LOG(DEBUG) << "c" << cell_of_picked_No[i] << " ";
    }
    LOG(DEBUG) << endl;
    // trace back to the valid state with min cut size
    for (size_t i = cell_of_picked_No.size(); i --> picked_No_after_min_cut;) {
        cell_t cell = cell_of_picked_No[i];
        bool old_group = group_of_cell[cell], new_group = !old_group;
        group_of_cell[cell] = new_group;
        unlocked_size_of_group[new_group] += size_of_cell[cell];
        size_of_group[old_group] -= size_of_cell[cell];
        size_of_group[new_group] += size_of_cell[cell];
        for (net_t net : nets_of_cell[cell]) {
            cell_count_of_net_by_group[net][old_group]--;
            cell_count_of_net_by_group[net][new_group]++;
        }
    }
    cut_size = min_cut_size;

    LOG(INFO) << "after trace back: group0_size " << size_of_group[0] 
        << " group1_size " << size_of_group[1] << endl;
    for (cell_t cell = 0; cell < num_of_cells(); cell++) {
        LOG(DEBUG) << "c" << cell << ":" << group_of_cell[cell] << " ";
    }
    LOG(DEBUG) << endl;

    if (min_cut_size < iteration_min_cut_size) {
        iteration_min_cut_size = min_cut_size;
        return 0;
    } else {
        LOG(INFO) << "cut size does not decrease, ending the loop." << endl;
        return 1;
    }
}

int FiducciaMattheysesPartitioner::move_one_cell(
    size_t& cut_size, size_t& min_cut_size, size_t& picked_No_after_min_cut,
    array<size_t, 2>& unlocked_size_of_group
) {
    while (unlocked_size_of_group[0] && head_cell_of_group_by_gain(0, max_gain_of_group[0]) == -1) max_gain_of_group[0]--;
    while (unlocked_size_of_group[1] && head_cell_of_group_by_gain(1, max_gain_of_group[1]) == -1) max_gain_of_group[1]--;

    bool permit_group_1_to_0 = size_of_group[0] <= max_group_size;
    bool permit_group_0_to_1 = size_of_group[0] >= min_group_size;

    // pick a cell to move
    cell_t cell_to_move;
    if ((permit_group_1_to_0 ^ permit_group_0_to_1) || 
    unlocked_size_of_group[0] == 0 || unlocked_size_of_group[1] == 0) {
        bool src = unlocked_size_of_group[0] == 0 ? 1 : unlocked_size_of_group[1] == 0 ? 0 : 
            permit_group_1_to_0 ? 1 : 0;
        if (unlocked_size_of_group[src] == 0) {
            LOG(INFO) << "Fail to find a cell in group " << src << ", ending the loop." << endl;
            return 1;
        }
        LOG(DEBUG) << "pick a cell from group " << src << endl;
        gain_t& max_gain = max_gain_of_group[src];
        cell_to_move = head_cell_of_group_by_gain(src, max_gain);
    } else { // pick a cell in either group
        cell_t cell0 = head_cell_of_group_by_gain(0, max_gain_of_group[0]);
        cell_t cell1 = head_cell_of_group_by_gain(1, max_gain_of_group[1]);
        if (max_gain_of_group[0] == max_gain_of_group[1]) {
            if (gain_update_ts_of_cell[cell0] == gain_update_ts_of_cell[cell1]) {
                cell_to_move = size_of_group[0] > size_of_group[1] ? cell0 : cell1; // for balance
            } else {
                cell_to_move = gain_update_ts_of_cell[cell0] > gain_update_ts_of_cell[cell1] ? cell0 : cell1;
            }
        } else {
            cell_to_move = max_gain_of_group[1] > max_gain_of_group[0] ? cell1 : cell0;
        }
    }

    LOG(DEBUG) << "group0_max_gain " << max_gain_of_group[0] 
        << " group1_max_gain " << max_gain_of_group[1] << endl;
    LOG(DEBUG) << "move c" << cell_to_move << " from " 
        << group_of_cell[cell_to_move] << " to " << !group_of_cell[cell_to_move]
        << " gain " << gain_of_cell[cell_to_move] << " cut_size " << cut_size
        << " -> " << cut_size - gain_of_cell[cell_to_move] << endl;

    gain_update_ts++;

    // update the gain of the affected cells
    for (size_t i = 0; i < nets_of_cell[cell_to_move].size(); i++) {
        net_t net = nets_of_cell[cell_to_move][i];

        bool src_group = group_of_cell[cell_to_move], dst_group = !src_group;
        size_t dst_group_cell_count = cell_count_of_net_by_group[net][dst_group];
        size_t src_group_cell_count = cell_count_of_net_by_group[net][src_group];

        if (dst_group_cell_count == 0) {
            for (cell_t cell : cells_of_net[net]) {
                if (!locked_of_cell[cell] && cell != cell_to_move) increase_gain_of_cell(cell, 1);
            }
        } else if (dst_group_cell_count == 1) {
            for (cell_t cell : cells_of_net[net]) {
                if (!locked_of_cell[cell] && cell != cell_to_move && group_of_cell[cell] == dst_group) {
                    increase_gain_of_cell(cell, -1);
                    break;
                }
            }
        }

        dst_group_cell_count++;
        src_group_cell_count--;

        if (src_group_cell_count == 0) {
            for (cell_t cell : cells_of_net[net]) {
                if (!locked_of_cell[cell] && cell != cell_to_move) increase_gain_of_cell(cell, -1);
            }
        } else if (src_group_cell_count == 1) {
            for (cell_t cell : cells_of_net[net]) {
                if (!locked_of_cell[cell] && cell != cell_to_move && group_of_cell[cell] == src_group) {
                    increase_gain_of_cell(cell, 1);
                    break;
                }
            }
        }
    }

    // move cell to another group
    bool old_group = group_of_cell[cell_to_move], new_group = !old_group;
    gain_t gain = gain_of_cell[cell_to_move];
    unlocked_size_of_group[old_group] -= size_of_cell[cell_to_move];
    size_of_group[old_group] -= size_of_cell[cell_to_move];
    erase_cell_of_gain(cell_to_move, gain, old_group);
    group_of_cell[cell_to_move] = new_group;
    size_of_group[new_group] += size_of_cell[cell_to_move];
    for (net_t net : nets_of_cell[cell_to_move]) {
        cell_count_of_net_by_group[net][old_group]--;
        cell_count_of_net_by_group[net][new_group]++;
    }
    cut_size -= gain_of_cell[cell_to_move];
    cell_of_picked_No.push_back(cell_to_move);
    // if not satisfy balance constraint, consider invalid
    if (cut_size < min_cut_size && is_balanced()) {
        min_cut_size = cut_size; // min_cut_size must be valid
        picked_No_after_min_cut = cell_of_picked_No.size();
    }
    locked_of_cell[cell_to_move] = true;

    RUN(TRACE) {
        for (int group = 0; group < 2; group++) {
            LOG(TRACE) << "group" << group << " gain table: " << endl;
            size_t cell_count = unlocked_size_of_group[group];
            for (gain_t gain = max_gain_of_group[group]; cell_count; gain--) {
                LOG(TRACE) << "gain " << gain << ": ";
                for (cell_t cell = head_cell_of_group_by_gain(group, gain); cell != -1; cell = next_of_cell[cell]) { 
                    LOG(TRACE) << "c" << cell << " ";
                    cell_count -= size_of_cell[cell];
                }
                LOG(TRACE) << endl;
            }
        }
    }
    
    return 0;
}

void FiducciaMattheysesPartitioner::increase_gain_of_cell(cell_t cell, gain_t gain_change) {
    if (gain_change == 0) return;
    bool group = group_of_cell[cell];
    gain_t old_gain = gain_of_cell[cell], new_gain = old_gain + gain_change;
    erase_cell_of_gain(cell, old_gain, group);
    gain_of_cell[cell] = new_gain;
    push_cell_to_front_of_gain(cell, new_gain, group);
    max_gain_of_group[group] = max(max_gain_of_group[group], new_gain);
    gain_update_ts_of_cell[cell] = gain_update_ts;
}

void FiducciaMattheysesPartitioner::push_cell_to_front_of_gain(cell_t cell, gain_t gain, bool group) {
    cell_t head_cell = head_cell_of_group_by_gain(group, gain);
    prev_of_cell[cell] = -1;
    next_of_cell[cell] = head_cell;
    if (head_cell != -1) prev_of_cell[head_cell] = cell;
    head_cell_of_group_by_gain(group, gain) = cell;
    gain_update_ts_of_cell[cell] = gain_update_ts;
}

void FiducciaMattheysesPartitioner::erase_cell_of_gain(cell_t cell, gain_t gain, bool group) {
    cell_t p = prev_of_cell[cell], n = next_of_cell[cell];
    if (p != -1) {
        next_of_cell[p] = n;
    } else {
        head_cell_of_group_by_gain(group, gain) = n;
    }
    if (n != -1) prev_of_cell[n] = p;
    prev_of_cell[cell] = -1;
    next_of_cell[cell] = -1;
}

std::tuple<
    PartitionInput, 
    std::vector<std::string> // name_of_cell: cell -> cell name
> PartitionInput::from_file(const std::string& input_file) {
    std::ifstream input(input_file);
    if (!input.is_open()) {
        LOG(FETAL) << "Error: Could not open file " << input_file << endl;
        exit(1);
    }

    PartitionInput partition_input;
    input >> partition_input.balance_degree;

    std::vector<std::string> name_of_cell;
    std::unordered_map<std::string, cell_t> cell_of_name;
    
    std::string tmp;
    while (input >> tmp) {
        input >> tmp; // net name
        partition_input.cells_of_net.push_back({});
        for (;;) {
            input >> tmp;
            if (tmp == ";") break;

            if (!cell_of_name.count(tmp)) {
                cell_of_name[tmp] = name_of_cell.size();
                name_of_cell.push_back(tmp);
            }
            cell_t cell = cell_of_name[tmp];
            partition_input.cells_of_net.back().insert(cell);
        }
    }

    return {partition_input, name_of_cell};
}

void PartitionOutput::to_file(const std::string& output_file, const std::vector<std::string>& name_of_cell) const {
    std::ofstream output(output_file);
    if (!output.is_open()) {
        LOG(FETAL) << "Error: Could not open file " << output_file << endl;
        exit(1);
    }

    output << "Cutsize = " << cut_size << endl;
    for (int group = 0; group < 2; group++) {
        output << "G" << group + 1 << " " << cells_of_group[group].size() << endl;
        for (cell_t cell : cells_of_group[group]) {
            output << name_of_cell[cell] << " ";
        }
        output << ";" << endl;
    }
}