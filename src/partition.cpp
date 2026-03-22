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
#include <set>
#include <chrono>
#include <numeric>
#include <exception>
#include "partition.h"

using std::vector; 
using std::string; 
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::unordered_map;
using std::list;
using std::max;
using std::min;
using std::array;
using std::chrono::high_resolution_clock;

#define BEGIN_END(container) (container).begin(), (container).end()

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

class FiducciaMattheysesPartitioner {
public:
    FiducciaMattheysesPartitioner(const PartitionInput& input, uint32_t repeat, double time_limit_seconds);
    void solve();
    PartitionOutput output;

private:
    PartitionOutput fiduccia_mattheyses();
    bool do_one_iteration(size_t& iteration_min_cut_size);
    bool move_one_cell(size_t& cut_size, size_t& min_cut_size, size_t& picked_No_after_min_cut, 
        array<size_t, 2>& unlocked_size_of_group);
    
    size_t get_cut_size() const;

    uint32_t repeat;
    high_resolution_clock::time_point start_time;
    double time_limit_seconds;

    double balance_degree;

    gain_t gain_offset; // max possible gain = max number of nets connected to a cell

    size_t num_of_nets() { return cells_of_net.size(); }
    size_t num_of_cells() { return nets_of_cell.size(); }

    vector<vector<cell_t>> cells_of_net; // net -> cell indices
    vector<vector<net_t>> nets_of_cell; // cell -> net indices

    vector<cell_t> cell_of_picked_No; // cell picked order sequence -> cell index
    uint32_t gain_update_ts;

    vector<bool> group_of_cell; // cell -> group index (0 or 1)
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

    array<size_t, 2> size_of_group; // group -> number of cells in this group
    array<gain_t, 2> max_gain_of_group;
};

FiducciaMattheysesPartitioner::FiducciaMattheysesPartitioner(
    const PartitionInput& input, uint32_t repeat, double time_limit_seconds
) {
    this->repeat = repeat;
    start_time = high_resolution_clock::now();
    this->time_limit_seconds = time_limit_seconds;
    balance_degree = input.balance_degree;
    
    gain_offset = 0;
    for (net_t net = 0; net < input.cells_of_net.size(); net++) {
        cells_of_net.push_back(vector<cell_t>(BEGIN_END(input.cells_of_net[net])));
        for (cell_t cell : cells_of_net[net]) {
            nets_of_cell.resize(max(nets_of_cell.size(), (size_t) cell + 1));
            nets_of_cell[cell].push_back(net);
        }
    }

    gain_offset = std::max_element(BEGIN_END(nets_of_cell), [](const vector<net_t>& a, const vector<net_t>& b) {
        return a.size() < b.size();
    })->size();

    cell_of_picked_No.reserve(num_of_cells());
    group_of_cell.reserve(num_of_cells());
    gain_of_cell.reserve(num_of_cells());
    gain_update_ts_of_cell.reserve(num_of_cells());
    locked_of_cell.reserve(num_of_cells());
    cell_count_of_net_by_group.reserve(num_of_nets());
    head_cell_of_group_by_offset_gain[0].reserve(2 * gain_offset + 1);
    head_cell_of_group_by_offset_gain[1].reserve(2 * gain_offset + 1);
    prev_of_cell.reserve(num_of_cells());
    next_of_cell.reserve(num_of_cells());


    LOG(INFO) << "balance_degree: " << balance_degree << endl;
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
}

void FiducciaMattheysesPartitioner::solve() {
    if ((1 - balance_degree) / 2 * num_of_cells() > num_of_cells() / 2) {
        LOG(FETAL) << "Error: balance_degree is too small, cannot satisfy the balance constraint." << endl;
        exit(1);
    }

    output.cut_size = get_cut_size();

    for (uint32_t i = 0; i < repeat; i++) {
        LOG(INFO) << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX repeat " 
            << i + 1 << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;

        auto solution = fiduccia_mattheyses();

        if (solution.cut_size < output.cut_size) {
            output = solution;
        }

        if (time_limit_seconds > 0) {
            std::chrono::duration<double> elapsed = high_resolution_clock::now() - start_time;
            if (elapsed.count() > time_limit_seconds) break;
        }
    }

    if (repeat > 1)
        LOG(INFO) << "\nfinal cut size: " << output.cut_size << endl;
}

PartitionOutput FiducciaMattheysesPartitioner::fiduccia_mattheyses() {
    group_of_cell.resize(num_of_cells());
    cell_count_of_net_by_group.assign(num_of_nets(), {0, 0});
    size_of_group = {0, 0};

    std::vector<size_t> random_cells(num_of_cells());
    std::iota(BEGIN_END(random_cells), 0);
    std::random_shuffle(BEGIN_END(random_cells));
    
    for (cell_t cell = 0; cell < num_of_cells(); cell++) {
        bool group = random_cells[cell] < num_of_cells() / 2;
        group_of_cell[cell] = group;
        size_of_group[group]++;
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

    LOG(INFO) << "\nmax possible gain: " << gain_offset << endl;

    size_t iteration_min_cut_size = get_cut_size();
    LOG(INFO) << "\ninitial cut size: " << iteration_min_cut_size << endl;
    for (size_t iter = 0;; iter++) {
        LOG(INFO) << "================================ iteration " << iter + 1 
            << " ================================" << endl;
        
        if (!do_one_iteration(iteration_min_cut_size)) {
            break;
        }
        
        if (time_limit_seconds > 0) {
            std::chrono::duration<double> elapsed = high_resolution_clock::now() - start_time;
            if (elapsed.count() > time_limit_seconds) break;
        }
    }

    PartitionOutput solution;
    solution.cut_size = get_cut_size();
    if (solution.cut_size != iteration_min_cut_size) {
        LOG(FETAL) << "Error: cut size " << iteration_min_cut_size 
            << " after trace back does not equal to the minimum cut size " << solution.cut_size 
            << " found in this iteration." << endl;
        exit(1);
    }
    for (cell_t cell = 0; cell < num_of_cells(); cell++) {
        solution.cells_of_group[group_of_cell[cell]].push_back(cell);
    }

    return solution;
}

bool FiducciaMattheysesPartitioner::do_one_iteration(size_t& iteration_min_cut_size) {
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
                    cell_count--;
                }
                LOG(DEBUG) << endl;
            }
        }
    }
    
    size_t cut_size = iteration_min_cut_size;
    size_t min_cut_size = iteration_min_cut_size;
    cell_of_picked_No.clear();
    cell_of_picked_No.reserve(num_of_cells());
    size_t picked_No_after_min_cut = 0; // {cut 5} {pick cell 1} {cut 2} {pick cell 0} {cut 6}, in this case picked_No_after_min_cut = 1
    array<size_t, 2> unlocked_size_of_group = {size_of_group[0], size_of_group[1]}; // group -> number of unlocked cells in this group
    auto num_of_unlocked_cells = [&]() { return unlocked_size_of_group[0] + unlocked_size_of_group[1]; };
    locked_of_cell.assign(num_of_cells(), false);

    while (num_of_unlocked_cells()) {
        LOG(DEBUG) << "---------------------------------- move " 
            << num_of_cells() - num_of_unlocked_cells() + 1
            << "th cell ----------------------------------" << endl;

        if (!move_one_cell(cut_size, min_cut_size, picked_No_after_min_cut, unlocked_size_of_group)) {
            break;
        }
    }
    
    LOG(INFO) << "min_cut_size = " << min_cut_size << " after " 
        << picked_No_after_min_cut << " cells moved" << endl;
    LOG(DEBUG) << "order of picked cells: ";
    RUN(DEBUG) for (size_t i = 0; i < cell_of_picked_No.size(); i++) {
        LOG(DEBUG) << "c" << cell_of_picked_No[i] << " ";
    }
    LOG(DEBUG) << endl;
    // trace back to the state with min cut size
    for (size_t i = cell_of_picked_No.size(); i --> picked_No_after_min_cut;) {
        cell_t cell = cell_of_picked_No[i];
        bool old_group = group_of_cell[cell], new_group = !old_group;
        group_of_cell[cell] = new_group;
        unlocked_size_of_group[new_group]++;
        size_of_group[old_group]--;
        size_of_group[new_group]++;
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
        return true;
    } else {
        LOG(INFO) << "cut size does not decrease, ending the loop." << endl;
        return false;
    }
}

bool FiducciaMattheysesPartitioner::move_one_cell(
    size_t& cut_size, size_t& min_cut_size, size_t& picked_No_after_min_cut,
    array<size_t, 2>& unlocked_size_of_group
) {
    bool permit_group_1_to_0 = size_of_group[0] + 1 <= (1 + balance_degree) / 2 * num_of_cells();
    bool permit_group_0_to_1 = size_of_group[0] - 1 >= (1 - balance_degree) / 2 * num_of_cells();

    // pick a cell to move
    cell_t cell_to_move;
    if (!permit_group_1_to_0 && !permit_group_0_to_1) {
        // WIP
        LOG(FETAL) << "Error: balance_degree is too large, cannot satisfy the balance constraint." << endl;
        exit(1);
    } else if (!permit_group_0_to_1 || unlocked_size_of_group[0] == 0) { // pick a cell in group 1
        if (unlocked_size_of_group[1] == 0) { // no cell to pick
            LOG(INFO) << "Fail to find a cell in group 1, ending the loop." << endl;
            return false;
        }
        LOG(DEBUG) << "group 0 too small, pick a cell in group 1" << endl;
        while (head_cell_of_group_by_gain(1, max_gain_of_group[1]) == -1) max_gain_of_group[1]--;
        cell_to_move = head_cell_of_group_by_gain(1, max_gain_of_group[1]);
    } else if (!permit_group_1_to_0 || unlocked_size_of_group[1] == 0) { // pick a cell in group 0
        if (unlocked_size_of_group[0] == 0) { // no cell to pick
            LOG(INFO) << "Fail to find a cell in group 0, ending the loop." << endl;
            return false;
        }
        LOG(DEBUG) << "group 1 too small, pick a cell in group 0" << endl;
        while (head_cell_of_group_by_gain(0, max_gain_of_group[0]) == -1) max_gain_of_group[0]--;
        cell_to_move = head_cell_of_group_by_gain(0, max_gain_of_group[0]);
    } else { // pick a cell in either group
        while (head_cell_of_group_by_gain(0, max_gain_of_group[0]) == -1) max_gain_of_group[0]--;
        while (head_cell_of_group_by_gain(1, max_gain_of_group[1]) == -1) max_gain_of_group[1]--;
        cell_t cell0 = head_cell_of_group_by_gain(0, max_gain_of_group[0]);
        cell_t cell1 = head_cell_of_group_by_gain(1, max_gain_of_group[1]);
        if (max_gain_of_group[0] == max_gain_of_group[1]) {
            cell_to_move = (gain_update_ts_of_cell[cell0] >= gain_update_ts_of_cell[cell1]) ? cell0 : cell1;
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

        bool FROM_group = group_of_cell[cell_to_move], TO_group = !FROM_group;
        size_t TO_group_cell_count = cell_count_of_net_by_group[net][TO_group];
        size_t FROM_group_cell_count = cell_count_of_net_by_group[net][FROM_group];

        if (TO_group_cell_count == 0) {
            for (cell_t cell : cells_of_net[net]) {
                if (!locked_of_cell[cell] && cell != cell_to_move) increase_gain_of_cell(cell, 1);
            }
        } else if (TO_group_cell_count == 1) {
            for (cell_t cell : cells_of_net[net]) {
                if (!locked_of_cell[cell] && cell != cell_to_move && group_of_cell[cell] == TO_group) {
                    increase_gain_of_cell(cell, -1);
                    break;
                }
            }
        }

        TO_group_cell_count++;
        FROM_group_cell_count--;

        if (FROM_group_cell_count == 0) {
            for (cell_t cell : cells_of_net[net]) {
                if (!locked_of_cell[cell] && cell != cell_to_move) increase_gain_of_cell(cell, -1);
            }
        } else if (FROM_group_cell_count == 1) {
            for (cell_t cell : cells_of_net[net]) {
                if (!locked_of_cell[cell] && cell != cell_to_move && group_of_cell[cell] == FROM_group) {
                    increase_gain_of_cell(cell, 1);
                    break;
                }
            }
        }
    }

    // move cell to another group
    bool old_group = group_of_cell[cell_to_move], new_group = !old_group;
    gain_t gain = gain_of_cell[cell_to_move];
    unlocked_size_of_group[old_group]--;
    size_of_group[old_group]--;
    erase_cell_of_gain(cell_to_move, gain, old_group);
    group_of_cell[cell_to_move] = new_group;
    size_of_group[new_group]++;
    for (auto net : nets_of_cell[cell_to_move]) {
        cell_count_of_net_by_group[net][old_group]--;
        cell_count_of_net_by_group[net][new_group]++;
    }
    cut_size -= gain_of_cell[cell_to_move];
    cell_of_picked_No.push_back(cell_to_move);
    if (cut_size < min_cut_size) {
        min_cut_size = cut_size;
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
                    cell_count--;
                }
                LOG(TRACE) << endl;
            }
        }
    }
    
    return true;
}

size_t FiducciaMattheysesPartitioner::get_cut_size() const {
    size_t cut_size = 0;
    for (net_t net = 0; net < cells_of_net.size(); net++) {
        if (cell_count_of_net_by_group[net][0] != 0 && cell_count_of_net_by_group[net][1] != 0) 
            cut_size++;
    }
    return cut_size;
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
    std::unordered_map<std::string, cell_t> index_of_cell_name;
    
    std::string tmp;
    while (input >> tmp) {
        input >> tmp; // net name
        std::set<cell_t> cells_in_this_net;
        for (;;) {
            input >> tmp;
            if (tmp == ";") break;

            cell_t cell;
            if (index_of_cell_name.find(tmp) == index_of_cell_name.end()) {
                cell = name_of_cell.size();
                name_of_cell.push_back(tmp);
                index_of_cell_name[tmp] = cell;
            } else {
                cell = index_of_cell_name[tmp];
            }
            cells_in_this_net.insert(cell);
        }
        partition_input.cells_of_net.emplace_back(BEGIN_END(cells_in_this_net));
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

PartitionOutput partition(
    const PartitionInput& input, 
    uint32_t repeat,
    double time_limit_seconds
) {
    FiducciaMattheysesPartitioner partitioner(input, repeat, time_limit_seconds);
    partitioner.solve();
    return partitioner.output;
}