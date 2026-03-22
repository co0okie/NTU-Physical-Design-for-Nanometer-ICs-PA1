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

#define BEGIN_END(container) (container).begin(), (container).end()

using std::vector, std::string, std::cout, std::endl, std::ifstream, std::ofstream, std::ostream, 
    std::unordered_map, std::list, std::max, std::min, std::array;

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

constexpr double time_limit_seconds = 60.0;

typedef int32_t cell_t; // cell index type
typedef int32_t net_t; // net index type
typedef int32_t gain_t; // gain type
constexpr gain_t MIN_GAIN = INT32_MIN;

#define RUN(level) if (level > LOG_LEVEL) {} else
#define LOG(level) RUN(level) std::cout

struct Solution {
    size_t cut_size;
    array<vector<cell_t>, 2> cells_of_group;
};

double balance_degree;

size_t num_of_nets = 0;
size_t num_of_cells = 0;

vector<string> name_of_cell; // cell -> cell name
vector<vector<net_t>> nets_of_cell; // cell -> net indices
vector<bool> group_of_cell; // cell -> group index (0 or 1)
vector<gain_t> gain_of_cell; // cell -> gain of this cell
vector<uint32_t> gain_update_ts_of_cell; // cell -> timestamp of the last gain update of this cell, to choose the latest updated cell
vector<bool> locked_of_cell; // cell -> {0: unlocked, 1: locked}

vector<string> name_of_net; // net -> net name
vector<vector<cell_t>> cells_of_net; // net -> cell indices
vector<array<size_t, 2>> cell_count_of_net_by_group; // net -> group -> number of cells in this group connected to this net

unordered_map<string, cell_t> index_of_cell_name; // cell name -> cell index
unordered_map<string, net_t> index_of_net_name; // net name -> net index
vector<cell_t> cell_of_picked_No; // cell picked order sequence -> cell index

array<vector<cell_t>, 2> head_cell_of_group_by_offset_gain; // group -> gain -> gain list head cell with this gain in this group
vector<cell_t> prev_of_cell; // cell -> previous cell in the same gain list
vector<cell_t> next_of_cell; // cell -> next cell in the same gain list

uint32_t gain_update_ts;

gain_t gain_offset = 0;

array<size_t, 2> size_of_group; // group -> number of cells in this group
array<gain_t, 2> max_gain_of_group;

size_t get_cut_size() {
    size_t cut_size = 0;
    for (net_t net = 0; net < cells_of_net.size(); net++) {
        if (cell_count_of_net_by_group[net][0] != 0 && cell_count_of_net_by_group[net][1] != 0) 
            cut_size++;
    }
    return cut_size;
}

Solution fiduccia_mattheyses();
void parse_input(ifstream& input_file);
void write_output(ofstream& output_file, const Solution& solution);
void increase_gain_of_cell(cell_t cell, gain_t gain_change);
bool move_one_cell(size_t& cut_size, size_t& min_cut_size, size_t& picked_No_after_min_cut, array<size_t, 2>& unlocked_size_of_group);
bool do_one_iteration(size_t& iteration_min_cut_size, size_t iter);

cell_t& head_cell_of_group_by_gain(bool group, gain_t gain) {
    return head_cell_of_group_by_offset_gain[group][gain + gain_offset];
}

void push_cell_to_front_of_gain(cell_t cell, gain_t gain, bool group) {
    cell_t head_cell = head_cell_of_group_by_gain(group, gain);
    prev_of_cell[cell] = -1;
    next_of_cell[cell] = head_cell;
    if (head_cell != -1) prev_of_cell[head_cell] = cell;
    head_cell_of_group_by_gain(group, gain) = cell;
    gain_update_ts_of_cell[cell] = gain_update_ts;
}

void erase_cell_of_gain(cell_t cell, gain_t gain, bool group) {
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

void increase_gain_of_cell(cell_t cell, gain_t gain_change) {
    if (gain_change == 0) return;
    bool group = group_of_cell[cell];
    gain_t old_gain = gain_of_cell[cell], new_gain = old_gain + gain_change;
    erase_cell_of_gain(cell, old_gain, group);
    gain_of_cell[cell] = new_gain;
    push_cell_to_front_of_gain(cell, new_gain, group);
    max_gain_of_group[group] = max(max_gain_of_group[group], new_gain);
    gain_update_ts_of_cell[cell] = gain_update_ts;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cout << "Usage: " << argv[0] << " <input file name> <output file name>" << endl;
        return 1;
    }

    ifstream input_file(argv[1]);
    if (!input_file.is_open()) {
        LOG(FETAL) << "Error: Could not open file " << argv[1] << endl;
        return 1;
    }
    ofstream output_file(argv[2]);
    if (!output_file.is_open()) {
        LOG(FETAL) << "Error: Could not open file " << argv[2] << endl;
        return 1;
    }

    parse_input(input_file);

    LOG(INFO) << "balance_degree: " << balance_degree << endl;
    LOG(INFO) << "\nnum_of_nets: " << num_of_nets << endl;
    for (net_t net = 0; net < num_of_nets; ++net) {
        LOG(DEBUG) << index_of_net_name[name_of_net[net]] << ":" << name_of_net[net] << "{";
        for (size_t j = 0; j < cells_of_net[net].size(); ++j) {
            if (j) { LOG(DEBUG) << " "; }
            LOG(DEBUG) << name_of_cell[cells_of_net[net][j]];
        }
        LOG(DEBUG) << "}\n";
    }
    LOG(DEBUG) << endl;
    LOG(INFO) << "\nnum_of_cells: " << num_of_cells << endl;
    for (cell_t cell = 0; cell < num_of_cells; ++cell) {
        LOG(DEBUG) << index_of_cell_name[name_of_cell[cell]] << ":" << name_of_cell[cell] << "{";
        for (size_t i = 0; i < nets_of_cell[cell].size(); ++i) {
            if (i) { LOG(DEBUG) << " "; }
            LOG(DEBUG) << name_of_net[nets_of_cell[cell][i]];
        }
        LOG(DEBUG) << "}\n";
    }
    LOG(DEBUG) << endl;

    if ((1 - balance_degree) / 2 * num_of_cells > num_of_cells / 2) {
        LOG(FETAL) << "Error: balance_degree is too small, cannot satisfy the balance constraint." << endl;
        return 1;
    }

    Solution best_global_solution;
    best_global_solution.cut_size = SIZE_MAX;

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < 10; i++) {
        auto solution = fiduccia_mattheyses();

        if (solution.cut_size < best_global_solution.cut_size) {
            best_global_solution = solution;
        }

        auto current_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = current_time - start_time;
        if (elapsed.count() > time_limit_seconds) break;
    }
    
    write_output(output_file, best_global_solution);

    return 0;    
}

void parse_input(ifstream& input_file) {
    input_file >> balance_degree;
    string tmp;
    while (input_file >> tmp) {
        input_file >> tmp;
        const net_t net = num_of_nets++;
        index_of_net_name[tmp] = net;
        name_of_net.push_back(tmp);
        std::set<cell_t> cells_in_this_net; // to avoid duplicate cells in a net
        for (;;) {
            input_file >> tmp;
            if (tmp == ";") break;

            cell_t cell;
            if (index_of_cell_name.find(tmp) == index_of_cell_name.end()) {
                cell = num_of_cells++;
                name_of_cell.push_back(tmp);
                index_of_cell_name[tmp] = cell;
                nets_of_cell.push_back({});
            } else {
                cell = index_of_cell_name[tmp];
            }
            auto result = cells_in_this_net.insert(cell);
            if (result.second) {
                // new cell in this net
                nets_of_cell[cell].push_back(net);
            }
        }
        cells_of_net.emplace_back(BEGIN_END(cells_in_this_net));
    }

    gain_offset = std::max_element(BEGIN_END(nets_of_cell), [](const vector<net_t>& a, const vector<net_t>& b) {
        return a.size() < b.size();
    })->size();
}

void write_output(ofstream& output_file, const Solution& solution) {
    output_file << "Cutsize = " << solution.cut_size << endl;
    for (int group = 0; group < 2; group++) {
        output_file << "G" << group + 1 << " " << solution.cells_of_group[group].size() << endl;
        for (cell_t cell : solution.cells_of_group[group]) {
            output_file << name_of_cell[cell] << " ";
        }
        output_file << ";" << endl;
    }
}

Solution fiduccia_mattheyses() {
    group_of_cell.resize(num_of_cells);
    cell_count_of_net_by_group.assign(num_of_nets, {0, 0});
    size_of_group = {0, 0};

    std::vector<size_t> random_cells(num_of_cells);
    std::iota(BEGIN_END(random_cells), 0);
    std::random_shuffle(BEGIN_END(random_cells));
    
    for (cell_t cell = 0; cell < num_of_cells; cell++) {
        bool group = random_cells[cell] < num_of_cells / 2;
        group_of_cell[cell] = group;
        size_of_group[group]++;
        for (net_t net : nets_of_cell[cell]) {
            cell_count_of_net_by_group[net][group]++;
        }
    }
    
    LOG(INFO) << "\ngroup0_size " << size_of_group[0] << " group1_size " << size_of_group[1] << endl;
    LOG(DEBUG) << "\n<cell_name>:<group_index> " << endl;
    for (cell_t cell = 0; cell < num_of_cells; cell++) {
        LOG(DEBUG) << name_of_cell[cell] << ":" << group_of_cell[cell] << " ";
    }
    LOG(DEBUG) << endl;
    LOG(DEBUG) << "\n<net_name>:<group 0 cell count>/<group 1 cell count> " << endl;
    for (net_t net = 0; net < num_of_nets; net++) {
        LOG(DEBUG) << name_of_net[net] << ":" << cell_count_of_net_by_group[net][0] 
            << "/" << cell_count_of_net_by_group[net][1] << " ";
    }
    LOG(DEBUG) << endl;

    size_t iteration_min_cut_size = get_cut_size();
    LOG(INFO) << "\ninitial cut size: " << iteration_min_cut_size << endl;
    for (size_t iter = 0;; iter++) {
        if (!do_one_iteration(iteration_min_cut_size, iter)) {
            break;
        }
    }

    Solution solution;
    solution.cut_size = get_cut_size();
    if (solution.cut_size != iteration_min_cut_size) {
        LOG(FETAL) << "Error: cut size " << iteration_min_cut_size 
            << " after trace back does not equal to the minimum cut size " << solution.cut_size 
            << " found in this iteration." << endl;
        exit(1);
    }
    for (cell_t cell = 0; cell < num_of_cells; cell++) {
        solution.cells_of_group[group_of_cell[cell]].push_back(cell);
    }

    return solution;
}

bool do_one_iteration(size_t& iteration_min_cut_size, size_t iter) {
    LOG(INFO) << "\n================================ iteration " << iter + 1 
        << " ================================" << endl;

    max_gain_of_group = {MIN_GAIN, MIN_GAIN};
    head_cell_of_group_by_offset_gain[0].assign(2 * gain_offset + 1, -1);
    head_cell_of_group_by_offset_gain[1].assign(2 * gain_offset + 1, -1);
    prev_of_cell.assign(num_of_cells, -1);
    next_of_cell.assign(num_of_cells, -1);
    gain_of_cell.resize(num_of_cells);
    gain_update_ts = 0;
    gain_update_ts_of_cell.assign(num_of_cells, gain_update_ts);
    
    for (cell_t cell = 0; cell < num_of_cells; cell++) {
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
                    LOG(DEBUG) << name_of_cell[cell] << " ";
                    cell_count--;
                }
                LOG(DEBUG) << endl;
            }
        }
    }
    
    size_t cut_size = iteration_min_cut_size;
    size_t min_cut_size = iteration_min_cut_size;
    cell_of_picked_No.clear();
    cell_of_picked_No.reserve(num_of_cells);
    size_t picked_No_after_min_cut = 0; // {cut 5} {pick cell 1} {cut 2} {pick cell 0} {cut 6}, in this case picked_No_after_min_cut = 1
    array<size_t, 2> unlocked_size_of_group = {size_of_group[0], size_of_group[1]}; // group -> number of unlocked cells in this group
    auto num_of_unlocked_cells = [&]() { return unlocked_size_of_group[0] + unlocked_size_of_group[1]; };
    locked_of_cell.assign(num_of_cells, false);

    while (num_of_unlocked_cells()) {
        LOG(DEBUG) << "---------------------------------- move " 
            << num_of_cells - num_of_unlocked_cells() + 1
            << "th cell ----------------------------------" << endl;

        if (!move_one_cell(cut_size, min_cut_size, picked_No_after_min_cut, unlocked_size_of_group)) {
            break;
        }
    }
    
    LOG(INFO) << "\nmin_cut_size = " << min_cut_size << " after " 
        << picked_No_after_min_cut << " cells moved" << endl;
    LOG(DEBUG) << "order of picked cells: ";
    RUN(DEBUG) for (size_t i = 0; i < cell_of_picked_No.size(); i++) {
        LOG(DEBUG) << name_of_cell[cell_of_picked_No[i]] << " ";
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
    for (cell_t cell = 0; cell < num_of_cells; cell++) {
        LOG(DEBUG) << name_of_cell[cell] << ":" << group_of_cell[cell] << " ";
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

bool move_one_cell(size_t& cut_size, size_t& min_cut_size, size_t& picked_No_after_min_cut, array<size_t, 2>& unlocked_size_of_group) {
    bool permit_group_1_to_0 = size_of_group[0] + 1 <= (1 + balance_degree) / 2 * num_of_cells;
    bool permit_group_0_to_1 = size_of_group[0] - 1 >= (1 - balance_degree) / 2 * num_of_cells;

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
    LOG(DEBUG) << "move " << name_of_cell[cell_to_move] << " from " 
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
                    LOG(TRACE) << name_of_cell[cell] << " ";
                    cell_count--;
                }
                LOG(TRACE) << endl;
            }
        }
    }
    
    return true;
}