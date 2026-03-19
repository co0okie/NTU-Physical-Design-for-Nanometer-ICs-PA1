#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <cmath>
#include <list>
#include <algorithm>
#include <cassert>

#define BEGIN_END(container) (container).begin(), (container).end()

using std::vector, std::string, std::cout, std::endl, std::ifstream, std::ofstream, std::ostream, 
    std::unordered_map, std::list, std::max, std::min;

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

typedef uint64_t cell_t; // cell index type
typedef uint64_t net_t; // net index type
typedef int64_t gain_t; // gain type

struct Solution {
    size_t cut_size;
    vector<cell_t> group0_cells;
    vector<cell_t> group1_cells;
};

double balance_degree;

size_t num_of_nets = 0;
size_t num_of_cells = 0;
vector<string> name_of_net; // net index -> net name
vector<string> name_of_cell; // cell index -> cell name
vector<vector<cell_t>> cells_of_net; // net index -> cell indices
vector<vector<net_t>> nets_of_cell; // cell index -> net indices
unordered_map<string, cell_t> index_of_cell_name; // cell name -> cell index
unordered_map<string, net_t> index_of_net_name; // net name -> net index

size_t get_cut_size(const vector<vector<cell_t>> cells_of_net, const vector<net_t> cell_count_in_group0_of_net) {
    size_t cut_size = 0;
    for (net_t net = 0; net < cells_of_net.size(); net++) {
        if (cell_count_in_group0_of_net[net] != 0 
        && cell_count_in_group0_of_net[net] != cells_of_net[net].size()) 
            cut_size++;
    }
    return cut_size;
}

Solution fiduccia_mattheyses();

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


    input_file >> balance_degree;
    string tmp;
    while (input_file >> tmp) {
        input_file >> tmp;
        const net_t net_index = num_of_nets++;
        index_of_net_name[tmp] = net_index;
        name_of_net.push_back(tmp);
        cells_of_net.push_back({});
        for (;;) {
            input_file >> tmp;
            if (tmp == ";") break;

            cell_t cell_index;
            if (index_of_cell_name.find(tmp) == index_of_cell_name.end()) {
                cell_index = num_of_cells++;
                name_of_cell.push_back(tmp);
                index_of_cell_name[tmp] = cell_index;
                nets_of_cell.push_back({});
            } else {
                cell_index = index_of_cell_name[tmp];
            }
            cells_of_net[net_index].push_back(cell_index);
            nets_of_cell[cell_index].push_back(net_index);
        }
    }

    LOG(INFO) << "balance_degree: " << balance_degree << endl;
    LOG(INFO) << "\nnum_of_nets: " << num_of_nets << endl;
    for (net_t net = 0; net < num_of_nets; ++net) {
        LOG(DEBUG) << index_of_net_name[name_of_net[net]] << ":" << name_of_net[net] << "{";
        for (size_t j = 0; j < cells_of_net[net].size(); ++j) {
            if (j) { LOG(DEBUG) << " "; }
            LOG(DEBUG) << name_of_cell[cells_of_net[net][j]];
        }
        LOG(DEBUG) << "} ";
    }
    LOG(DEBUG) << endl;
    LOG(INFO) << "\nnum_of_cells: " << num_of_cells << endl;
    for (cell_t cell = 0; cell < num_of_cells; ++cell) {
        LOG(DEBUG) << index_of_cell_name[name_of_cell[cell]] << ":" << name_of_cell[cell] << "{";
        for (size_t i = 0; i < nets_of_cell[cell].size(); ++i) {
            if (i) { LOG(DEBUG) << " "; }
            LOG(DEBUG) << name_of_net[nets_of_cell[cell][i]];
        }
        LOG(DEBUG) << "} ";
    }
    LOG(DEBUG) << endl;

    if ((1 - balance_degree) / 2 * num_of_cells > num_of_cells / 2) {
        LOG(FETAL) << "Error: balance_degree is too small, cannot satisfy the balance constraint." << endl;
        return 1;
    }

    auto solution = fiduccia_mattheyses();
    
    output_file << "Cutsize = " << solution.cut_size << endl;
    output_file << "G1 " << solution.group0_cells.size() << endl;
    for (size_t i = 0; i < solution.group0_cells.size(); i++) {
        output_file << name_of_cell[solution.group0_cells[i]] << " ";
    }
    output_file << ";" << endl;
    output_file << "G2 " << solution.group1_cells.size() << endl;
    for (size_t i = 0; i < solution.group1_cells.size(); i++) {
        output_file << name_of_cell[solution.group1_cells[i]] << " ";
    }
    output_file << ";" << endl;

    return 0;    
}

Solution fiduccia_mattheyses() {
    vector<bool> group_of_cell(num_of_cells); // cell index -> group index (0 or 1)
    vector<size_t> cell_count_in_group0_of_net(num_of_nets, 0); // net index -> number of cells in group 0 connected to this net
    auto cell_count_in_group1_of_net = [&](net_t net) {
        return cells_of_net[net].size() - cell_count_in_group0_of_net[net]; 
    };
    size_t group0_size = num_of_cells / 2;
    auto group1_size = [&]() { return num_of_cells - group0_size; };
    for (cell_t cell = 0; cell < num_of_cells; cell++) {
        group_of_cell[cell] = cell * 2 < num_of_cells;
        if (group_of_cell[cell] == 0) {
            for (size_t j = 0; j < nets_of_cell[cell].size(); j++) {
                cell_count_in_group0_of_net[nets_of_cell[cell][j]]++;
            }
        }
    }
    
    LOG(INFO) << "\ngroup0_size " << group0_size << " group1_size " << group1_size() << endl;
    LOG(DEBUG) << "\n<cell_name>:<group_index> " << endl;
    for (cell_t cell = 0; cell < num_of_cells; cell++) {
        LOG(DEBUG) << name_of_cell[cell] << ":" << group_of_cell[cell] << " ";
    }
    LOG(DEBUG) << endl;
    LOG(DEBUG) << "\n<net_name>:<group 0 cell count>/<group 1 cell count> " << endl;
    for (net_t net = 0; net < num_of_nets; net++) {
        LOG(DEBUG) << name_of_net[net] << ":" << cell_count_in_group0_of_net[net] 
            << "/" << cell_count_in_group1_of_net(net) << " ";
    }
    LOG(DEBUG) << endl;

    size_t iteration_min_cut_size = get_cut_size(cells_of_net, cell_count_in_group0_of_net);
    LOG(INFO) << "\ninitial cut size: " << iteration_min_cut_size << endl;
    for (size_t iter = 0;; iter++) {
        LOG(INFO) << "\n================================ iteration " << iter + 1 
            << " ================================" << endl;

        gain_t group0_max_gain = INT64_MIN;
        gain_t group1_max_gain = INT64_MIN;
        const gain_t MAX_GAIN = num_of_nets; // the maximum gain of a cell is when all its nets are cut nets and it moves to the other group, then all its nets become uncut nets, gain = number of nets connected to this cell <= num_of_nets
        vector<list<cell_t>> group0_cells_of_offset_gain(2 * MAX_GAIN + 1); // gain -> list of cell indices in group 0 with this gain, offset by num_of_nets to allow negative gain
        vector<list<cell_t>> group1_cells_of_offset_gain(2 * MAX_GAIN + 1); // gain -> list of cell indices in group 1 with this gain, offset by num_of_nets to allow negative gain
        auto group0_cells_of_gain = [&](gain_t gain) -> list<cell_t>& { return group0_cells_of_offset_gain[gain + MAX_GAIN]; };
        auto group1_cells_of_gain = [&](gain_t gain) -> list<cell_t>& { return group1_cells_of_offset_gain[gain + MAX_GAIN]; };
        vector<list<cell_t>::iterator> iterator_of_cell(num_of_cells); // cell index -> iterator to the cell's position in cells_of_gain
        vector<gain_t> gain_of_cell(num_of_cells); // cell index -> gain of this cell
        for (cell_t cell = 0; cell < num_of_cells; cell++) {
            gain_t gain = 0;
            for (size_t j = 0; j < nets_of_cell[cell].size(); j++) {
                const net_t net = nets_of_cell[cell][j];
                if (group_of_cell[cell] == 0) { // cell i at group 0, try to move it to group 1
                    if (cell_count_in_group0_of_net[net] == 1) gain++;
                    else if (cell_count_in_group0_of_net[net] == cells_of_net[net].size()) gain--;
                } else { // cell i at group 1, try to move it to group 0
                    if (cell_count_in_group0_of_net[net] == 0) gain--;
                    else if (cell_count_in_group0_of_net[net] == cells_of_net[net].size() - 1) gain++;
                }
            }
            if (group_of_cell[cell] == 0) {
                group0_cells_of_gain(gain).push_front(cell);
                iterator_of_cell[cell] = group0_cells_of_gain(gain).begin();
                group0_max_gain = max(group0_max_gain, gain);
            } else {
                group1_cells_of_gain(gain).push_front(cell);
                iterator_of_cell[cell] = group1_cells_of_gain(gain).begin();
                group1_max_gain = max(group1_max_gain, gain);
            }
            gain_of_cell[cell] = gain;
        }
        
        RUN(DEBUG) {
            LOG(DEBUG) << "group0 gain table: " << endl;
            size_t cell_count = group0_size;
            for (gain_t gain = group0_max_gain; cell_count; gain--) {
                LOG(DEBUG) << "gain " << gain << ": ";
                for (auto cell : group0_cells_of_gain(gain)) { 
                    LOG(DEBUG) << name_of_cell[cell] << " ";
                    cell_count--;
                }
                LOG(DEBUG) << endl;
            }
            LOG(DEBUG) << "group1 gain table: " << endl;
            cell_count = group1_size();
            for (gain_t gain = group1_max_gain; cell_count; gain--) {
                LOG(DEBUG) << "gain " << gain << ": ";
                for (auto cell : group1_cells_of_gain(gain)) { 
                    LOG(DEBUG) << name_of_cell[cell] << " ";
                    cell_count--;
                }
                LOG(DEBUG) << endl;
            }
        }
        
        // pick the cell with the maximum gain and move it to the other group, then update the gain of the affected cells, repeat until every cell has been moved once
        size_t cut_size = iteration_min_cut_size;
        size_t min_cut_size = iteration_min_cut_size;
        vector<cell_t> cell_of_picked_No; // the cell index of the cell picked at each step, in order
        cell_of_picked_No.reserve(num_of_cells);
        size_t picked_No_after_min_cut = 0; // {cut 5} {pick cell 1} {cut 2} {pick cell 0} {cut 6}, in this case picked_No_after_min_cut = 1
        size_t unlocked_group0_size = group0_size; // the number of unlocked cells in group 0
        vector<bool> locked_of_cell(num_of_cells, false); // cell index -> {0: unlocked, 1: locked}
        for (size_t num_of_unlocked_cells = num_of_cells; num_of_unlocked_cells;) {
            LOG(DEBUG) << "---------------------------------- move " << num_of_cells - num_of_unlocked_cells + 1
                << "th cell ----------------------------------" << endl;
            auto unlocked_group1_size = [&]() { return num_of_unlocked_cells - unlocked_group0_size; };
        
            const bool permit_group_1_to_0 = group0_size + 1 <= (1 + balance_degree) / 2 * num_of_cells;
            const bool permit_group_0_to_1 = group0_size - 1 >= (1 - balance_degree) / 2 * num_of_cells;
        
            // pick a cell to move
            cell_t cell_to_move;
            if (!permit_group_1_to_0 && !permit_group_0_to_1) {
                // WIP
                LOG(FETAL) << "Error: balance_degree is too large, cannot satisfy the balance constraint." << endl;
                exit(1);
            } else if (!permit_group_0_to_1 || unlocked_group0_size == 0) { // pick a cell in group 1
                if (unlocked_group1_size() == 0) { // no cell to pick
                    LOG(INFO) << "Fail to find a cell in group 1, ending the loop." << endl;
                    break;
                }
                LOG(DEBUG) << "group 0 too small, pick a cell in group 1" << endl;
                while (group1_cells_of_gain(group1_max_gain).empty()) group1_max_gain--;
                cell_to_move = group1_cells_of_gain(group1_max_gain).front();
            } else if (!permit_group_1_to_0 || unlocked_group1_size() == 0) { // pick a cell in group 0
                if (unlocked_group0_size == 0) { // no cell to pick
                    LOG(INFO) << "Fail to find a cell in group 0, ending the loop." << endl;
                    break;
                }
                LOG(DEBUG) << "group 1 too small, pick a cell in group 0" << endl;
                while (group0_cells_of_gain(group0_max_gain).empty()) group0_max_gain--;
                cell_to_move = group0_cells_of_gain(group0_max_gain).front();
            } else { // pick a cell in either group
                while (group0_cells_of_gain(group0_max_gain).empty()) group0_max_gain--;
                while (group1_cells_of_gain(group1_max_gain).empty()) group1_max_gain--;
                cell_to_move = (
                    (group0_max_gain > group1_max_gain) ? 
                    group0_cells_of_gain(group0_max_gain) :
                    group1_cells_of_gain(group1_max_gain)
                ).front();
            }
        
            LOG(DEBUG) << "group0_max_gain " << group0_max_gain << " group1_max_gain " << group1_max_gain << endl;
            LOG(DEBUG) << "move " << name_of_cell[cell_to_move] << " from " 
                << group_of_cell[cell_to_move] << " to " << !group_of_cell[cell_to_move]
                << " gain " << gain_of_cell[cell_to_move] << " cut_size " << cut_size
                << " -> " << cut_size - gain_of_cell[cell_to_move] << endl;
        
            // update the gain of the affected cells
            for (size_t i = 0; i < nets_of_cell[cell_to_move].size(); i++) {
                net_t net = nets_of_cell[cell_to_move][i];
                bool from_1_to_0 = group_of_cell[cell_to_move], from_0_to_1 = !from_1_to_0;
                bool become_pair = (from_0_to_1 && cell_count_in_group1_of_net(net) == 1)
                                || (from_1_to_0 && cell_count_in_group0_of_net[net] == 1);
                bool escape_pair = (from_0_to_1 && cell_count_in_group0_of_net[net] == 2)
                                || (from_1_to_0 && cell_count_in_group1_of_net(net) == 2);
                bool become_empty = (from_0_to_1 && cell_count_in_group0_of_net[net] == 1)
                                || (from_1_to_0 && cell_count_in_group1_of_net(net) == 1);
                bool escape_empty = (from_0_to_1 && cell_count_in_group1_of_net(net) == 0)
                                || (from_1_to_0 && cell_count_in_group0_of_net[net] == 0);
                
                // no gain change
                if (!escape_pair && !become_pair && !become_empty && !escape_empty) continue;
        
                LOG(TRACE) << "[" << name_of_net[net] << " "
                    << cell_count_in_group0_of_net[net] << "/" 
                    << cell_count_in_group1_of_net(net) << " "
                    << (become_pair ? "become_pair " : "") 
                    << (escape_pair ? "escape_pair " : "") 
                    << (become_empty ? "become_empty " : "") 
                    << (escape_empty ? "escape_empty " : "");
        
                for (size_t j = 0; j < cells_of_net[net].size(); j++) {
                    cell_t cell = cells_of_net[net][j];
                    if (cell == cell_to_move || locked_of_cell[cell]) continue;
        
                    int gain_change = 0;
                    if (become_pair && ((from_0_to_1 && group_of_cell[cell] == 1) 
                    || (from_1_to_0 && group_of_cell[cell] == 0))) {
                        gain_change -= 1;
                    } else if (escape_empty) {
                        gain_change += 1;
                    }
                    if (escape_pair && ((from_0_to_1 && group_of_cell[cell] == 0) 
                    || (from_1_to_0 && group_of_cell[cell] == 1))) {
                        gain_change += 1;
                    } else if (become_empty) {
                        gain_change -= 1;
                    }
        
                    if (gain_change == 0) continue;
        
                    LOG(TRACE) << "(" << name_of_cell[cell] 
                        << " gain " << gain_of_cell[cell] << " -> " 
                        << gain_of_cell[cell] + gain_change << ")" << std::flush;
                    if (group_of_cell[cell] == 0) {
                        group0_cells_of_gain(gain_of_cell[cell]).erase(iterator_of_cell[cell]);
                        gain_of_cell[cell] += gain_change;
                        group0_cells_of_gain(gain_of_cell[cell]).push_front(cell);
                        iterator_of_cell[cell] = group0_cells_of_gain(gain_of_cell[cell]).begin();
                        group0_max_gain = max(group0_max_gain, gain_of_cell[cell]);
                    } else {
                        group1_cells_of_gain(gain_of_cell[cell]).erase(iterator_of_cell[cell]);
                        gain_of_cell[cell] += gain_change;
                        group1_cells_of_gain(gain_of_cell[cell]).push_front(cell);
                        iterator_of_cell[cell] = group1_cells_of_gain(gain_of_cell[cell]).begin();
                        group1_max_gain = max(group1_max_gain, gain_of_cell[cell]);
                    }
                }
        
                LOG(TRACE) << "]" << std::flush;
            }
            LOG(TRACE) << endl;
        
            // move cell to another group
            if (group_of_cell[cell_to_move] == 0) {
                group_of_cell[cell_to_move] = 1;
                group0_size--;
                unlocked_group0_size--;
                for (auto net : nets_of_cell[cell_to_move]) {
                    cell_count_in_group0_of_net[net]--;
                }
            } else {
                group_of_cell[cell_to_move] = 0;
                group0_size++;
                for (auto& net : nets_of_cell[cell_to_move]) {
                    cell_count_in_group0_of_net[net]++;
                }
            }
            num_of_unlocked_cells--;
            cut_size -= gain_of_cell[cell_to_move];
            cell_of_picked_No.push_back(cell_to_move);
            if (cut_size < min_cut_size) {
                min_cut_size = cut_size;
                picked_No_after_min_cut = cell_of_picked_No.size();
            }
            auto& cells = (group_of_cell[cell_to_move] == 0) ? 
                group0_cells_of_gain(gain_of_cell[cell_to_move]) : 
                group1_cells_of_gain(gain_of_cell[cell_to_move]);
            cells.erase(iterator_of_cell[cell_to_move]);
            locked_of_cell[cell_to_move] = true;
        
            RUN(TRACE) {
                LOG(TRACE) << "group0 gain table: " << endl;
                size_t cell_count = unlocked_group0_size;
                for (gain_t gain = group0_max_gain; cell_count; gain--) {
                    LOG(TRACE) << "gain " << gain << ": ";
                    for (auto cell : group0_cells_of_gain(gain)) { 
                        LOG(TRACE) << name_of_cell[cell] << " ";
                        cell_count--;
                    }
                    LOG(TRACE) << endl;
                }
                LOG(TRACE) << "group1 gain table: " << endl;
                cell_count = unlocked_group1_size();
                for (gain_t gain = group1_max_gain; cell_count; gain--) {
                    LOG(TRACE) << "gain " << gain << ": ";
                    for (auto cell : group1_cells_of_gain(gain)) { 
                        LOG(TRACE) << name_of_cell[cell] << " ";
                        cell_count--;
                    }
                    LOG(TRACE) << endl;
                }
            }
        }
        
        LOG(INFO) << "\nmin_cut_size = " << min_cut_size << " after " 
            << picked_No_after_min_cut << " cells moved" << endl;
        LOG(DEBUG) << "order of picked cells: ";
        RUN(DEBUG) for (size_t i = 0; i < num_of_cells; i++) {
            LOG(DEBUG) << name_of_cell[cell_of_picked_No[i]] << " ";
        }
        LOG(DEBUG) << endl;
        // trace back to the state with min cut size
        for (size_t i = picked_No_after_min_cut; i < num_of_cells; i++) {
            cell_t cell = cell_of_picked_No[i];
            group_of_cell[cell] = !group_of_cell[cell];
            if (group_of_cell[cell] == 0) {
                group0_size++;
                unlocked_group0_size++;
                for (auto net : nets_of_cell[cell]) {
                    cell_count_in_group0_of_net[net]++;
                }
            } else {
                group0_size--;
                for (auto& net : nets_of_cell[cell]) {
                    cell_count_in_group0_of_net[net]--;
                }
            }
        }
        cut_size = min_cut_size;

        LOG(INFO) << "after trace back: group0_size " << group0_size << " group1_size " 
            << group1_size() << endl;
        for (cell_t cell = 0; cell < num_of_cells; cell++) {
            LOG(DEBUG) << name_of_cell[cell] << ":" << group_of_cell[cell] << " ";
        }
        LOG(DEBUG) << endl;

        if (min_cut_size < iteration_min_cut_size) {
            iteration_min_cut_size = min_cut_size;
        } else {
            LOG(INFO) << "cut size does not decrease, ending the loop." << endl;
            break;
        }
    }

    Solution solution;
    solution.cut_size = get_cut_size(cells_of_net, cell_count_in_group0_of_net);
    for (cell_t cell = 0; cell < num_of_cells; cell++) {
        if (group_of_cell[cell] == 0) solution.group0_cells.push_back(cell);
        else solution.group1_cells.push_back(cell);
    }

    return solution;
}
