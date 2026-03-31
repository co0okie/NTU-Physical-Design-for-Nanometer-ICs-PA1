#ifndef PARTITION_H
#define PARTITION_H

#include <cstdint>
#include <vector>
#include <unordered_set>
#include <array>
#include <string>
#include <tuple>
#include "constant.h"

struct PartitionInput {
    double balance_degree;
    std::vector<std::unordered_set<cell_t>> cells_of_net; // net -> cell indices

    static std::tuple<
        PartitionInput, 
        std::vector<std::string> // name_of_cell: cell -> cell name
    > from_file(const std::string& input_file);
};

struct PartitionOutput {
    uint64_t cut_size;
    std::array<std::vector<cell_t>, 2> cells_of_group; // group -> cell indices

    void to_file(const std::string& output_file, const std::vector<std::string>& name_of_cell) const;
};

PartitionOutput partition(
    const PartitionInput& input, 
    double time_limit_seconds = 60 * 60
);

#endif