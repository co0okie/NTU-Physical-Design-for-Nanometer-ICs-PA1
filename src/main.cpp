#include <iostream>
#include "partition.h"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cout << "Usage: " << argv[0] << " <input file name> <output file name>" << endl;
        return 1;
    }

    auto input = PartitionInput::from_file(argv[1]);
    PartitionOutput output = partition(std::get<0>(input), 16, 60 * 59);
    output.to_file(argv[2], std::get<1>(input));

    return 0;
}