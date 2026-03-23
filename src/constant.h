#ifndef CONSTANT_H
#define CONSTANT_H

#include <cstdint>

#define BEGIN_END(container) (container).begin(), (container).end()

typedef int32_t cell_t; // cell index type
typedef int32_t net_t; // net index type
typedef int32_t gain_t; // gain type
constexpr gain_t MIN_GAIN = INT32_MIN;

#endif