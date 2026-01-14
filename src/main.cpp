#include "fft.h"
#include <iostream>

int main() {
    std::vector<float> realData = {1.0f, 2.0f, 3.0f, 2.0f};
    fft::CArray result = fft::fft(realData);

    for (const auto& val : result) {
        std::cout << std::abs(val) << std::endl;
    }
    return 0;
}