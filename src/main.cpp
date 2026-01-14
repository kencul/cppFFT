#include "fft.h"
#include <iostream>
#include <cmath>

#define M_PI 3.14159265358979323846f

using namespace std;

int main() {
    // Parameters
    const int sampleRate = 44100;
    const float frequency = 50.0f;
    const float dur = 5.f;
    const int bufferSize = 1024;

    size_t numSamples = static_cast<size_t>(sampleRate * dur);

    // genearate sine wave
    vector<float> buffer(bufferSize);

    for (size_t i = 0; i < bufferSize; ++i) {
        buffer[i] = sin(2.0f * M_PI * frequency * (static_cast<float>(i) / static_cast<float>(sampleRate)));
    }

    // vector<float> realData = {1.0f, 2.0f, 3.0f, 2.0f};
    fft::CArray result = fft::fft(buffer);
    vector<double> mags = fft::magnitude(result);

    for (const auto& mag : mags) {
        std::cout << mag << std::endl;
    }
    return 0;
}