#include "fft.h"
#include <iostream>
#include <cmath>

#define M_PI 3.14159265358979323846f

using namespace std;

int main() {
    // Parameters
    const int sampleRate = 44100;
    const float frequency =  67.f;
    const float dur = 10.f;
    const int bufferSize = pow(2, 15); // must be power of 2

    size_t numSamples = static_cast<size_t>(sampleRate * dur);

    // genearate sine wave
    vector<float> buffer(bufferSize);

    for (size_t i = 0; i < bufferSize; ++i) {
        buffer[i] = sin(2.0f * M_PI * frequency * (static_cast<float>(i) / static_cast<float>(sampleRate)));
    }

    buffer = fft::applyHanningWindow(buffer);
    fft::CArray result = fft::fft(buffer);
    vector<double> power = fft::power(result);
    auto freqScale = fft::frequencyScale(bufferSize, sampleRate);

    auto maxDb = max_element(power.begin(), power.end());
    int index = maxDb - power.begin();
    double freqAtMaxDb = freqScale[index];
    if(maxDb != power.end()) {
        std::cout << format("Max bin at {} Hz at {} dB", freqAtMaxDb, *maxDb) << std::endl;
    }

    return 0;
}