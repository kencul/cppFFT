#ifndef FFT_H
#define FFT_H

#include <vector>
#include <complex>

namespace fft {

    // Type aliases
    using Complex = std::complex<double>;
    using CArray = std::vector<Complex>;

    /**
     * @brief Perform Fast Fourier Transform (FFT) on a complex-valued array.
     * @param x Reference to the input/output array of complex numbers.
     */
    void fft(CArray &buffer);

    /**
     * @brief Perform Fast Fourier Transform (FFT) on a real-valued vector.
     * @param x Input vector of real numbers (floats).
     * @return Vector of complex numbers representing the FFT result.
     */
    std::vector<Complex> fft(const std::vector<float>& buffer);

    std::vector<double> magnitude(const CArray& spectrum);

    std::vector<double> power(const CArray& spectrum);
}

#endif // FFT_H