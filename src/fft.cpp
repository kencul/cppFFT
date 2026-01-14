#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

using namespace std;

namespace fft {
    // Define a complex number type. The complex<double> class holds a double for both the real and imaginary parts. 
    // Complex variables can be manipulated using standard arithmetic operators to perform complex number operations, so it looks simpler.
    using Complex = complex<double>; 
    // Defining a type alias CArray for a vector of Complex numbers to simplify code readability.
    using CArray = vector<Complex>;
    
    // Get PI value
    const double PI = acos(-1);
    
    /*
    Implementation of the Cooley-Tukey algorithm for Fast Fourier Transform (FFT) in C++.
    REF: https://www.w3computing.com/articles/how-to-implement-a-fast-fourier-transform-fft-in-cpp/
    */
    /// @brief Perform Fast Fourier Transform (FFT) on a complex-valued array.
    /// @param x Reference to the input/output array of complex numbers. Must be a length that is a power of 2.
    void fft(CArray &buffer) {
        const size_t N = buffer.size();
        if (N <= 1) return;
        
        // Ensure N is a power of 2
        // The FFT algorithm requires the input size to be a power of two to function
        if (N > 1 && (N & (N - 1)) != 0) {
            cerr << "Error: FFT requires input size to be a power of two. Got: " << N << endl;
            return;
        }

        // Bit-reversed addressing permutation
        /*
        Reversing the bits of the indices of the input array sorts the data in the order they need to be processed for the FFT algorithm.
        The bit representation of the indices reversed represents the grouping of data when the data is divided into even and odd groups.
        By reversing the bits of the indices, the data can be processed by iterating through the array, making the memory access patterns more efficient, thus improving the performance of the FFT algorithm.
        This for loop achieves this by making j the mirror image of i in binary form.
        */
        size_t j = 0;
        for (size_t i = 1; i < N; ++i) {
            /*
            As N (size of array) is always a power of 2 in FFT, the binary representation of N will have only one '1' bit followed by '0's.
            By bitshifting by 1, we get the most significant bit (MSB) position for the current size of the array.
            The bit variable is used as a mask to iterate through all the bits of the index, and flip them accordingly to achieve the bit-reversed order.
            */
            size_t bit = N >> 1;

            // Keep flipping 1's into 0's until a 0 is found
            // This works because the reversed bit is being carried over to the next bit position, essentially counting in binary but mirrored.
            while (j & bit) {
                // bitwise XOR operation
                // Reverse the current bit (set the 1 to 0)wh
                j ^= bit;
                // Iterate to next bit
                bit >>= 1;
            }
            // Flip the first 0 found to a 1
            j ^= bit;

            // j now represents the mirror image of i in binary form, so swap the values at indices i and j
            // Only do the swap if j is farther in the array than i to avoid swapping them back
            if (i < j) {
                swap(buffer[i], buffer[j]);
            }
        }

        // Iterative FFT
        /*
        This is the merging part of the FFT algorithm, where smaller DFTs are combined to form larger DFTs.
        This works because a DFT of a single sample is just the sample itself, so each sample is treated as a DFT of size 1.
        These samples are then combined in pairs to form DFTs of size 2, which are then combined to form DFTs of size 4, and so on, until the entire array is processed.
        The samples must be skewed by the appropriate twiddle factors (complex exponentials) to ensure correct phase relationships during the combination, as each sample was originally offset in time.
        */
        // Start combining from length 2 to N, doubling each time (len <<= 1 means bitshift left, or len = len * 2)
        for (size_t len = 2; len <= N; len <<= 1) {
            // Angle for a full circle divided by the current length
            double angle = -2 * PI / len;
            // Precompute a complex number with magnitude 1 representing a single step rotation
            // Multiplying by this applies one step of the rotation
            Complex wlen(cos(angle), sin(angle));
            // Iterate over each segment len, combining all pairs
            for (size_t i = 0; i < N; i += len) {
                Complex w(1);
                // Iterate through the first half of the segment, where j is the index within the segment
                // The second half of the segment is symmetrical and can be derived from the first half due to conjugate symmetry of the DFT
                // Compute the first half of the segment and use it to compute the second half with twiddle factors
                for (size_t j = 0; j < len / 2; ++j) {
                    // Combine odd and even halfs
                    // Even component
                    Complex u = buffer[i + j];
                    // Apply twiddle factor to odd component
                    Complex v = buffer[i + j + len / 2] * w;
                    // Update the value at the current index with the sum of even and odd components
                    buffer[i + j] = u + v;
                    // Update the value at the index offset by len/2 with the symmetrical opposite
                    buffer[i + j + len / 2] = u - v;
                    // Update the twiddle factor for the next position
                    w *= wlen;
                }
            }
        }
    }

    /// @brief Perform Fast Fourier Transform (FFT) on a real-valued vector.
    /// @param x Input vector of real numbers (floats). Must be a length that is a power of 2.
    /// @return Vector of complex numbers representing the FFT result.
    std::vector<Complex> fft(const std::vector<float>& buffer) {
        // Convert real-valued input to complex-valued array
        CArray data(buffer.size());
        transform(buffer.begin(), buffer.end(), data.begin(),
                  [](float val) { return Complex(static_cast<double>(val), 0); });
        // Perform FFT
        fft(data);
        // Return the result as a vector of complex numbers
        return vector<Complex>(data.begin(), data.end());
    }

    // Convert FFT spectrum to magnitude
    /// @brief Convert FFT spectrum to magnitude.
    /// @param spectrum Input complex spectrum from FFT.
    /// @return Vector of magnitudes corresponding to the input spectrum.
    std::vector<double> magnitude(const CArray& spectrum) {
        vector<double> mags(spectrum.size());
        transform(spectrum.begin(), spectrum.end(), mags.begin(),
                  [](const Complex& val) { return abs(val); });
        return mags;
    }

    // Convert FFT spectrum to dB scale
    /// @brief Convert FFT spectrum to dB scale.
    /// @param spectrum Input complex spectrum from FFT.
    /// @return Vector of magnitudes in dB corresponding to the input spectrum.
    std::vector<double> power(const CArray& spectrum) {
        std::vector<double> mags = magnitude(spectrum);
        std::transform(mags.begin(), mags.end(), mags.begin(),
                    [](double val) { return 20.0 * std::log10(val); });
        return mags;
    }
}