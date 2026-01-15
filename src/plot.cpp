#include <matplot/matplot.h>
#include <vector>

void plotFFT(const std::vector<double>& freqScale, const std::vector<double>& magnitudes) {
    using namespace matplot;

    // Create a plot
    bar(freqScale, magnitudes);
    
    // This displays the window and stays open until you close it
    show(); 
    return;
}