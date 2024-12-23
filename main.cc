#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include <chrono>
#include <vector>
#include "xxHash64.h"
#include <cstdint>
#include <bitset>
#include <bit>
#include <random>
#include <ctime>
#include <cmath>

int recordinality_no_hash(std::istream& stream, int k, int seed) {
    std::string word;
    std::vector<std::string> currentMaxVals(k, "");
    int R = 0;

    auto start = std::chrono::high_resolution_clock::now();
    while (stream >> word){
        std::string wordCopy = word;
        for (int i = 0; i < k; i++)
            if (word > currentMaxVals[i])
                std::swap(word, currentMaxVals[i]);
            else if (word == currentMaxVals[i])
                break;
        if (word != wordCopy)
            R++;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    return k*pow(1+1.0/k, R-k+1)-1;
}

int recordinality(std::istream& stream, int k, int seed) {
    std::string word;
    std::vector<uint64_t> currentMaxVals(k, 0);
    int R = 0;

    auto start = std::chrono::high_resolution_clock::now();
    while (stream >> word){
        uint64_t hashValue = XXHash64::hash(word.data(), word.size(), seed);
        uint64_t hashValueCopy = hashValue;
        for (int i = 0; i < k; i++)
            if (hashValue > currentMaxVals[i]) 
                std::swap(hashValue, currentMaxVals[i]);
            else if (hashValue == currentMaxVals[i])
                break;
        if (hashValue != hashValueCopy)
            R++;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    return k*pow(1+1.0/k, R-k+1)-1;
}

int minCount(std::istream& stream, int log_2_m, int seed) {
    std::string word;
    uint64_t m = std::pow(2, log_2_m);
    std::vector<long double> streamMins(m, 1.0);

    auto start = std::chrono::high_resolution_clock::now();
    while (stream >> word){
        uint64_t hashValue = XXHash64::hash(word.data(), word.size(), seed);
        long double floatHash = hashValue/pow(2.0, 64);

        auto i = hashValue%m;
        streamMins[i] = std::min(streamMins[i], floatHash);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    float minSums = 0;
    for (auto stmin : streamMins) minSums += stmin;
    return m*(m-1)/minSums;
}

int kMinVals(std::istream& stream, int k, int seed) {
    std::string word;
    std::vector<long double> minValues(k, 1.0);
    auto start = std::chrono::high_resolution_clock::now();
    while (stream >> word){
        long double hashValue = XXHash64::hash(word.data(), word.size(), seed)/(pow(2.0, 64));
        for (int i = 0; i < k; i++)
            if (hashValue < minValues[i])
                std::swap(hashValue, minValues[i]);
            else if (hashValue == minValues[i])
                break;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    return int((k-1)/minValues[k-1]);
}

int hyperLogLog(std::istream& stream, int log_2_m, int seed) {
    std::string word;
    uint64_t m = std::pow(2, log_2_m);
    float alpha = 0.7213/(1+1.079/m);
    if (m == 64) alpha = 0.709;
    else if (m == 32) alpha = 0.697;
    else if (m == 16) alpha = 0.673;
    std::vector<uint64_t> bmaps(m, 0-1);
    
    auto start = std::chrono::high_resolution_clock::now();
    while (stream >> word){
        uint64_t hashValue = XXHash64::hash(word.data(), word.size(), seed);
        uint64_t firstOne = hashValue & (-hashValue);
        auto i = hashValue>>(64-log_2_m);
        bmaps[i] &= -(firstOne<<1);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    float harm_mean = 0;
    for (int bmap : bmaps) harm_mean += 1.0/(bmap & (-bmap));
    return alpha*(m*m/harm_mean);
}

int PCSA(std::istream& stream, int seed) {
    std::string word;
    uint64_t log_2_m = 6;
    uint64_t m = std::pow(2, log_2_m);
    std::vector<uint64_t> bmaps(m, 0-1);

    auto start = std::chrono::high_resolution_clock::now();
    while (stream >> word){
        uint64_t hashValue = XXHash64::hash(word.data(), word.size(), seed);
        uint64_t firstOne = hashValue & (-hashValue);

        auto i = hashValue>>(64-log_2_m);
        bmaps[i] &= ~(firstOne);
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> duration = end - start;
    float R_mean = 0;
    for (int bmap : bmaps) {
        R_mean += std::__countr_zero(bmap);
    }
    return m*int(powf(2, R_mean/m))/0.77351;
}


int probabilisticCounting(std::istream& stream, int seed) {
    std::string word;
    uint64_t bmap = 0-1;
    auto start = std::chrono::high_resolution_clock::now();
    while (stream >> word){
        uint64_t hashValue = XXHash64::hash(word.data(), word.size(), seed);
        uint64_t firstOne = hashValue & (-hashValue);
        bmap &= ~(firstOne);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    return int(float(
        (bmap & (-bmap))
    )/0.77351);
}

int dummy(std::istream& stream, int seed) {
    std::string word;
    auto start = std::chrono::high_resolution_clock::now();
    int i = 0;
    while (stream >> word){i++;}
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    return 0;
}

int main() {
    srand(time(0));
    const std::string directory = "datasets";

    // Iterate over all files in the directory
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_regular_file()) { // Only process regular files
            const std::string filePath = entry.path().string();
            const std::string extension = entry.path().extension().string();

            if (extension == ".txt") {
                std::ifstream file(filePath);
                std::stringstream buffer;
                buffer << file.rdbuf();
                std::string fileContents = buffer.str();

                std::cout << "Processing data file: " << filePath << "\n";
                std::string datFile = filePath.substr(0, filePath.find_last_of('.'))+".dat";;
                std::ifstream resultFile(datFile);
                
                int true_card = 0;
                std::string word;
                while (resultFile >> word) true_card++;
                true_card /= 2;

                int estim_card = 0;
                int error_average = 0;
                // Open the file for reading
                int reps = 1;
                // Process the file

                auto start = std::chrono::high_resolution_clock::now();
                for (int i = 0; i < reps; i++) {
                    int seed = rand();
                    std::istringstream simstream(fileContents);
                    //int estim = 0;
                    //int estim = probabilisticCounting(simstream, seed);
                    //int estim = PCSA(simstream, seed);
                    //int estim = hyperLogLog(simstream, 9, seed);
                    //int estim = kMinVals(simstream, 128, seed);
                    //int estim = minCount(simstream, 8, seed);
                    //int estim = recordinality(simstream, 512, seed);
                    int estim = recordinality_no_hash(simstream, 512, seed);
                    //int estim = dummy(simstream, seed);
                    error_average += std::abs(estim-true_card);
                    estim_card += estim;
                }
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = end - start;
                //std::cout << "Time spent : " << duration.count() << " seconds\n";

                std::cout << "Finished processing: " << filePath << ", in " << duration.count() << " seconds" << "\n";
                std::cout 
                    << "True cardinality: " << true_card 
                    << ", Mean of Estimated cardinalities: " << estim_card/reps 
                    << ", Mean error: " << error_average/reps
                    << "\n\n";
            }
        }
    }
    
    /*
    int estim_card = 0;
    int error_average = 0;
    int true_card = 16;
    // Open the file for reading
    int reps = 10000;
    // Process the file
    for (int i = 0; i < reps; i++) {
        int seed = rand();
        std::istringstream simulatedStream("stream strean s d f e f e s f a s f s g hg r h v d hgjh er  s a");
        int estim = probabilisticCounting(simulatedStream, seed);
        //int estim = hyperLogLog(simulatedStream, seed);
        error_average += std::abs(estim-true_card);
        estim_card += estim;
    }
    
    std::cout 
        << "True cardinality: " << 1 
        << ", Mean of Estimated cardinalities: " << estim_card/reps 
        << ", Mean error: " << error_average/reps
        << "\n\n";
    */
    return 0;
}
