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


int hyperLogLog(std::istream& stream, int seed) {
    std::string word;
    uint64_t log_2_m = 10;
    uint64_t m = std::pow(2, log_2_m);
    float alpha = 0.7213/(1+1.079/m);
    if (m == 64) alpha = 0.709;
    else if (m == 32) alpha = 0.697;
    else if (m == 16) alpha = 0.673;

    std::vector<uint64_t> bmaps(std::pow(2, log_2_m), 0-1);

    //std::vector<int> Rs(std::pow(2, log_2_m), 0);
    
    auto start = std::chrono::high_resolution_clock::now();
    while (stream >> word){
        uint64_t hashValue = XXHash64::hash(word.data(), word.size(), seed);
        //std::hash<std::string> hasher;
        //auto hashValue = hasher(word);
        //std::cout << "hashValue: " << std::bitset<64>(hashValue) << "\n";
        //std::cout << "Count Right Zeros: " <<  std::__countr_zero(hashValue) << "\n";
        //std::cout << "firstOne : " << std::bitset<64>(hashValue & (-hashValue)) << "\n";
        //std::cout << "Neg      : " << std::bitset<64>(-(hashValue & (-hashValue))) << "\n";
        //std::cout << "Inv      : " << std::bitset<64>(~(hashValue & (-hashValue))) << "\n";
        
        uint64_t firstOne = hashValue & (-hashValue);

        // Update the bitmap directly
        auto i = hashValue>>(64-log_2_m);
        //Rs[i] = std::max(Rs[i], std::__countr_zero(hashValue));
        bmaps[i] &= -(firstOne<<1);
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> duration = end - start;
    
    //std::cout << "Bitmap             : " << std::bitset<64>(bmap) << "\n";
    //std::cout << "Bitmap first one   : " << std::bitset<64>(bmap & (-bmap)) << "\n";
    //std::cout << "2^R                : " << int(float(((bmap & (-bmap)) << 1))/0.77351) << "\n";
    //std::cout << "Time spent hashing : " << duration.count() << " seconds\n";
    //std::cout << "Number of words    : " << i << "\n";
    float R_harm_mean = 0;
    float harm_mean = 0;
    //for (auto R : Rs) R_harm_mean += 1.0/R; 
    for (int bmap : bmaps) {
        /*
        std::cout << "Bitmap             : " << std::bitset<64>(bmap) << "\n";
        std::cout << "Bitmap first one   : " << std::bitset<64>(bmap & (-bmap)) << "\n";
        std::cout << "Value              : " << (bmap & (-bmap)) << "\n";
        std::cout << "Inverse            : " << 1.0/(bmap & (-bmap)) << "\n";
        std::cout << "Count Right Zeros  : " << std::__countr_zero((bmap & (-bmap))) << "\n";
        */
       //std::cout<<R<<std::endl;
       //R_harm_mean += 1.0/(std::__countr_zero(bmap));
       harm_mean += 1.0/(bmap & (-bmap));
    }
    //std::cout<<"Sum: "<<R_harm_mean<<std::endl;
    //return m*int(pow(2, m/R_harm_mean));
    return alpha*(m*m/harm_mean);
}

int PCSA(std::istream& stream, int seed) {
    std::string word;
    uint64_t log_2_m = 6;
    uint64_t m = std::pow(2, log_2_m);

    std::vector<uint64_t> bmaps(std::pow(2, log_2_m), 0-1);
    
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
    uint64_t bmap_ones = 0;
    auto start = std::chrono::high_resolution_clock::now();
    while (stream >> word){
        uint64_t hashValue = XXHash64::hash(word.data(), word.size(), seed);
        //std::hash<std::string> hasher;
        //auto hashValue = hasher(word);
        //std::cout << "hashValue: " << std::bitset<64>(hashValue) << "\n";
        //std::cout << "firstOne : " << std::bitset<64>(hashValue & (-hashValue)) << "\n";
        //std::cout << "Inv shift: " << std::bitset<64>(~(hashValue & (-hashValue))) << "\n";
        uint64_t firstOne = hashValue & (-hashValue);
        // Update the bitmap directly
        bmap &= ~(firstOne);
        //bmap_ones |= firstOne;
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> duration = end - start;


    //std::cout << "Bitmap             : " << std::bitset<64>(bmap) << "\n";
    //std::cout << "Bitmap ones        : " << std::bitset<64>(bmap_ones) << "\n";
    //std::cout << "Bitmap first one   : " << std::bitset<64>(bmap & (-bmap)) << "\n";
    //std::cout << "2^R                : " << int(float(bmap & (-bmap))/0.77351) << "\n";
    //std::cout << "Time spent hashing : " << duration.count() << " seconds\n";
    //std::cout << "Number of words    : " << i << "\n";
    //std::cout << "count_right_ones: " << (std::pow(2, std::__countr_one(bmap_ones))) <<", value: " << (bmap & (-bmap))<< std::endl;
    return int(float(
        (bmap & (-bmap))
        //std::pow(2, std::__countr_one(bmap_ones))
    //));
    )/0.77351);
}

int dummy(std::istream& stream, int seed) {
    std::string word;
    auto start = std::chrono::high_resolution_clock::now();
    int i = 0;
    while (stream >> word){i++;}
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    //std::cout << "Time spent : " << duration.count() << " seconds\n";
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
                int reps = 100;
                // Process the file

                auto start = std::chrono::high_resolution_clock::now();
                for (int i = 0; i < reps; i++) {
                    int seed = rand();
                    std::ifstream file(filePath);
                    //int estim = 0;
                    int estim = probabilisticCounting(file, seed);
                    //int estim = PCSA(file, seed);
                    //int estim = hyperLogLog(file, seed);
                    //int estim = dummy(file, seed);
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
