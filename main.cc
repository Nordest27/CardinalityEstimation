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
#include <tuple>

std::tuple<int, std::chrono::duration<double>> recordinality_no_hash(std::istream& stream, int k) {
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

return {k*pow(1+1.0/k, R-k+1)-1, duration};
}

std::tuple<int, std::chrono::duration<double>> recordinality(std::istream& stream, int k, int seed) {
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

    return {k*pow(1+1.0/k, R-k+1)-1, duration};
}

std::tuple<int, std::chrono::duration<double>> minCount(std::istream& stream, int log_2_m, int seed) {
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
    return {m*(m-1)/minSums, duration};
}

std::tuple<int, std::chrono::duration<double>> kMinVals(std::istream& stream, int k, int seed) {
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
    return {int((k-1)/minValues[k-1]), duration};
}

std::tuple<int, std::chrono::duration<double>> hyperLogLog(std::istream& stream, int log_2_m, int seed) {
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
    return {alpha*(m*m/harm_mean), duration};
}

std::tuple<int, std::chrono::duration<double>> PCSA(std::istream& stream, uint64_t log_2_m, int seed) {
    std::string word;
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
    return {m*int(powf(2, R_mean/m))/0.77351, duration};
}


std::tuple<int, std::chrono::duration<double>> probabilisticCounting(std::istream& stream, int seed) {
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
    return {int(float((bmap & (-bmap)))/0.77351), duration};
}

std::tuple<int, std::chrono::duration<double>> dummy(std::istream& stream, int seed) {
    std::string word;
    auto start = std::chrono::high_resolution_clock::now();
    int i = 0;
    while (stream >> word){i++;}
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    return {0, duration};
}

std::string generateText(int n, int N, float alpha, int seed) {
    std::mt19937 gen(std::random_device{}());
    std::string result = "";
    std::vector<float> probs = std::vector<float>(n, 0);
    float cn_sum = 0;
    for (int i = 1; i <= n; ++i) cn_sum += powf(i, -alpha);
    float cn = 1.0/cn_sum;
    for (int i = 1; i <= n; ++i) {
        probs[i] = cn/(powf(i, alpha));
    }
    std::discrete_distribution<std::size_t> dist{probs.begin(), probs.end()};
    for (int i = 0; i < N; ++i)
        result.append(std::to_string(dist(gen)) + " ");
    return result;
}

int main() {
    srand(time(0));
    const std::string directory = "datasets";

    int reps = 100;
    int m = 512;
    // Iterate over all files in the directory
    std::ofstream output("all-files-results-k-"+std::to_string(m)+".csv"), times("all-files-times-k-"+std::to_string(m)+".csv");
    output<<"Card, PC, PCSA, HLL, KMV, MC, RC"<<std::endl;
    times<<"Read, PC, PCSA, HLL, KMV, MC, RC"<<std::endl;
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

                auto start = std::chrono::high_resolution_clock::now();
                
                for (int i = 0; i < reps; i++) {
                    int seed = rand();
                    output << filePath << ", ";
                    times << filePath << ", ";

                    // Declare `simstream` once, then reset its content for each function.
                    std::istringstream simstream(fileContents);

                    // Dummy
                    simstream.str(fileContents); simstream.clear();
                    auto [dummy_estim, dummy_duration] = dummy(simstream, seed);
                    output << true_card << ", ";
                    times << dummy_duration.count() << ", ";

                    // Probabilistic Counting
                    simstream.str(fileContents); simstream.clear();
                    auto [pc_estim, pc_duration] = probabilisticCounting(simstream, seed);
                    output << pc_estim << ", ";
                    times << pc_duration.count() << ", ";

                    // PCSA
                    simstream.str(fileContents); simstream.clear();
                    auto [pcsa_estim, pcsa_duration] = PCSA(simstream, std::log2(m), seed);
                    output << pcsa_estim << ", ";
                    times << pcsa_duration.count() << ", ";

                    // HyperLogLog
                    simstream.str(fileContents); simstream.clear();
                    auto [hll_estim, hll_duration] = hyperLogLog(simstream, std::log2(m), seed);
                    output << hll_estim << ", ";
                    times << hll_duration.count() << ", ";

                    // K-Min Values
                    simstream.str(fileContents); simstream.clear();
                    auto [kmv_estim, kmv_duration] = kMinVals(simstream, m, seed);
                    output << kmv_estim << ", ";
                    times << kmv_duration.count() << ", ";

                    // MinCount
                    simstream.str(fileContents); simstream.clear();
                    auto [mc_estim, mc_duration] = minCount(simstream, std::log2(m), seed);
                    output << mc_estim << ", ";
                    times << mc_duration.count() << ", ";

                    // Recordinality
                    simstream.str(fileContents); simstream.clear();
                    auto [rc_estim, rc_duration] = recordinality(simstream, m, seed);
                    output << rc_estim << "\n";
                    times << rc_duration.count() << "\n";

                    error_average += std::abs(rc_estim - true_card);
                    estim_card += rc_estim;
                }
                /*
                // Recordinality No Hash (separate declaration as it doesn't need `seed`)
                std::istringstream simstream(fileContents);
                auto [rcnh_estim, rcnh_duration] = recordinality_no_hash(simstream, m);
                
                output << filePath << ", ";
                times << filePath << ", ";
                output << rcnh_estim << "\n";
                times << rcnh_duration.count() << "\n";
                */
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
    std::ofstream output_syn("synthetic-"+std::to_string(m)+".csv"), times_syn("synthetic-times-"+std::to_string(m)+".csv");
    output_syn<<"n, N, PC, PCSA, HLL, KMV, MC, RC"<<std::endl;
    times_syn<<"n, N, Read, PC, PCSA, HLL, KMV, MC, RC"<<std::endl;
    for (int n=1000; n <= 10000; n+=1000)
        for (int N=n*2; N <= 10*n; N += n) {
            std::string fileContents=generateText(n, N, 0, 0);
            for (int i = 0; i < reps; i++) {
                int seed = rand();

                // Declare `simstream` once, then reset its content for each function.
                std::istringstream simstream(fileContents);

                // Dummy
                simstream.str(fileContents); simstream.clear();
                auto [dummy_estim, dummy_duration] = dummy(simstream, seed);
                output_syn << n << ", " << N << ", ";
                times_syn << n << ", " << N << ", " << dummy_duration.count() << ", ";

                // Probabilistic Counting
                simstream.str(fileContents); simstream.clear();
                auto [pc_estim, pc_duration] = probabilisticCounting(simstream, seed);
                output_syn << pc_estim << ", ";
                times_syn << pc_duration.count() << ", ";

                // PCSA
                simstream.str(fileContents); simstream.clear();
                auto [pcsa_estim, pcsa_duration] = PCSA(simstream, std::log2(m), seed);
                output_syn << pcsa_estim << ", ";
                times_syn << pcsa_duration.count() << ", ";

                // HyperLogLog
                simstream.str(fileContents); simstream.clear();
                auto [hll_estim, hll_duration] = hyperLogLog(simstream, std::log2(m), seed);
                output_syn << hll_estim << ", ";
                times_syn << hll_duration.count() << ", ";

                // K-Min Values
                simstream.str(fileContents); simstream.clear();
                auto [kmv_estim, kmv_duration] = kMinVals(simstream, m, seed);
                output_syn << kmv_estim << ", ";
                times_syn << kmv_duration.count() << ", ";

                // MinCount
                simstream.str(fileContents); simstream.clear();
                auto [mc_estim, mc_duration] = minCount(simstream, std::log2(m), seed);
                output_syn << mc_estim << ", ";
                times_syn << mc_duration.count() << ", ";

                // Recordinality
                simstream.str(fileContents); simstream.clear();
                auto [rc_estim, rc_duration] = recordinality(simstream, m, seed);
                output_syn << rc_estim << "\n";
                times_syn << rc_duration.count() << "\n";
            }
        }
    return 0;
}
