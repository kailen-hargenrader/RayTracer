#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <cmath>
#include <string>
#include <cstdlib>

#include "../../s2899221/Code/raytracer.h"

static double mean(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double s = std::accumulate(v.begin(), v.end(), 0.0);
    return s / static_cast<double>(v.size());
}

static double stdev(const std::vector<double>& v, double mu) {
    if (v.size() < 2) return 0.0;
    double sumsq = 0.0;
    for (double x : v) { double d = x - mu; sumsq += d * d; }
    return std::sqrt(sumsq / static_cast<double>(v.size() - 1));
}

int main(int argc, char** argv) {
    const std::string json_path = (argc > 1) ? argv[1] : std::string("s2899221/ASCII/50_rand_1_cam.json");
    const std::string camera_id = (argc > 2) ? argv[2] : std::string("1");
    const int trials = (argc > 3) ? std::max(1, std::atoi(argv[3])) : 10;

    std::vector<double> times_unacc; times_unacc.reserve(trials);
    std::vector<double> times_acc;   times_acc.reserve(trials);

    for (int i = 0; i < trials; ++i) {
        double t0 = time_trace_all_one_pass(json_path, camera_id);
        times_unacc.push_back(t0);
    }
    for (int i = 0; i < trials; ++i) {
        double t1 = time_trace_all_accelerated_one_pass(json_path, camera_id);
        times_acc.push_back(t1);
    }

    const double mu_un = mean(times_unacc);
    const double sd_un = stdev(times_unacc, mu_un);
    const double mu_ac = mean(times_acc);
    const double sd_ac = stdev(times_acc, mu_ac);

    // Write timings file in same folder as this executable source
    std::ofstream out("tests/raytracer/timing.txt");
    if (!out) {
        std::cerr << "Failed to open output file tests/raytracer/timing.txt\n";
        return 1;
    }
    out << "unaccelerated_times\n";
    for (double t : times_unacc) out << t << "\n";
    out << "accelerated_times\n";
    for (double t : times_acc) out << t << "\n";
    out << "mean_unaccelerated\n" << mu_un << "\n";
    out << "std_unaccelerated\n" << sd_un << "\n";
    out << "mean_accelerated\n" << mu_ac << "\n";
    out << "std_accelerated\n" << sd_ac << "\n";
    out.close();

    std::cout << "Wrote tests/raytracer/timing.txt\n";
    std::cout << "Mean (unaccelerated): " << mu_un << ", std: " << sd_un << "\n";
    std::cout << "Mean (accelerated):   " << mu_ac << ", std: " << sd_ac << "\n";
    return 0;
}


