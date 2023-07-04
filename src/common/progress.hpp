#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <atomic>
#include <thread>
#include <chrono>
#include <iomanip>
#include <sstream>

namespace progress_meter {

class ProgressMeter {
public:
    std::string banner;
    std::atomic<uint64_t> total;
    std::atomic<uint64_t> completed;
    std::chrono::time_point<std::chrono::steady_clock> start_time;
    std::thread logger;
    ProgressMeter(uint64_t _total, const std::string& _banner)
        : total(_total), banner(_banner) {
        start_time = std::chrono::steady_clock::now();
        completed = 0;
        logger = std::thread(
            [&](void) {
                do_print();
                auto last = 0;
                while (completed < total) {
                    auto curr = completed - last;
                    if (curr > 0) {
                        do_print();
                        last = completed;
                    }
                    std::this_thread::sleep_for(std::chrono::milliseconds(500));
                }
            });
    };
    void do_print(void) {
        auto curr = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = curr-start_time;
        double rate = completed / elapsed_seconds.count();
        double seconds_to_completion = (completed > 0 ? (total - completed) / rate : 0);
        std::cerr << "\r" << banner << " "
                  << std::defaultfloat
                  << std::setfill(' ')
                  << std::setw(5)
                  << std::fixed
                  << std::setprecision(2)
                  << 100.0 * ((double)completed / (double)total) << "%"
                  << " @ "
                  << std::setw(4) << std::scientific << rate << " bp/s "
                  << "elapsed: " << print_time(elapsed_seconds.count()) << " "
                  << "remain: " << print_time(seconds_to_completion);
    }
    void finish(void) {
        completed.store(total);
        logger.join();
        do_print();
        std::cerr << std::endl;
    }
    std::string print_time(const double& _seconds) {
        int days = 0, hours = 0, minutes = 0, seconds = 0;
        distribute_seconds(days, hours, minutes, seconds, _seconds);
        std::stringstream buffer;
        buffer << std::setfill('0') << std::setw(2) << days << ":"
               << std::setfill('0') << std::setw(2) << hours << ":"
               << std::setfill('0') << std::setw(2) << minutes << ":"
               << std::setfill('0') << std::setw(2) << seconds;
        return buffer.str();
    }
    void distribute_seconds(int& days, int& hours, int& minutes, int& seconds, const double& input_seconds) {
        const int cseconds_in_day = 86400;
        const int cseconds_in_hour = 3600;
        const int cseconds_in_minute = 60;
        const int cseconds = 1;
        days = std::floor(input_seconds / cseconds_in_day);
        hours = std::floor(((int)input_seconds % cseconds_in_day) / cseconds_in_hour);
        minutes = std::floor((((int)input_seconds % cseconds_in_day) % cseconds_in_hour) / cseconds_in_minute);
        seconds = ((((int)input_seconds % cseconds_in_day) % cseconds_in_hour) % cseconds_in_minute) / cseconds; // + (input_seconds - std::floor(input_seconds));
        //std::cerr << input_seconds << " seconds is " << days << " days, " << hours << " hours, " << minutes << " minutes, and " << seconds << " seconds." << std::endl;
    }
    void increment(const uint64_t& incr) {
        completed += incr;
    }
};

}
