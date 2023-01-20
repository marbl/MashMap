//
// Created by guarracino on 28/07/21.
//

#include <algorithm>
#include "utils.hpp"

namespace mashmap {
    bool is_a_number(const std::string& s) {
        return !s.empty() && s.find_first_not_of("0123456789.") == std::string::npos && std::count(s.begin(), s.end(), '.') < 2;
    };

    int64_t handy_parameter(const std::string& value) {
        uint64_t str_len = value.length();
        uint8_t exp = 0;
        if (value[str_len-1] == 'k' || value[str_len-1] == 'K') {
            exp = 3;
            --str_len;
        } else if (value[str_len-1] == 'm' || value[str_len-1] == 'M') {
            exp = 6;
            --str_len;
        } else if (value[str_len-1] == 'g' || value[str_len-1] == 'G') {
            exp = 9;
            --str_len;
        }

        const std::string tmp = value.substr(0, str_len);
        return is_a_number(tmp) ? (int)(stof(tmp) * pow(10, exp)) : -1;
    }

}
