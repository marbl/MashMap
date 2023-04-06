/**
 * @file    align_utils.hpp
 * @brief   Critical type defintions for generating alignments
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef ALIGN_UTILS_HPP 
#define ALIGN_UTILS_HPP

#include <vector>
#include <cstdint>
#include <algorithm>
#include <string>
#include <cstring>
#include <memory>

namespace align
{
  char* alignment_to_cigar(
        const std::string& edit_cigar,
        int& target_aligned_length,
        int& query_aligned_length,
        int& matches,
        int& mismatches,
        int& insertions,
        int& inserted_bp,
        int& deletions,
        int& deleted_bp) {
    // the edit cigar contains a character string of ops
    // here we compress them into the standard cigar representation
    const int start_idx = 0;
    const int end_idx = edit_cigar.size() - 1;

    auto *cigar = new std::vector<char>();
    char lastMove = 0; // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;

    // std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int i = start_idx; i <= end_idx; i++) {
        // if new sequence of same moves started
        if (i == end_idx || (edit_cigar[i] != lastMove && lastMove != 0)) {
            // calculate matches, mismatches, insertions, deletions
            switch (lastMove) {
                case 'M':
                    matches += numOfSameMoves;
                    query_aligned_length += numOfSameMoves;
                    target_aligned_length += numOfSameMoves;
                    break;
                case 'X':
                    mismatches += numOfSameMoves;
                    query_aligned_length += numOfSameMoves;
                    target_aligned_length += numOfSameMoves;
                    break;
                case 'I':
                    ++insertions;
                    inserted_bp += numOfSameMoves;
                    query_aligned_length += numOfSameMoves;
                    break;
                case 'D':
                    ++deletions;
                    deleted_bp += numOfSameMoves;
                    target_aligned_length += numOfSameMoves;
                    break;
                default:
                    break;
            }

            // Write number of moves to cigar string.
            int numDigits = 0;
            for (; numOfSameMoves; numOfSameMoves /= 10) {
                cigar->push_back('0' + numOfSameMoves % 10);
                numDigits++;
            }
            std::reverse(cigar->end() - numDigits, cigar->end());
            // Write code of move to cigar string.
            // reassign 'M' to '=' for convenience
            lastMove = lastMove == 'M' ? '=' : lastMove;
            cigar->push_back(lastMove);
            // If not at the end, start new sequence of moves.
            if (i < end_idx) {
                numOfSameMoves = 0;
            }
        }
        if (i < end_idx) {
            lastMove = edit_cigar[i];
            numOfSameMoves++;
        }
    }
    cigar->push_back(0); // Null character termination.

    char *cigar_ = (char *)malloc(cigar->size() * sizeof(char));
    std::memcpy(cigar_, &(*cigar)[0], cigar->size() * sizeof(char));
    delete cigar;

    return cigar_;
  }
}

#endif
