#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <bit>
#include <bitset>

using namespace std;

const uint64_t one = 1;
const uint64_t all_ones = -1;

class MatchVectors {
    private:
        string seq;
        uint64_t matchA = 0;
        uint64_t matchC = 0;
        uint64_t matchG = 0;
        uint64_t matchT = 0;
    public:
        MatchVectors(const string& seq) {
            this->seq = seq;
            for (int i = 0; i < seq.length(); i++) {
                // set corresponding bit of matchvector
                // first character corresponds to lowest order bit of the bitvector
                get(seq[i]) |= one << i;
            }
        }

        uint64_t& get(char c) {
            // get the match bitvector for the given character
            switch (c) {
                case 'A':
                    return matchA;
                case 'C':
                    return matchC;
                case 'G':
                    return matchG;
                case 'T':
                    return matchT;
                default:
                    throw invalid_argument("only characters 'A', 'C', 'G', and 'T' are accepted");
                    exit(EXIT_FAILURE);
            }
        }

        friend ostream& operator<<(ostream& os, MatchVectors& mv) {
            os << "sequence:       " << mv.seq << endl;
            os << "match vector A: " << bitset<64>(mv.matchA) << endl;
            os << "match vector C: " << bitset<64>(mv.matchC) << endl;
            os << "match vector G: " << bitset<64>(mv.matchG) << endl;
            os << "match vector T: " << bitset<64>(mv.matchT) << endl << endl;
            return os;
        }
};

void print_bitvec(uint64_t bv, string name) {
    cout << bitset<64>(bv) << " " << name << endl;
}

int align_bitpal(const string& X, const string& Y) {
    // create match vectors, these are created from the horizontal sequence
    MatchVectors match_vectors = MatchVectors(Y);

    // create deltaV and deltaH bitvectors
    uint64_t deltaV4 = 0;
    uint64_t deltaV3 = 0;
    uint64_t deltaV2 = 0;
    uint64_t deltaV1 = 0;
    uint64_t deltaV0 = 0;
    uint64_t deltaVmin1 = 0;
    uint64_t deltaVmin2 = 0;
    uint64_t deltaVmin3 = 0;

    uint64_t deltaH4 = 0;
    uint64_t deltaH3 = 0;
    uint64_t deltaH2 = 0;
    uint64_t deltaH1 = 0;
    uint64_t deltaH0 = 0;
    uint64_t deltaHmin1 = 0;
    uint64_t deltaHmin2 = 0;
    uint64_t deltaHmin3 = all_ones;

    // main loop
    for (char c: X) {
        uint64_t matches = match_vectors.get(c);

        // step 1, compute deltaV_max
        uint64_t deltaV4_initial = deltaHmin3 & matches;
        deltaV4 = ((deltaV4_initial + deltaHmin3) ^ deltaHmin3) ^ deltaV4_initial;
        uint64_t deltaV4_or_match = deltaV4 | matches;

        // step 2, compute remaining deltaV_high (only one in this case)
        uint64_t deltaHmin3_remaining = deltaHmin3 ^ (deltaV4 >> 1);
        deltaV3 = (((deltaHmin2 & deltaV4_or_match) << 1) + deltaHmin3_remaining) ^ deltaHmin3_remaining;
        uint64_t deltaV3_not_match = deltaV3 & ~matches;

        // step 3, compute deltaV_low
        uint64_t deltaV_low = ~(deltaV4_or_match | deltaV3);

        deltaV2 = ((deltaHmin1 & deltaV4_or_match) | (deltaHmin2 & deltaV3_not_match) | (deltaHmin3 & deltaV_low)) << 1;
        deltaV1 = ((deltaH0 & deltaV4_or_match) | (deltaHmin1 & deltaV3_not_match) | (deltaHmin2 & deltaV_low)) << 1;
        deltaV0 = ((deltaH1 & deltaV4_or_match) | (deltaH0 & deltaV3_not_match) | (deltaHmin1 & deltaV_low)) << 1;
        deltaVmin1 = ((deltaH2 & deltaV4_or_match) | (deltaH1 & deltaV3_not_match) | (deltaH0 & deltaV_low)) << 1;
        deltaVmin2 = ((deltaH3 & deltaV4_or_match) | (deltaH2 & deltaV3_not_match) | (deltaH1 & deltaV_low)) << 1;

        // step 4, compute deltaV_min
        deltaVmin3 = all_ones ^ (deltaV4 | deltaV3 | deltaV2 | deltaV1 | deltaV0 | deltaVmin1 | deltaVmin2);

        // step 5, compute deltaH values

        // convert all deltaH_low values to deltaH_mid
        deltaH2 |= deltaHmin3 | deltaHmin2 | deltaHmin1 | deltaH0 | deltaH1;
        deltaH1 = 0;
        deltaH0 = 0;
        deltaHmin1 = 0;
        deltaHmin2 = 0;
        deltaHmin3 = 0;
        // add match positions to previous deltaH_max
        deltaH4 |= matches;
        // remove match positions from all other previous deltaH vectors
        deltaH3 &= ~matches;
        deltaH2 &= ~matches;
        // all other vectors are 0

        // compute current deltaH values
        uint64_t deltaH_low = ~(deltaH4 | deltaH3);
        deltaHmin2 = (deltaH4 & deltaV3) | (deltaH3 & deltaV2) | (deltaH_low & deltaV1);
        deltaHmin1 = (deltaH4 & deltaV2) | (deltaH3 & deltaV1) | (deltaH_low & deltaV0);
        deltaH0 = (deltaH4 & deltaV1) | (deltaH3 & deltaV0) | (deltaH_low & deltaVmin1);
        deltaH1 = (deltaH4 & deltaV0) | (deltaH3 & deltaVmin1) | (deltaH_low & deltaVmin2);
        deltaH2 = (deltaH4 & deltaVmin1) | (deltaH3 & deltaVmin2) | (deltaH_low & deltaVmin3);
        deltaH3 = (deltaH4 & deltaVmin2) | (deltaH3 & deltaVmin3);
        deltaH4 = deltaH4 & deltaVmin3;
        deltaHmin3 = all_ones ^ (deltaH4 | deltaH3 | deltaH2 | deltaH1 | deltaH0 | deltaHmin1 | deltaHmin2);
    }

    uint64_t seq_len = ~(all_ones << Y.length());

    // compute final result
    return -3*X.length()    // sum of left column of deltaV matrix (always G)
        // sum of bottom row of deltaH matrix
        + popcount(deltaHmin3 & seq_len) * -3
        + popcount(deltaHmin2 & seq_len) * -2
        + popcount(deltaHmin1 & seq_len) * -1
        + popcount(deltaH1 & seq_len) * 1
        + popcount(deltaH2 & seq_len) * 2
        + popcount(deltaH3 & seq_len) * 3
        + popcount(deltaH4 & seq_len) * 4;
}


void readSequences(const string& filename, vector<string>& sequences) {
    ifstream ifs(filename.c_str());
    if (!ifs) {
        throw runtime_error("Could not open file: " + filename);
    }
    string line;
    while (ifs) {
        getline(ifs, line);
        if (line.empty()) {
            continue;
        }
        if (line.front() == '>') {
            sequences.push_back(string());
            continue;
        }
        sequences.back().append(line);
    }
}

int main(int argc, char** argv) {
    if (argc != 2) {
        cout << "Usage: bitpal input.fasta\n\n";
        return EXIT_FAILURE;
    }

    vector<string> sequences;
    readSequences(argv[1], sequences);

    if (sequences.size() != 2) {
        cerr << "Input FASTA file should contain only two sequences\n";
        return EXIT_FAILURE;
    }

    // X is vertical, Y is horizontal
    // we want the longest sequence to be the horizontal one, 
    // because we know that the two sequences are always 64 characters or less
    string X = sequences[0];
    string Y = sequences[1];
    if (X.length() > Y.length()) {
        // swap sequences
        string tmp = X;
        X = Y;
        Y = tmp;
    }

    int align_score = align_bitpal(X, Y);

    cout << "alignment score: " << align_score << endl;

    return EXIT_SUCCESS;
}