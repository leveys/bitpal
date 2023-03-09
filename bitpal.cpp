#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>

#define M 1
#define I -1
#define G -3

#define MIN G
#define MID I - G
#define MAX M - G

using namespace std;

class MatchVectors {
    private:
        uint64_t matchA;
        uint64_t matchC;
        uint64_t matchG;
        uint64_t matchT;
    public:
        MatchVectors(const string& seq) {
            for (int i = 0; i < seq.length(); i++) {
                // set ith bit of matchvector
                get(seq[i]) |= 1 << i;
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
};

void align_bitpal(const string& X, const string& Y) {
    // create match vectors, these are created from the horizontal sequence
    MatchVectors matches = MatchVectors(Y);

    // create deltaV and deltaH bitvectors
    uint64_t deltaV[MAX - MIN + 1];
    uint64_t deltaH[MAX - MIN + 1];

    // initialize deltaH_min to all ones
    deltaH[MIN] = -1;

    // main loop
    for (char c: X) {
        uint64_t curr_matches = matches.get(c);

        // step 1
        uint64_t initmax = deltaH[MIN] & curr_matches;
        deltaV[MAX] = ((initmax + deltaH[MIN]) ^ deltaH[MIN]) ^ initmax;

        // step 2
        uint64_t DH_min_remaining = deltaH[MIN] ^ (deltaV[MAX] >> 1);
        for (int i = MAX-1; i > MID; i--) {
            uint64_t initcurrent = deltaH[MIN+MAX-i] & deltaV[MAX];
            for (int j = i; j < MAX; j++) {
                initcurrent |= deltaH[MIN+MAX-j] & (deltaV[j] & ~curr_matches);
            }
            deltaV[i] = ((initcurrent << 1) + DH_min_remaining) ^ DH_min_remaining;
        }

        // step 3
        uint64_t not_DV_high;
        for (int i = MID; i <= MAX; )
        for (int i = MID; i >= MIN; i--) {

        }
    }
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

    align_bitpal(X, Y);

    return EXIT_SUCCESS;
}