/**
 * Seed-and-Extend Algorithm using K-mer Hashing
 * 
 * This algorithm is used for fast sequence alignment by:
 * 1. Finding exact k-mer matches (seeds) using hash tables
 * 2. Extending seeds to find longer alignments
 * 
 * This is the foundational approach used in BLAST, MAQ, and SOAP.
 */

#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

struct Seed {
    int query_pos;
    int ref_pos;
};

struct AlignmentHit {
    int query_start;
    int query_end;
    int ref_start;
    int ref_end;
    int score;
    string query_seq;
    string ref_seq;
};

/**
 * Build k-mer index from reference sequence.
 */
unordered_map<string, vector<int>> build_kmer_index(const string& sequence, int k) {
    unordered_map<string, vector<int>> index;
    
    for (size_t i = 0; i <= sequence.length() - k; i++) {
        string kmer = sequence.substr(i, k);
        index[kmer].push_back(i);
    }
    
    return index;
}

/**
 * Find exact k-mer matches between query and reference.
 */
vector<Seed> find_seeds(const string& query, 
                       const unordered_map<string, vector<int>>& index, 
                       int k) {
    vector<Seed> seeds;
    
    for (size_t i = 0; i <= query.length() - k; i++) {
        string kmer = query.substr(i, k);
        
        auto it = index.find(kmer);
        if (it != index.end()) {
            for (int ref_pos : it->second) {
                seeds.push_back({(int)i, ref_pos});
            }
        }
    }
    
    return seeds;
}

/**
 * Extend a seed match in both directions.
 */
AlignmentHit extend_seed(const string& seq1, const string& seq2,
                        int seed_pos1, int seed_pos2,
                        int match_score = 1, int mismatch_penalty = -1) {
    int len1 = seq1.length();
    int len2 = seq2.length();
    
    // Extend right
    int score = 0;
    int max_score = 0;
    int max_right = 0;
    
    int i = 0;
    while (seed_pos1 + i < len1 && seed_pos2 + i < len2) {
        if (seq1[seed_pos1 + i] == seq2[seed_pos2 + i]) {
            score += match_score;
        } else {
            score += mismatch_penalty;
        }
        
        if (score > max_score) {
            max_score = score;
            max_right = i + 1;
        }
        
        // Stop if score drops too much below max
        if (score < max_score - 5) {
            break;
        }
        
        i++;
    }
    
    // Extend left
    score = max_score;
    int max_left = 0;
    
    i = 1;
    while (seed_pos1 - i >= 0 && seed_pos2 - i >= 0) {
        if (seq1[seed_pos1 - i] == seq2[seed_pos2 - i]) {
            score += match_score;
        } else {
            score += mismatch_penalty;
        }
        
        if (score > max_score) {
            max_score = score;
            max_left = i;
        }
        
        // Stop if score drops too much below max
        if (score < max_score - 5) {
            break;
        }
        
        i++;
    }
    
    int start1 = seed_pos1 - max_left;
    int end1 = seed_pos1 + max_right;
    int start2 = seed_pos2 - max_left;
    int end2 = seed_pos2 + max_right;
    
    return {
        start1, end1, start2, end2, max_score,
        seq1.substr(start1, end1 - start1),
        seq2.substr(start2, end2 - start2)
    };
}

/**
 * Perform seed-and-extend alignment between query and reference.
 */
vector<AlignmentHit> seed_and_extend(const string& reference, 
                                    const string& query, 
                                    int k = 11) {
    // Build k-mer index
    auto index = build_kmer_index(reference, k);
    
    // Find seeds
    vector<Seed> seeds = find_seeds(query, index, k);
    
    // Extend each seed
    vector<AlignmentHit> alignments;
    for (const auto& seed : seeds) {
        AlignmentHit hit = extend_seed(query, reference, seed.query_pos, seed.ref_pos);
        alignments.push_back(hit);
    }
    
    // Sort by score descending
    sort(alignments.begin(), alignments.end(), 
         [](const AlignmentHit& a, const AlignmentHit& b) {
             return a.score > b.score;
         });
    
    return alignments;
}

int main() {
    // Example usage
    string reference = "ACGTACGTACGTAAACCCGGGTTTACGTACGT";
    string query = "ACGTAAACCCGGG";
    
    cout << "Seed-and-Extend Alignment" << endl;
    cout << "Reference: " << reference << endl;
    cout << "Query: " << query << endl;
    cout << "\nK-mer length: 5" << endl;
    
    vector<AlignmentHit> alignments = seed_and_extend(reference, query, 5);
    
    cout << "\nFound " << alignments.size() << " alignment(s):" << endl;
    
    // Show top 5
    for (size_t i = 0; i < min(size_t(5), alignments.size()); i++) {
        const auto& aln = alignments[i];
        cout << "\nAlignment " << (i + 1) << ":" << endl;
        cout << "  Query position: " << aln.query_start << "-" << aln.query_end << endl;
        cout << "  Reference position: " << aln.ref_start << "-" << aln.ref_end << endl;
        cout << "  Score: " << aln.score << endl;
        cout << "  Query sequence: " << aln.query_seq << endl;
        cout << "  Ref sequence:   " << aln.ref_seq << endl;
    }
    
    return 0;
}
