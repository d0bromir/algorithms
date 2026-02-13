/**
 * Burrows-Wheeler Transform (BWT) and FM-Index
 * 
 * The BWT is a reversible transformation that reorganizes a string to make it more compressible.
 * The FM-index uses the BWT to enable fast exact pattern matching and is foundational to
 * BWA, Bowtie, and HISAT2 aligners.
 */

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>

using namespace std;

/**
 * Compute the Burrows-Wheeler Transform of a text.
 */
string burrows_wheeler_transform(string text) {
    // Add sentinel character if not present
    if (text.empty() || text.back() != '$') {
        text += '$';
    }
    
    int n = text.length();
    
    // Create all rotations
    vector<string> rotations;
    for (int i = 0; i < n; i++) {
        rotations.push_back(text.substr(i) + text.substr(0, i));
    }
    
    // Sort rotations lexicographically
    sort(rotations.begin(), rotations.end());
    
    // BWT is the last column of sorted rotations
    string bwt;
    for (const auto& rotation : rotations) {
        bwt += rotation.back();
    }
    
    return bwt;
}

/**
 * Reverse the Burrows-Wheeler Transform.
 */
string inverse_bwt(const string& bwt) {
    int n = bwt.length();
    
    // Create table with indices
    vector<pair<char, int>> table;
    for (int i = 0; i < n; i++) {
        table.push_back({bwt[i], i});
    }
    
    // Sort to get first column
    sort(table.begin(), table.end());
    
    // Follow the links to reconstruct
    string result;
    int idx = 0;  // Start with the row containing '$'
    
    for (int i = 0; i < n; i++) {
        result += table[idx].first;
        idx = table[idx].second;
    }
    
    // Rotate to put '$' at the end
    return result.substr(1) + result[0];
}

/**
 * FM-Index for efficient pattern matching using BWT.
 */
class FMIndex {
private:
    string text;
    string bwt;
    map<char, int> C;  // Number of characters lexicographically smaller
    map<char, vector<int>> occ;  // Occurrence array
    
    void build_index() {
        int n = bwt.length();
        
        // Count occurrences of each character
        map<char, int> counts;
        for (char c : bwt) {
            counts[c]++;
        }
        
        // Compute C array
        int total = 0;
        for (auto& p : counts) {
            C[p.first] = total;
            total += p.second;
        }
        
        // Build occurrence array (Occ)
        for (auto& p : counts) {
            occ[p.first] = vector<int>(n + 1, 0);
        }
        
        for (int i = 0; i < n; i++) {
            for (auto& p : occ) {
                p.second[i + 1] = p.second[i];
            }
            occ[bwt[i]][i + 1]++;
        }
    }
    
    pair<int, int> backward_search(const string& pattern) {
        int top = 0;
        int bottom = bwt.length() - 1;
        
        // Process pattern from right to left
        for (int i = pattern.length() - 1; i >= 0; i--) {
            char c = pattern[i];
            
            if (C.find(c) == C.end()) {
                return {1, 0};  // Pattern not found
            }
            
            // Update range using LF mapping
            top = C[c] + occ[c][top];
            bottom = C[c] + occ[c][bottom + 1] - 1;
            
            if (top > bottom) {
                return {1, 0};  // Pattern not found
            }
        }
        
        return {top, bottom};
    }

public:
    FMIndex(string input_text) {
        // Add sentinel if not present
        if (input_text.empty() || input_text.back() != '$') {
            input_text += '$';
        }
        
        text = input_text;
        bwt = burrows_wheeler_transform(input_text);
        build_index();
    }
    
    /**
     * Count occurrences of pattern in the text.
     */
    int count(const string& pattern) {
        auto [top, bottom] = backward_search(pattern);
        if (top > bottom) {
            return 0;
        }
        return bottom - top + 1;
    }
    
    /**
     * Find all positions where pattern occurs in text.
     */
    vector<int> locate(const string& pattern) {
        auto [top, bottom] = backward_search(pattern);
        
        if (top > bottom) {
            return {};
        }
        
        // Create suffix array for position lookup
        int n = text.length();
        vector<pair<string, int>> rotations;
        for (int i = 0; i < n; i++) {
            rotations.push_back({text.substr(i) + text.substr(0, i), i});
        }
        sort(rotations.begin(), rotations.end());
        
        vector<int> positions;
        for (int i = top; i <= bottom; i++) {
            positions.push_back(rotations[i].second);
        }
        
        sort(positions.begin(), positions.end());
        return positions;
    }
};

int main() {
    // Example usage
    string text = "ACGTACGTACGT";
    
    cout << "Burrows-Wheeler Transform and FM-Index" << endl;
    cout << "Original text: " << text << endl;
    
    // BWT
    string bwt = burrows_wheeler_transform(text);
    cout << "BWT: " << bwt << endl;
    
    // Inverse BWT
    string original = inverse_bwt(bwt);
    cout << "Inverse BWT: " << original << endl;
    cout << "Reconstruction correct: " << (original == text + "$" ? "true" : "false") << endl;
    
    // FM-Index search
    cout << "\n--- FM-Index Search ---" << endl;
    FMIndex fm_index(text);
    
    vector<string> patterns = {"ACG", "CGT", "TAC", "XYZ"};
    for (const auto& pattern : patterns) {
        int count = fm_index.count(pattern);
        cout << "\nPattern: " << pattern << endl;
        cout << "Count: " << count << endl;
        
        if (count > 0) {
            vector<int> positions = fm_index.locate(pattern);
            cout << "Positions: ";
            for (size_t i = 0; i < positions.size(); i++) {
                cout << positions[i];
                if (i < positions.size() - 1) cout << ", ";
            }
            cout << endl;
        }
    }
    
    return 0;
}
