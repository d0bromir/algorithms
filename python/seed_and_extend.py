"""
Seed-and-Extend Algorithm using K-mer Hashing

This algorithm is used for fast sequence alignment by:
1. Finding exact k-mer matches (seeds) using hash tables
2. Extending seeds to find longer alignments

This is the foundational approach used in BLAST, MAQ, and SOAP.
"""

from collections import defaultdict


class KmerIndex:
    """Index for storing k-mer positions in a reference sequence."""
    
    def __init__(self, sequence, k=11):
        """
        Build k-mer index from a reference sequence.
        
        Args:
            sequence: Reference sequence (string)
            k: K-mer length (default: 11)
        """
        self.sequence = sequence
        self.k = k
        self.index = defaultdict(list)
        self._build_index()
    
    def _build_index(self):
        """Build hash table index of k-mer positions."""
        for i in range(len(self.sequence) - self.k + 1):
            kmer = self.sequence[i:i + self.k]
            self.index[kmer].append(i)
    
    def find_seeds(self, query):
        """
        Find exact k-mer matches between query and reference.
        
        Args:
            query: Query sequence (string)
        
        Returns:
            list: List of (query_pos, ref_pos) tuples for each seed match
        """
        seeds = []
        for i in range(len(query) - self.k + 1):
            kmer = query[i:i + self.k]
            if kmer in self.index:
                for ref_pos in self.index[kmer]:
                    seeds.append((i, ref_pos))
        return seeds


def extend_seed(seq1, seq2, seed_pos1, seed_pos2, match_score=1, mismatch_penalty=-1):
    """
    Extend a seed match in both directions.
    
    Args:
        seq1: First sequence (string)
        seq2: Second sequence (string)
        seed_pos1: Starting position in seq1
        seed_pos2: Starting position in seq2
        match_score: Score for matching characters
        mismatch_penalty: Penalty for mismatching characters
    
    Returns:
        tuple: (start1, end1, start2, end2, score)
    """
    # Extend right
    score = 0
    max_score = 0
    max_right = 0
    
    i = 0
    while (seed_pos1 + i < len(seq1) and 
           seed_pos2 + i < len(seq2)):
        if seq1[seed_pos1 + i] == seq2[seed_pos2 + i]:
            score += match_score
        else:
            score += mismatch_penalty
        
        if score > max_score:
            max_score = score
            max_right = i + 1
        
        # Stop if score drops too much below max
        if score < max_score - 5:
            break
        
        i += 1
    
    # Extend left
    score = max_score
    max_left = 0
    
    i = 1
    while (seed_pos1 - i >= 0 and 
           seed_pos2 - i >= 0):
        if seq1[seed_pos1 - i] == seq2[seed_pos2 - i]:
            score += match_score
        else:
            score += mismatch_penalty
        
        if score > max_score:
            max_score = score
            max_left = i
        
        # Stop if score drops too much below max
        if score < max_score - 5:
            break
        
        i += 1
    
    start1 = seed_pos1 - max_left
    end1 = seed_pos1 + max_right
    start2 = seed_pos2 - max_left
    end2 = seed_pos2 + max_right
    
    return start1, end1, start2, end2, max_score


def seed_and_extend(reference, query, k=11):
    """
    Perform seed-and-extend alignment between query and reference.
    
    Args:
        reference: Reference sequence (string)
        query: Query sequence (string)
        k: K-mer length (default: 11)
    
    Returns:
        list: List of alignment hits with (query_start, query_end, ref_start, ref_end, score)
    """
    # Build k-mer index
    kmer_index = KmerIndex(reference, k)
    
    # Find seeds
    seeds = kmer_index.find_seeds(query)
    
    # Extend each seed
    alignments = []
    for query_pos, ref_pos in seeds:
        start1, end1, start2, end2, score = extend_seed(
            query, reference, query_pos, ref_pos
        )
        alignments.append({
            'query_start': start1,
            'query_end': end1,
            'ref_start': start2,
            'ref_end': end2,
            'score': score,
            'query_seq': query[start1:end1],
            'ref_seq': reference[start2:end2]
        })
    
    # Sort by score descending
    alignments.sort(key=lambda x: x['score'], reverse=True)
    
    return alignments


if __name__ == "__main__":
    # Example usage
    reference = "ACGTACGTACGTAAACCCGGGTTTACGTACGT"
    query = "ACGTAAACCCGGG"
    
    print("Seed-and-Extend Alignment")
    print(f"Reference: {reference}")
    print(f"Query: {query}")
    print(f"\nK-mer length: 5")
    
    alignments = seed_and_extend(reference, query, k=5)
    
    print(f"\nFound {len(alignments)} alignment(s):")
    for i, aln in enumerate(alignments[:5], 1):  # Show top 5
        print(f"\nAlignment {i}:")
        print(f"  Query position: {aln['query_start']}-{aln['query_end']}")
        print(f"  Reference position: {aln['ref_start']}-{aln['ref_end']}")
        print(f"  Score: {aln['score']}")
        print(f"  Query sequence: {aln['query_seq']}")
        print(f"  Ref sequence:   {aln['ref_seq']}")
