"""
Burrows-Wheeler Transform (BWT) and FM-Index

The BWT is a reversible transformation that reorganizes a string to make it more compressible.
The FM-index uses the BWT to enable fast exact pattern matching and is foundational to
BWA, Bowtie, and HISAT2 aligners.
"""


def burrows_wheeler_transform(text):
    """
    Compute the Burrows-Wheeler Transform of a text.
    
    Args:
        text: Input string (must end with unique character like '$')
    
    Returns:
        str: BWT of the input text
    """
    # Add sentinel character if not present
    if not text.endswith('$'):
        text = text + '$'
    
    # Create all rotations
    rotations = [text[i:] + text[:i] for i in range(len(text))]
    
    # Sort rotations lexicographically
    rotations.sort()
    
    # BWT is the last column of sorted rotations
    bwt = ''.join(rotation[-1] for rotation in rotations)
    
    return bwt


def inverse_bwt(bwt):
    """
    Reverse the Burrows-Wheeler Transform.
    
    Args:
        bwt: BWT string
    
    Returns:
        str: Original text
    """
    n = len(bwt)
    
    # Create table with indices
    table = [(bwt[i], i) for i in range(n)]
    
    # Sort to get first column
    table.sort()
    
    # Follow the links to reconstruct
    result = []
    idx = 0  # Start with the row containing '$'
    
    for _ in range(n):
        result.append(table[idx][0])
        idx = table[idx][1]
    
    # Remove sentinel and return
    text = ''.join(result)
    return text[1:] + text[0]  # Rotate to put '$' at the end


class FMIndex:
    """
    FM-Index for efficient pattern matching using BWT.
    
    The FM-index enables backward search for exact pattern matching in O(m) time
    where m is the pattern length, independent of the text length.
    """
    
    def __init__(self, text):
        """
        Build FM-index from text.
        
        Args:
            text: Input text (string)
        """
        # Add sentinel if not present
        if not text.endswith('$'):
            text = text + '$'
        
        self.text = text
        self.bwt = burrows_wheeler_transform(text)
        self._build_index()
    
    def _build_index(self):
        """Build auxiliary data structures for FM-index."""
        n = len(self.bwt)
        
        # Count occurrences of each character
        self.counts = {}
        for char in self.bwt:
            self.counts[char] = self.counts.get(char, 0) + 1
        
        # Compute C array - number of characters lexicographically smaller
        self.C = {}
        sorted_chars = sorted(self.counts.keys())
        total = 0
        for char in sorted_chars:
            self.C[char] = total
            total += self.counts[char]
        
        # Build occurrence array (Occ)
        # Occ[char][i] = number of occurrences of char in bwt[0:i]
        self.occ = {char: [0] * (n + 1) for char in self.counts}
        
        for i, char in enumerate(self.bwt):
            for c in self.occ:
                self.occ[c][i + 1] = self.occ[c][i]
            self.occ[char][i + 1] += 1
    
    def count(self, pattern):
        """
        Count occurrences of pattern in the text.
        
        Args:
            pattern: Pattern to search for
        
        Returns:
            int: Number of occurrences
        """
        top, bottom = self._backward_search(pattern)
        if top > bottom:
            return 0
        return bottom - top + 1
    
    def _backward_search(self, pattern):
        """
        Perform backward search using FM-index.
        
        Args:
            pattern: Pattern to search for
        
        Returns:
            tuple: (top, bottom) range in suffix array
        """
        top = 0
        bottom = len(self.bwt) - 1
        
        # Process pattern from right to left
        for i in range(len(pattern) - 1, -1, -1):
            char = pattern[i]
            
            if char not in self.C:
                return (1, 0)  # Pattern not found
            
            # Update range using LF mapping
            top = self.C[char] + self.occ[char][top]
            bottom = self.C[char] + self.occ[char][bottom + 1] - 1
            
            if top > bottom:
                return (1, 0)  # Pattern not found
        
        return (top, bottom)
    
    def locate(self, pattern):
        """
        Find all positions where pattern occurs in text.
        
        Args:
            pattern: Pattern to search for
        
        Returns:
            list: List of positions where pattern occurs
        """
        top, bottom = self._backward_search(pattern)
        
        if top > bottom:
            return []
        
        # For simplicity, we reconstruct positions
        # In practice, this would use sampled suffix array
        positions = []
        
        # Create suffix array
        rotations = [(self.text[i:] + self.text[:i], i) for i in range(len(self.text))]
        rotations.sort()
        suffix_array = [rot[1] for rot in rotations]
        
        for i in range(top, bottom + 1):
            positions.append(suffix_array[i])
        
        return sorted(positions)
    
    def search(self, pattern):
        """
        Search for exact matches of pattern.
        
        Args:
            pattern: Pattern to search for
        
        Returns:
            dict: Search results with count and positions
        """
        count = self.count(pattern)
        positions = self.locate(pattern) if count > 0 else []
        
        return {
            'pattern': pattern,
            'count': count,
            'positions': positions
        }


if __name__ == "__main__":
    # Example usage
    text = "ACGTACGTACGT"
    
    print("Burrows-Wheeler Transform and FM-Index")
    print(f"Original text: {text}")
    
    # BWT
    bwt = burrows_wheeler_transform(text)
    print(f"BWT: {bwt}")
    
    # Inverse BWT
    original = inverse_bwt(bwt)
    print(f"Inverse BWT: {original}")
    print(f"Reconstruction correct: {original == text + '$'}")
    
    # FM-Index search
    print("\n--- FM-Index Search ---")
    fm_index = FMIndex(text)
    
    patterns = ["ACG", "CGT", "TAC", "XYZ"]
    for pattern in patterns:
        result = fm_index.search(pattern)
        print(f"\nPattern: {pattern}")
        print(f"Count: {result['count']}")
        if result['positions']:
            print(f"Positions: {result['positions']}")
