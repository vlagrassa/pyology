# Built-in libraries
import argparse
import numpy    as np
import sys

# Local files
from ._constants   import BLOSUM62_r, PAM250_r, SIMPLE_r, protein_codes
from ._show_matrix import *


# Constant bit filters for moving left, up, and down
# Since any given cell can theoretically move in all three directions, traceback
# values are represented as a mask of the appropriate filters
# 
# E.g. If the diagonal and up scores are identical and better than the left score,
#      that cell will get the value 0b010 | 0b100 = 0b110 = 6
# 
MOVE_L = 0b001
MOVE_U = 0b010
MOVE_D = 0b100

# Unicode "arrow" values associated with each possible combination of moves
ARROWS = [ '.', '\u2190', '\u2191', '\u2518', '\u2196', '\u29a3', '\u29ad', '\u05e6']



## Parsing functions for string arguments ##

# Map type string arguments to boolean representing whether to use local algorithm
def parse_type(type_string):
    if   type_string == "global":      return False
    elif type_string == "local":       return True

# Map matrix names to the actual dictionaries
def parse_matrix(matrix_name):
    if   matrix_name == "blosum62":    return BLOSUM62_r
    elif matrix_name == "pam250":      return PAM250_r
    elif matrix_name == "simple":      return SIMPLE_r



# Define default argument strings and values
# Defined up here to standardize them across the command-line version and the import-able "align" function

DEF_TYPE_STR   = "global"
DEF_MATRIX_STR = "simple"

DEFAULT_TYPE    = parse_type(DEF_TYPE_STR)
DEFAULT_MATRIX  = parse_matrix(DEF_MATRIX_STR)
DEFAULT_PENALTY = -1
    
    
# Create the arguments to parse
parser = argparse.ArgumentParser()
parser.add_argument(
    '--seq1', metavar="FILENAME", required=True,
    help="FASTA formatted file containing first sequence. If input starts with a \"*\" character, interpreted as a raw sequence instead."
)
parser.add_argument(
    '--seq2', metavar="FILENAME", required=True,
    help="FASTA formatted file containing second sequence. If input starts with a \"*\" character, interpreted as a raw sequence instead."
)
parser.add_argument(
    '--type',
    choices=['global', 'local'], default=DEF_TYPE_STR,
    help="Specify global or local match (default " + DEF_TYPE_STR + ")."
)
parser.add_argument(
    '--matrix',
    choices=['simple', 'blosum62', 'pam250'], default=DEF_MATRIX_STR,
    help="Name of scoring matrix (default \"" + DEF_MATRIX_STR + "\": +1 for match, -1 for mismatch)."
)
parser.add_argument(
    '--gap_penalty', metavar="PENALTY",
    type=int, default=DEFAULT_PENALTY,
    help="Penalty for introducing a gap in either sequence (default " + str(DEFAULT_PENALTY) + ")."
)
parser.add_argument(
    '--one_result', action="store_true",
    help="If set, return the first max-scoring path found (by default finds all)."
)
parser.add_argument(
    '--show_tables', action="store_true",
    help="If set, display the scoring and traceback matrices (hidden by default)."
)



# Find the indices in the matrix with the maximum score, to start the traceback
def find_max_idxs(matrix):
    result = np.where(matrix == np.amax(matrix))
    return list(zip(result[0], result[1]))



# Generator function that yields valid paths through the matrix
# Paths are not returned in any particular order, since they all have the same score
# Paths are already "reversed" - the sequence of steps moves down and rightward
def traceback(score_matrix, arrow_matrix, coords):
    
    # Unpack current location and value in the matrix
    row, col = coords
    curr_val = int(arrow_matrix[row, col])
    
    # Base case
    # On reaching a 0 value, yield the current coordinates as the origin of the alignment
    # and an empty list to start recursively building up the paths from
    if curr_val == 0:
        yield ( coords, [] )
    
    # Recursive step
    # Recursively call traceback from each next possible cell in the path
    else:
    
        # Pair each movement option with the starting coords if we were to take that movement
        opts = [
            (move, (row + d[0], col + d[1])) for (move, d) in traceback.OPTIONS if curr_val & move
        ]
        
        # Get the corresponding scores for each possible path, along with the maximum score
        scores = [
            score_matrix[new_coords] for (_, new_coords) in opts
        ]
        max_score = max(scores)

        # For each move that would go to a max-square cell...
        for (move, new_coords), score in zip(opts, scores):
            if score == max_score:
                
                # ...recursively call traceback starting from that cell, appending the move that
                # got us there to the end of the path
                for ( origin, partial_path ) in traceback( score_matrix, arrow_matrix, new_coords ):
                    yield ( origin, partial_path + [move] )



# Static OPTIONS variable for the traceback function
# Pairs each possible move with the coordinate offset associated with that move
traceback.OPTIONS = [
    ( MOVE_D, (-1, -1) ),
    ( MOVE_L, ( 0, -1) ),
    ( MOVE_U, (-1,  0) ),
]


# Generator function that produces pairs of characters from the given sequences
def path_to_alignment(source, path, seq1, seq2):

    # Start the (i,j) indices from the source location
    i, j = source

    # Loop through each step in the path, updating both sequences appropriately
    for step in path:

        # Upward step
        # Consume a character from the first (vertical) sequence only
        if step == MOVE_U:
            yield (seq1[i], '-')
            i += 1

        # Left step
        # Consume a character from the second (horizontal) sequence only
        if step == MOVE_L:
            yield ('-', seq2[j])
            j += 1

        # Diagonal step
        # Consume a character from each sequence
        if step == MOVE_D:
            yield (seq1[i], seq2[j])
            i += 1
            j += 1



# The main logic of the program
def align(
    seq1, seq2,
    is_local      = DEFAULT_TYPE,
    score         = DEFAULT_MATRIX,
    gap_penalty   = DEFAULT_PENALTY,
    one_result    = False,
    return_tables = False,
    whole_seqs    = False,
    no_mismatch   = False,
):
    
    # The arrow matrix uses bit-masked values to determine where to go
    # If the first  bit is 1, moving diagonal produces an optimal score
    # If the second bit is 1, moving upward   produces an optimal score
    # If the third  bit is 1, moving leftward produces an optimal score
    # If all bits are 0, we are at the start of the sequence
    #  - 000: the source (top left)
    #  - 001: move left
    #  - 010: move up
    #  - 100: move diagonal

    from ._sequence import Seq

    # Ensure sequences are Seq objects
    if not isinstance(seq1, Seq):
        seq1 = Seq(seq1)
    if not isinstance(seq2, Seq):
        seq2 = Seq(seq2)

    # Create a scoring matrix from the "score" parameter
    if isinstance(score, str):
        matrix = parse_matrix(score)
    elif isinstance(score, tuple):
        keys = "CSTPAGNDEQHRKMILVFYW"
        matrix = { ord(key1): { ord(key2): score[0] if key1 == key2 else score[1] for key2 in keys } for key1 in keys }
    else:
        matrix = score

    # Compute the number of columns and rows
    rows = len(seq1) + 1    # Sequence 1 is vertical
    cols = len(seq2) + 1    # Sequence 2 is horizontal
    
    # Create the result matrix and the arrow (directions) matrix
    # I store the directions as chars rather than ints - it's a bit more space efficient,
    # since I only need three bits for each entry
    score_matrix = np.zeros((rows, cols))
    arrow_matrix = np.full( (rows, cols), '0' )
    
    # Initialize the first row and column using numpy
    # If using local alignment, we want the first row/col to be filled with zeroes
    #    - In the score matrix, this represents a potential start point
    #    - In the arrow matrix, this represents a location to terminate the traceback
    if not is_local:
        
        # The scores just compound the gap penalty
        score_matrix[:, 0] = np.arange(rows) * gap_penalty
        score_matrix[0, :] = np.arange(cols) * gap_penalty
        
        # The arrows just point left and upward
        # Explicitly slice from 1 to keep the top left corner as a terminating 0
        arrow_matrix[1:, 0 ] = np.full(rows - 1, MOVE_U)
        arrow_matrix[0 , 1:] = np.full(cols - 1, MOVE_L)

    # Loop through each cell in the matrices, skipping the first row/col
    for i in range(1, rows):
        for j in range(1, cols):
            
            # Compute the score from the left and upper cells (add gap penalty)
            score_l = score_matrix[i  , j-1] + gap_penalty
            score_u = score_matrix[i-1, j  ] + gap_penalty

            if no_mismatch:
                score_d = min(score_l, score_u) - 1

            else:
                # Compute the score from the diagonal cell (add dictionary matrix value)
                # score_d = score_matrix[i-1, j-1] + matrix[seq1[i-1]][seq2[j-1]]
                base_1 = seq1.get_raw(i-1)
                base_2 = seq2.get_raw(j-1)
                if base_1 == ord('-') or base_2 == ord('-'):
                    add_score = 0
                else:
                    add_score = matrix[base_1][base_2]

                score_d = score_matrix[i-1, j-1] + add_score
            
            # Compute the maximum possible score and store it in the result array
            score_matrix[i, j] = max_score = max(score_l, score_u, score_d)
            
            # Special case: in a local search, all would-be negative values map to 0
            if is_local and max_score < 0:
                score_matrix[i, j] = 0
                arrow_matrix[i, j] = 0
            
            # For each direction, if that direction would give the max score,
            # add it to the result with bit masking
            else:
                direction = 0
                if score_l == max_score:    direction |= MOVE_L
                if score_u == max_score:    direction |= MOVE_U
                if score_d == max_score:    direction |= MOVE_D
                arrow_matrix[i, j] = direction
            
    
    # Get the start indices
    #   - for a local search, start from the max value(s)
    #   - for a global search, use the bottom-right corner
    if is_local:    start_idxs = find_max_idxs(score_matrix)
    else:           start_idxs = [(rows - 1, cols - 1)]
    
    # Get the total alignment score (will be the same for all possible start idxs)
    align_score = score_matrix[start_idxs[0]]
    
    # Initialize empty list to store alignments
    alignments = []

    # Loop through all possible starting indices
    # Will only be multiple indices if local search finds multiple cells w max value
    for start_idx in start_idxs:
        for (source, path) in traceback(score_matrix, arrow_matrix, start_idx):

            # If whole sequences are desired, add everything before the matching region
            if whole_seqs:
                s1 = seq1[0 : source[0]] >> source[1]
                s2 = seq2[0 : source[1]] >> source[0]

            # If not, initialize empty Sequences
            else:
                s1 = Seq()
                s2 = Seq()

            # Convert sequences to an alignment, appending each new pair of chars
            for (c1, c2) in path_to_alignment(source, path, seq1, seq2):
                s1 += c1
                s2 += c2

            # Add the rest of the sequences
            if whole_seqs:
                s1 += seq1[start_idx[0] :] << (len(seq2) - start_idx[1])
                s2 += seq2[start_idx[1] :] << (len(seq1) - start_idx[0])

            # Add the alignment to the list
            alignments += [(s1, s2)]
            
            # If the user just wants the first result, break out of both loops
            if one_result:
                alignments = alignments[0]
                break
        
        # Break out of the outer loop if the inner is broken
        else: continue
        break

    # Return the score and alignment(s), as well as the matrices if requested
    if return_tables:    return align_score, alignments, score_matrix, arrow_matrix
    else:                return align_score, alignments



if __name__ == "__main__":

    # Parse the command line arguments
    a = parser.parse_args()

    # Parse the first sequence
    if a.seq1[0] == '*':    seq1 = a.seq1[1:]
    else:                   seq1 = next(read_fasta.read_fasta_sequences(a.seq1))['sequence']

    # Parse the second sequence
    if a.seq2[0] == '*':    seq2 = a.seq2[1:]
    else:                   seq2 = next(read_fasta.read_fasta_sequences(a.seq2))['sequence']

    # Run the main align function
    score, alignments, score_matrix, arrow_matrix = align(
        seq1, seq2, parse_type(a.type), parse_matrix(a.matrix), a.gap_penalty, a.one_result, return_tables=True
    )
    
    print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=")
    
    # Print the alignments
    print("")
    for i, (s1, s2) in enumerate(alignments):
        print(f"  {i:3}:  {''.join(s1)}")
        print(f"        {''.join(s2)}")
        print("")

    # Print the overall alignment score
    print(f"  Alignment Score: {score}")
    print("")
    
    print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=")
    
    # Print the raw tables, if desired
    if a.show_tables:
    
        row_labels = [label for label in "-"+seq1]
        col_labels = [label for label in "-"+seq2]

        print("")
        print("Scoring array:")
        show_matrix.display_table(score_matrix, row_labels, col_labels)
        print("Traceback array:")
        show_matrix.display_table([ [ARROWS[int(x)] for x in row] for row in arrow_matrix ], row_labels, col_labels)
        print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=")
