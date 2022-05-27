import os
import time
from functools       import partial
from multiprocessing import Pool

from ._constants import frame_values



def find_orfs(self, frame=None, min_length=3, processes=6, verbose=False, verbose_output=False):
    """
    Find all open reading frames in a given sequence.

    If a reading frame is provided, searches that frame only.

    If not, searches all six reading frames in parallel.
    """

    # Search all frames
    if frame is None:

        # Start timer
        start_time = time.time()

        # Initialize empty list to track ORFs
        orfs = []

        # Initialize var to sum time of each computation
        sum_time = 0

        # Protect the temp file from garbage collection (???)
        # Honestly I have no idea why this happens but without this the whole thing breaks
        f = open(self.filename)

        # Spawn processes to read frames concurrently
        with Pool(processes=processes) as pool:

            # Loop through frame values, assigning a new process to read each one
            for frame, orf_list, t in pool.imap_unordered(
                partial(find_orfs, self, min_length=min_length, verbose_output=True),
                frame_values
            ):
                if verbose:
                    print(f"  - Read frame {frame:2} in {t} seconds.")
                orfs.extend(orf_list)
                sum_time += t

        # Safe to close the temporary handle now (for some reason)
        f.close()

        # Stop timer
        end_time = time.time()

        if verbose:
            print("")
            print(f"Total computing time: ~{sum_time} seconds.")
            print(f"Actual elapsed time:   {end_time - start_time} seconds.")

        # Return the built up list of ORFs
        return orfs


    # Search a specific frame
    else:

        # Start timer
        start_time = time.time()

        # Create generator for amino acids, reading from specified frame
        seq = self.translate(frame)

        # Initialize loop variables
        orfs   = []
        starts = []

        # Loop through amino acids one at a time
        for idx, char in enumerate(seq):

            # On start marker, add to list of current ORFs
            if char == 'M':
                starts.append(idx)

            # On stop marker, close all open ORFs and clear the list
            elif char == '*':
                orfs  += [ (frame, s, idx) for s in starts if idx - s >= min_length ]
                starts = []

        # Stop timer
        end_time = time.time()

        # Return the frame (for convenience), all found orfs, and the runtime
        if verbose_output:
            return frame, orfs, end_time - start_time
        else:
            return orfs
