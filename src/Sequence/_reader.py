import os

from ._sequence     import Seq
from ._constants    import protein_codes, frame_values
from ._file_reading import read_byte_file
from ._transform    import invert_gen


# The default window size for a reader
DEFAULT_WINDOW_SIZE = 100



class SeqCombiner:

    def __init__(self, op, *seqs):
        self.op   = op
        self.seqs = seqs

    def scan(self, window_size=DEFAULT_WINDOW_SIZE):

        # Map the sequences to generators with the same window size
        mapped_seqs = map( lambda x: x.as_DNA(window_size=window_size), self.seqs )

        # Zip the the next output from each sequence together and apply the operation
        return (
            self.op(*lines) for lines in zip(*mapped_seqs)
        )

    def __iter__(self):
        return self.scan()


# class SeqGen:




class SeqReader:


    ### Initialization and Deletion ###

    def __init__(self, name, desc, filename, num_lines, len_lines, window_size=DEFAULT_WINDOW_SIZE):

        # Set sequence info
        self.name = name
        self.desc = desc

        # Store info about the tempfile where the sequence is stored
        self.filename  = filename
        self.num_lines = num_lines
        self.len_lines = len_lines

        # Set the default window size
        self.window_size = window_size

        # self.is_local = False
        # self.sequence = None

        # Compute the length of the final line
        with open(self.filename, "r+b") as f:
            f.seek((self.len_lines + len(b'\n')) * (self.num_lines - 1))
            last_line = f.readline()
            self.last_line_len = len(last_line)


    def __del__(self):
        # print(f"Deleting file {self.filename}...")
        if os.path.exists(self.filename):
            os.remove(self.filename)



    ### Imported Methods ###

    from ._orf import find_orfs



    ### String Representation ###

    def __repr__(self):
        return f"SeqReader({self.name})"



    ### FASTA Header ###

    def get_ID(self):
        return self.name

    def get_desc(self):
        return self.desc

    def get_header_line(self, newline=False):
        if len(self.desc) > 0:
            return f">{self.name} {self.desc}" + ("\n" if newline else "")
        else:
            return f">{self.name}"  + ("\n" if newline else "")



    ### Sequence Length ###

    def __len__(self):
        return ((self.num_lines - 1) * self.len_lines) + self.last_line_len



    ### Transcription and Translation ###

    def as_DNA(self, complement=False, offset=0, window_size=None):

        # Get the default window size, if none is provided
        if window_size is None:
            window_size = self.window_size

        # Read chunks from the temp file and map to Seq objects
        generator = map( Seq, read_byte_file(
            self.filename, self.len_lines, self.num_lines,
            window_size=window_size, reverse=complement,
            offset=offset,
        ))

        # Invert the bases coming out of the generator, if applicable
        if complement:    yield from invert_gen(generator)
        else:             yield from            generator


    def as_RNA(self, complement=False, offset=0, window_size=None):
        yield from Seq.dna_to_rna(
            self.as_DNA(complement, offset, window_size),
            as_generator=True
        )

    def translate(self, frame=1):

        # Only accept valid frame values
        assert frame in frame_values

        # Offset is one closer to 0 than the frame identifier
        offset = frame - (1 if frame > 0 else -1)

        yield from Seq.translate(
            self.as_RNA(
                complement  = frame < 0,
                offset      = abs(offset),
                window_size = 3,
            ),
            as_generator=True
        )

    def __and__(self, other):
        return SeqCombiner( Seq.__and__, self, other )

    def __rand__(self, other):
        return SeqCombiner( Seq.__and__, other, self )

    def __or__(self, other):
        return SeqCombiner( Seq.__or__,  self, other )

    def __ror__(self, other):
        return SeqCombiner( Seq.__or__, other, self )

    def __xor__(self, other):
        return SeqCombiner( Seq.__xor__, self, other )

    def __rxor__(self, other):
        return SeqCombiner( Seq.__xor__, other, self )


    # @staticmethod
    # def local_reader(sequence, window_size=None):
    #     reader = SeqReader("local_sequence", "", None, 0, window_size=window_size)
    #     reader.set_local_sequence(sequence)
    #     return reader



def parse_id_and_desc(line):
    id_and_desc = line.split(' ', 1)

    seq_id   = id_and_desc[0].strip()[1:]
    seq_desc = id_and_desc[1].strip() if len(id_and_desc) > 1 else ""

    return seq_id, seq_desc



def read_first(filename, **kwargs):
    try:       return next(read_all(filename, **kwargs))
    except:    raise Exception(f"No such sequence found in file {filename}.")


# Produce SeqReaders for all matching sequences in the file
def read_all(filename, name=None, id_filter=None, size=None):

    from tempfile import NamedTemporaryFile

    # Default value for match function is equality
    if id_filter is None:
        id_filter = lambda a, b: b is None or a == b

    # Initialize vars to track most recent Reader and positions in file
    tempfile  = None
    num_lines = 0
    len_lines = 80

    line_buffer = ""

    def reset_loop_vars():
        nonlocal tempfile, num_lines, line_buffer

        tempfile    = None
        num_lines   = 0
        line_buffer = ""


    def write_buffer(line):
        nonlocal tempfile, num_lines, len_lines, line_buffer

        if tempfile is not None:

            line_buffer += line.strip()

            while len(line_buffer) > len_lines:

                tempfile.write(bytearray(line_buffer[0:len_lines], encoding="utf-8"))
                tempfile.write(b"\n")

                line_buffer = line_buffer[len_lines:]
                num_lines += 1


    def empty_buffer():
        nonlocal tempfile, num_lines, len_lines, line_buffer

        if len(line_buffer) > 0:
            tempfile.write(bytearray(line_buffer, encoding="utf-8"))
            num_lines += 1



    # Read the file
    with open(filename, "r") as f:

        # Ensure first line is a valid FASTA identifier line
        line = f.readline()
        if line[0] != '>':
            raise Exception("Invalid FASTA file: First line must begin with '>'.")

        # Loop through the lines, skipping blank space
        while line:
            if len(line.strip()) == 0: continue

            # Start of new sequence
            if line[0] == '>':

                # If we were building a previous sequence, finish it and reset
                if tempfile is not None:
                    empty_buffer()
                    # tempfile.seek(0)
                    yield SeqReader(seq_id, seq_desc, tempfile.name, num_lines, len_lines)
                    reset_loop_vars()

                # Parse out the sequence ID and description
                seq_id, seq_desc = parse_id_and_desc(line)

                # If ID matches filter, start building a new sequence
                if id_filter(seq_id, name):
                    tempfile = NamedTemporaryFile(delete=False)

            # Content line of existing sequence
            else:
                write_buffer(line)

            # Move to next line
            line = f.readline()

    # Yield the final sequence we were building
    if tempfile is not None:
        empty_buffer()
        # tempfile.seek(0)
        yield SeqReader(seq_id, seq_desc, tempfile.name, num_lines, len_lines)
