import os

def byte_file_lines(filename, len_lines, num_lines, reverse=False, offset=0):
    """
    Read a file as bytes in forward or backward order.
    
    Assumes all lines (except potentially the last) are the same length.
    """

    # Open the file in read mode & byte mode
    file = open(filename, "r+b")

    # Generate the lines in forward order, where offset defines the start
    if not reverse:
        file.seek(offset)
        for line in file:
            yield line.strip()

    # Generate the lines in reverse order, where offset defines the end
    else:

        # Make sure to count newlines
        line_offset = len(b'\n')

        # Start the index at the beginning of the final line
        line_start_idx = (len_lines + line_offset) * (num_lines - 1)

        # Seek to offset bits before the end, and get the absolute position
        file.seek(-offset, os.SEEK_END)
        offset_abs_idx = file.tell()

        # Back up the index to the start of the line containing the start offset
        while offset_abs_idx < line_start_idx:
            line_start_idx -= len_lines + line_offset

        # If the offset is in the middle of the line, return the first part of the line
        if line_start_idx < offset_abs_idx:

            # Seek to the start of the line containing the offset and read it
            file.seek(line_start_idx)
            next_line = file.readline().strip()

            # Cut out everything after the offset
            yield next_line[0: offset_abs_idx - line_start_idx]

            # Move back to start of previous line
            line_start_idx -= len_lines + line_offset

        # Move backwards line by line until no lines are left
        while line_start_idx >= 0:
            file.seek(line_start_idx)

            # Generate the current line
            yield file.readline().strip()

            # Move back to start of previous line
            line_start_idx -= len_lines + line_offset



def read_byte_file(filename, len_lines, num_lines, window_size=None, reverse=False, offset=0):
    """
    Read a file in sets of window_size bytes, ignoring line breaks.

    Assumes all lines (except potentially the last) are the same length.
    """

    # If no size specified, yield line-by-line
    if window_size is None:
        for line in byte_file_lines(filename, len_lines, num_lines, reverse=reverse, offset=offset):
            if reverse:    yield line[::-1]
            else:          yield line

    # Otherwise, yield in chunks of given size
    else:

        # Initialize buffer to store unyielded chars from previous line(s)
        buffer = []

        # Loop through each line starting from the offset in the file
        for line in byte_file_lines(filename, len_lines, num_lines, reverse=reverse, offset=offset):

            # Reverse the line, if applicable
            if reverse:
                line = line[::-1]

            # Move into buffer and yield until line runs out
            while len(line) > 0:

                # Compute how many more characters we can fit in the buffer
                num_to_add = window_size - len(buffer)

                # Move that many characters from line to buffer
                buffer += line[: num_to_add]
                line    = line[num_to_add :]

                # If buffer is big enough, yield and clear it
                if len(buffer) >= window_size:
                    yield buffer
                    buffer = []

        # If anything is left in the buffer, yield it
        if len(buffer) > 0:
            yield buffer
