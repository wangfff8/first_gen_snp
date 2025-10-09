import argparse
import sys


def parse_fasta(filepath: str) -> dict[str, str]:
    """
    Parses a FASTA file and returns a dictionary of sequence headers and sequences.
    """
    seqs_dict: dict[str, str] = {}
    current_header: str | None = None
    current_sequence: list[str] = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_header:
                        seqs_dict[current_header] = "".join(current_sequence).upper()
                    current_header = line[1:]
                    current_sequence = []
                else:
                    current_sequence.append(line)
            if current_header:
                seqs_dict[current_header] = "".join(current_sequence).upper()
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: Input file '{filepath}' not found.")
    except Exception as e:
        raise IOError(f"Error reading file '{filepath}': {e}")
    return seqs_dict

def colorize_multiple_alignment_html(seqs_dict: dict[str, str], line_length: int = 80) -> str:
    """
    Colorizes a multiple sequence alignment based on column conservation and returns HTML.
    - Green for conserved positions (all non-gap characters are identical).
    - Red for non-conserved positions (non-gap characters are different).
    - White for gap positions.
    All characters are converted to uppercase.
    The output is formatted with fixed line_length characters per line, with sequence IDs on the left,
    and sequences are outputted in an interleaved manner.
    """
    if not seqs_dict:
        return ""

    headers = list(seqs_dict.keys())
    sequences = list(seqs_dict.values())

    alignment_length = len(sequences[0])
    if not all(len(seq) == alignment_length for seq in sequences):
        raise ValueError("All sequences in the alignment must have the same length.")

    num_sequences = len(sequences)

    html_output_lines = []
    
    # Add CSS for preformatted text and colors
    html_output_lines.append("<!DOCTYPE html>")
    html_output_lines.append("<html>")
    html_output_lines.append("<head>")
    html_output_lines.append("<meta charset=\"UTF-8\">")
    html_output_lines.append("<title>Colored MAFFT Alignment</title>")
    html_output_lines.append("<style>")
    html_output_lines.append("body { background-color: #333; color: #eee; }") # Dark background for better contrast
    html_output_lines.append("pre { font-family: 'Courier New', Courier, monospace; white-space: pre; margin: 0; padding: 10px; }")
    html_output_lines.append(".green { color: #00ff00; }") # Brighter green
    html_output_lines.append(".red { color: #ff0000; }")   # Brighter red
    html_output_lines.append(".gap { color: #ffffff; background-color: #555; }") # Gaps in white on dark grey for visibility
    html_output_lines.append(".header { color: #00ffff; font-weight: bold; }") # Cyan for headers
    html_output_lines.append("</style>")
    html_output_lines.append("</head>")
    html_output_lines.append("<body>")
    html_output_lines.append("<pre>")

    # Store colored characters for each sequence
    colored_sequences_chars: list[list[str]] = [[] for _ in range(num_sequences)]

    for i in range(alignment_length): # Iterate through each column
        column_chars = [seq[i] for seq in sequences]
        
        is_conserved = True
        first_non_gap_char = None
        has_non_gap = False

        for char in column_chars:
            if char != '-':
                has_non_gap = True
                if first_non_gap_char is None:
                    first_non_gap_char = char
                elif char != first_non_gap_char:
                    is_conserved = False
                    break
        
        for seq_idx in range(num_sequences): # Apply color to each character in the current column
            char = sequences[seq_idx][i]
            if char == '-':
                colored_sequences_chars[seq_idx].append(f'<span class="gap">{char}</span>')
            elif is_conserved and has_non_gap:
                colored_sequences_chars[seq_idx].append(f'<span class="green">{char}</span>')
            else:
                colored_sequences_chars[seq_idx].append(f'<span class="red">{char}</span>')
    
    # Format and reconstruct colored sequences with interleaved output
    # Calculate max_raw_header_len for alignment
    max_raw_header_len = max(len(h) for h in headers) if headers else 0
    header_padding = 2 # Number of spaces between header and sequence segment
    
    # Calculate the total visual width for the header part (including padding)
    total_header_visual_width = max_raw_header_len + header_padding

    for i in range(0, alignment_length, line_length):
        for seq_idx, header in enumerate(headers):
            # Pad the raw header text to the total desired visual width
            display_header_text = header.ljust(total_header_visual_width)
            
            # Slice the list of colored character spans directly
            segment_spans = colored_sequences_chars[seq_idx][i : i + line_length]
            segment_html = "".join(segment_spans)
            
            # Append the formatted line: header (left-justified with padding inside span) + sequence segment
            html_output_lines.append(f'<span class="header">{display_header_text}</span>{segment_html}')
        html_output_lines.append("") # Add an empty line for better separation between blocks

    html_output_lines.append("</pre>")
    html_output_lines.append("</body>")
    html_output_lines.append("</html>")
    
    return "\n".join(html_output_lines)

def main():
    parser = argparse.ArgumentParser(
        description="Colorizes a multiple sequence alignment (e.g., MAFFT output) based on column conservation and outputs an HTML file."
    )
    parser.add_argument(
        "input_file",
        help="Path to the input FASTA-formatted multiple sequence alignment file."
    )
    parser.add_argument(
        "-o", "--output",
        default="colored_alignment.html",
        help="Name of the output HTML file. Default is 'colored_alignment.html'."
    )
    parser.add_argument(
        "-l", "--line-length",
        type=int,
        default=80,
        help="Number of characters per line in the output. Default is 80."
    )
    args = parser.parse_args()

    try:
        seqs_dict = parse_fasta(args.input_file)
        html_output = colorize_multiple_alignment_html(seqs_dict, args.line_length)

        with open(args.output, "w", encoding="utf-8") as f:
            f.write(html_output)
        print(f"Colored alignment saved to '{args.output}'.")
    except FileNotFoundError as e:
        print(e, file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Data error: {e}", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(e, file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
