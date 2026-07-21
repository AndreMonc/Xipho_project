header = None
length = 0
with open('xiph_elegans_ref.fa') as fasta:
    for line in fasta:
        # Trim newline
        line = line.rstrip()
        if line.startswith('>'):
            # If we captured one before, print it now
            if header is not None:
                print(header, length)
                length = 0
            header = line[1:]
        else:
            length += len(line)
# Don't forget the last one
if length:
    print(header, length)
