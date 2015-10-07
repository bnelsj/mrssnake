import sys
import argparse
import pysam

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("outfile")
    parser.add_argument("template_file")

    args = parser.parse_args()

    header = pysam.view("-H", args.template_file)
    with open(args.outfile, "w") as writer:

        for line in header:
            writer.write(line)
        input = open(args.infile, "r")
        first = input.readline()
        if first.startswith("Chunker:"):
            while True:
                try:
                    writer.write(first)
                except BrokenPipeError:
                    break
        else:
            sys.stderr.write("Parser: got first line\n")
            sys.stderr.flush()
            writer.write(first)
            for line in input:
                if not line.startswith("Chunker:"):
                    writer.write(line)
            input.close()
        sys.stderr.write("Parser: finished parsing input\n")
        sys.stderr.flush()
    sys.exit()
