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
        
        input = open(args.infile, "r")
        first = input.readline()
        if first.startswith("ERROR:"):
            while True:
                try:
                    writer.write(first)
                except BrokenPipeError:
                    break
        else:
            for line in header:
                writer.write(line)
            writer.write(first)
            for line in input:
                writer.write(line)
        sys.stderr.write("Finished parsing input\n")
        input.close()
