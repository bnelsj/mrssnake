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
        for line in input:
            writer.write(line)
        input.close()
