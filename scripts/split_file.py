#!/usr/bin/env python3
import argparse
from math import ceil

def parse_args():
    parser = argparse.ArgumentParser(
        description="Split file into batches")
    parser.add_argument("-i",
                        "--input",
                        type=str,
                        required=True,
                        help="Input file")
    parser.add_argument("-o",
                        "--output-files",
                        nargs='+',
                        type=str,
                        required=True,
                        help="Output file(s) list")
    parser.add_argument("-l",
                        "--lines-per-record",
                        type=int,
                        default=2,
                        help="Number of lines per record in input")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    print(args)
    assert(args.lines_per_record >= 1)
    lines = open(args.input, 'r').readlines()
    assert(len(lines)%args.lines_per_record==0)
    record_cnt  = len(lines)/args.lines_per_record
    outer_count = len(args.output_files)
    records_per_outer  = ceil(record_cnt / outer_count)
    line_num = 0
    try:
        for outer in args.output_files:
            outer = open(outer, 'w+')
            for _ in range(records_per_outer):
                for _ in range(args.lines_per_record):
                    print(lines[line_num], end='', file=outer)
                    line_num+=1
    except:
        assert(line_num==len(lines))

if __name__ == "__main__":
    main()
