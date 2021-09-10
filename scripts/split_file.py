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
    parser.add_argument("-op",
                        "--out-files-pattern",
                        type=str,
                        required=True,
                        help="Output file(s) list")
    parser.add_argument("-b",
                        "--batch-count",
                        type=int,
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
    record_cnt  = len(lines)//args.lines_per_record


    old_bid = -1
    for idx,line in enumerate(lines):
        bid = int(args.batch_count*((idx//args.lines_per_record)/(len(lines)//args.lines_per_record)))
        if old_bid == -1:
            old_bid = bid
            outer = open(args.out_files_pattern.format(bid), 'w+')
            bid_cnt = 0
        if old_bid != bid:
            assert bid_cnt%args.lines_per_record==0
            outer = open(args.out_files_pattern.format(bid), 'w+')
            old_bid = bid
            bid_cnt = 0
        bid_cnt += 1
        outer.write(lines[idx])
    outer.close()
    assert bid_cnt%args.lines_per_record==0

if __name__ == "__main__":
    main()
