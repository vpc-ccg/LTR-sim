import sys
import numpy as np
import copy as cp
from random import choices
from collections import Counter
from math import ceil

def write_new_variation(header, seq, fout, ls=60):
	fout.write(header)
	fout.write('\n')

	i = 0
	while i < len(seq):
		if i + ls < len(seq):
			fout.write(seq[i:i+ls])
		else:
			fout.write(seq[i:])
		fout.write('\n')

		i += ls


# makes a list of size part containing integer values that add up to num
# the difference between the numbers in the list is not greater than one
# the list is sorted in reverse order
def split_number(num, part):
	small_n = int(num / part)
	great_n = small_n + 1

	great_cnt = num - small_n * part
	small_cnt = part - great_cnt

	final = [ great_n ] * great_cnt + [ small_n ] * small_cnt

	return final

def sample_with_replacement(seq_len, slope, intercept, depth, points=1000):
	if seq_len > 10000:
		slope //= 1
	elif 75000 > seq_len > 5000:
		slope //= 1.5
	elif 5000 > seq_len > 3000:
		slope //= 3
	elif 3000 > seq_len > 1500:
		slope //= 6
	elif 1500 > seq_len > 0:
		slope //= 10
	all_lens = []
	cur_len = seq_len
	for i in range(points):
		if i % 2 == 0:
			percentage = cur_len/seq_len
			print('{:6.2f}% {:5}nt: {}'.format(percentage*100, cur_len, '*'*ceil(50*percentage)))
		# print('-'*200*int())
		all_lens.append(cur_len)
		cur_len = max(intercept, cur_len - slope)

	lens = choices(all_lens, k=depth)
	return lens

def degrade(header, seq, slope, intercept, fout):
	header_list = header.strip().split()
	depth = None
	for field in header_list:
		if 'depth=' in field:
			depth = int(field.split('=')[1])
			break
	if depth == None:
		raise Exception('Contig with header <{}> does not have depth comment'.format(header))

	seq_len = len(seq)

	lens = sample_with_replacement(seq_len, slope, intercept, depth)

	i = 0
	len_counts = Counter(lens)
	for length in sorted(len_counts):
		dep = len_counts[length]
		new_depth = 'depth={}'.format(dep)

		new_header_list = cp.deepcopy(header_list)
		new_header_list[0] += '_D{}'.format(i)
		new_header_list[1] = '{} offset={}'.format(new_depth, -1*(seq_len - length))
		new_header = ' '.join(new_header_list)

		new_seq = seq[ : length]
		write_new_variation(new_header, new_seq, fout)
		i += 1


def process_all(trans_in, trans_out, slope, intercept):
	fout = open(trans_out, 'w')

	i = 0
	header = ''
	seq = ''

	with open(trans_in) as f:
		for l in f:
			if l[0] == '>':
				if seq != '':
					degrade(header, seq, slope, intercept, fout)
					seq = ''
					i += 1
				header = l
			else:
				seq += l.strip()

	degrade(header, seq, slope, intercept, fout)

def usage():
	print('\nUsage: python {} transcriptome_FASTA output_FASTA slope intercept'.format(sys.argv[0]))

def main():
	print('This script generates degraded variations of the given transcripts')
	print(sys.argv)
	args = sys.argv[1:]

	if len(args) != 4:
		usage()
		exit(1)

	try:
		trans_in = args[0]
		trans_out = args[1]
		slope = int(args[2])
		intercept = int(args[3])
	except ValueError as err:
		print('Error in arguments!\n{}'.format(err))
		exit(1)

	print('Input FASTA: {}\nOutput FASTA: {}\nLinear function: y = {} * x + {}'.format(trans_in, trans_out, slope, intercept))

	process_all(trans_in, trans_out, slope, intercept)

if __name__ == '__main__':
	main()
