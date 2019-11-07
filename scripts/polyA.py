import sys
import numpy as np
import scipy.stats as stats

def count_transcripts(trans_in):
	tr_cnt = 0
	with open(trans_in) as f:
		for l in f:
			if l[0] == '>':
				tr_cnt += 1

	print('Contig count: {}'.format(tr_cnt))
	return tr_cnt

def generate_polyA_len_list(cnt, mean, std, mincut, maxcut):
	a, b = (mincut - mean) / std, (maxcut - mean) / std

	lens = stats.truncnorm.rvs(a, b, loc=mean, scale=std, size=cnt, random_state=42)
	lens = np.rint(lens).astype(int)

	print('Mean: {}'.format(np.mean(lens)))
	print('STD: {}'.format(np.std(lens)))
	print('Min: {}'.format(np.min(lens)))
	print('Max: {}'.format(np.max(lens)))

	return lens

def add_tail(seq, lsize, tail_len, line_size=100000000000):
	if line_size - lsize > tail_len:
		new_seq = seq.strip() + 'A' * tail_len + '\n'
	else:
		new_seq = seq.strip() + 'A' * (line_size - lsize) + '\n'
		remain_len = tail_len - (line_size - lsize)
		full_lines_cnt = int(remain_len / line_size)
		new_seq += ('A' * line_size + '\n') * full_lines_cnt
		remain_len -= line_size * full_lines_cnt
		if remain_len > 0:
			new_seq += 'A' * remain_len + '\n'

	return new_seq

def add_tail_all(trans_in, trans_out, lens):
	fout = open(trans_out, 'w')

	i = 0
	header = ''
	seq = ''
	lsize = 0

	with open(trans_in) as f:
		for l in f:
			if l[0] == '>':
				if seq != '':
					#print(header.strip())
					#print(lens[i])
					new_seq = add_tail(seq, lsize, lens[i])
					fout.write(header)
					fout.write(new_seq)
					seq = ''
					i += 1
				header = l
			else:
				seq += l
				lsize = len(l.strip())

	#print(header.strip())
	#print(lens[i])
	new_seq = add_tail(seq, lsize, lens[i])
	fout.write(header)
	fout.write(new_seq)

def usage():
	print('\nUsage: python {} transcriptome_FASTA output_FASTA mean std mincut maxcut'.format(sys.argv[0]))

def main():
	print('This script adds poly-A sequences to the tail of transcripts')
	args = sys.argv[1:]

	if len(args) != 7:
		usage()
		exit(1)

	try:
		trans_in = args[0]
		trans_out = args[1]
		mean = int(args[2])
		std = int(args[3])
		mincut = int(args[4])
		maxcut = int(args[5])
		seed = int(args[6])
	except ValueError as err:
		print('Error in arguments!\n{}'.format(err))
		exit(1)
	np.random.seed(seed)
	print('Input FASTA: {}\nOutput FASTA: {}\n(Expected) Mean: {}, std: {}, mincut: {}, maxcut: {}'.format(trans_in, trans_out, mean, std, mincut, maxcut))

	tr_cnt = count_transcripts(trans_in)
	polyA_lens = generate_polyA_len_list(tr_cnt, mean, std, mincut, maxcut)

	add_tail_all(trans_in, trans_out, polyA_lens)

if __name__ == '__main__':
	main()
