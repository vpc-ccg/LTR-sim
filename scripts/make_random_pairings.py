import sys
import random

def get_base_count(exons, bp, forward):

    bases = 0
    if forward:
        for e in exons:
            if e[2] > bp:
                if e[1] < bp:
                    bases = bases + bp - e[1]
                break
            else:
                bases = bases + e[2]-e[1]
    else:
        for e in reversed(exons):
            if e[1] < bp:
                if e[2] > bp:
                    bases = bases + e[2] - bp
                break
            else:
                bases = bases + e[2]-e[1]
                   
    return bases

def main(argc, argv):
    fastq_file = argv[1]
    
    pairing_ratio = float(argv[2])

    output_file = argv[3]

    rp_index_file = argv[4]

    rp_index = 0
    
    rp_seqs  = []
    rp_quals = [] 
    rp_flag = False
    rp_headers = []

    with open(fastq_file, 'r') as hand, open(output_file, 'w') as whand, open(rp_index_file, 'w') as ihand:
        line = hand.readline()
        while line:
            header = line.rstrip()
            seq = hand.readline().rstrip()
            plus = hand.readline().rstrip()
            qual = hand.readline().rstrip()
            rand_val = random.random()
            if rand_val <= pairing_ratio: #FUSE
                rp_flag = True
                rp_headers.append(header)
                rp_seqs.append(seq)
                rp_quals.append(qual)
            else: #SKIP and print if a rp exists 
                if rp_flag: #pritn rp read
                    
                    rp_headers.append(header)
                    rp_seqs.append(seq)
                    rp_quals.append(qual)
                    new_id = "RP{0:011}".format(rp_index)
                    rp_index+=1
                    print("@{}".format(new_id), file=whand)
                    for s in rp_seqs:
                        print(s, end="", file=whand)
                    print("\n+", file=whand)
                    for q in rp_quals:
                        print(q, end="", file=whand)
                    print(file=whand)

                    #PRINT INDEX FOR FUSED READS                    
                    print("{}\t{}".format(new_id, len(rp_headers)),file=ihand)
                    for rhead, rseq in zip(rp_headers,rp_seqs):
                        print( "{}\t{}".format(len(rseq),rhead), file=ihand)
                    #RESET LISTS

                    rp_seqs  = []
                    rp_quals = [] 
                    rp_flag = False
                    rp_headers = []

                else: #print non-rp read
                    print(header,file=whand)
                    print(seq,file=whand)
                    print(plus,file=whand)
                    print(qual,file=whand)

            line = hand.readline()

if __name__ == "__main__":
    exit(main(len(sys.argv),sys.argv))
