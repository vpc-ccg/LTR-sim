import sys

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
    cdna_fasta_file = argv[1]
    gtf_file = argv[2]
    fusion_info_file = argv[3]
    output_file = argv[4]
    fusion_index_file = argv[5]
    fusion_modes = {
        "WILDTYPES",
        "DISTRIBUTE"
    }
    fuse_strategy = "WILDTYPES"

    fusions = []
    fusion_genes = {}
    with open(fusion_info_file, 'r') as hand:
        for line in hand:
            if line[0] == "#":
                continue
            fields = line.rstrip().split("\t")

            g1 = fields[8]
            g2 = fields[9]
            
            chr1 = fields[0]
            chr2 = fields[3]

            pos1 = int(fields[1])
            pos2 = int(fields[4])
            genotype = fields[7]
            sv_type = fields[6]
            
            if len(fields) > 10:
                t1= fields[10]
                t2 = fields[11]
                
                fusions.append( (g1,g2,chr1,chr2,pos1,pos2,genotype,sv_type,True,t1,t2))
#22 24806290 24806291 22 26467150 26467151 inversion 1/0 ENSG00000167037 ENSG00000100099
            else:
                fusions.append( (g1,g2,chr1,chr2,pos1,pos2,genotype,sv_type,False))
            fusion_genes[g1] = fusions[-1] 
            fusion_genes[g2] = fusions[-1]

#1       havana  transcript      12010   13670   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000450305"; transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-201"; transcript_source "havana"; transcript_biotype "transcribed_unprocessed_pseudogene"; tag "basic"; transcript_support_level "NA";
    gene_to_strand = {}
    gene_to_transcript = {}
    transcript_to_exons = {}
    gene_id_to_name = {}
    with open(gtf_file,'r') as hand:
        for line in hand:
            if line[0] == "#":
                continue
            fields = line.rstrip().split("\t")

            info = fields[8]
            info_fields =  [x.strip() for x in fields[8].split(";")[:-1]]
            info_dic = {x.split()[0]:x.split()[1][1:-1] for x in info_fields}
            #print(info_dic)
            gene = info_dic["gene_id"]

            if fields[2] == "transcript":
                transcript = info_dic["transcript_id"]
                if gene not in gene_to_strand:
                    gene_to_strand[gene] = fields[6]
                if gene not in gene_to_transcript:
                    gene_to_transcript[gene] = []
                gene_to_transcript[gene].append(transcript)
            elif fields[2] == "exon":
                transcript = info_dic["transcript_id"]
                if transcript not in transcript_to_exons:
                    transcript_to_exons[transcript] = []
                transcript_to_exons[transcript].append((fields[0],int(fields[3]),int(fields[4])))
            elif fields[2] == "gene":
                gene_id_to_name[gene] = info_dic["gene_name"]

            #transcripts[]
    """
>>> st[1:16]
'ENSG00000223997'
>>> st[18:36]
'_ENST00000415118.1'
>>> st[19:36]
'ENST00000415118.1'
>>> st[19:34]
'ENST00000415118
    """
#>ENST00000525964_D0 depth=8 offset=0 cdna chromosome:GRCh38:11:134152847:134224216:-1 gene:ENSG00000151503.12 gene_biotype:protein_coding transcript_biotype:nonsense_mediated_decay gene_symbol:NCAPD3 description:non-SMC condensin II complex subunit D3 [Source:HGNC Symbol;Acc:HGNC:28952]
    fused_seqs = {}
    with open(cdna_fasta_file, 'r') as hand, open(output_file, 'w') as whand, open(fusion_index_file, 'w') as ihand:
        line = hand.readline()
        while line:
            header = line.rstrip()
            seq = hand.readline()
            seq = seq.rstrip()
#            gene=header[1:16]
#            transcript=header[19:34]
            hfil = header.split(" ")
            transcript = hfil[0][1:16]
            depth_str = hfil[1][6:]
            depth = int(depth_str)
            gene = hfil[4][5:20]
            if gene in fusion_genes:
                #print(header.rstrip(),"to_fuse",sep="\t",file=whand)
                #print(seq.rstrip(),file=whand)
                fused_seqs[transcript] = (header,seq,depth,gene)
                if fusion_genes[gene][6] != "1/1": #Genotype
                    print(hfil[0],"depth={}".format(int(depth/2)),sep=" ",end=" ",file=whand)
                    print(" ".join(hfil[2:]),file=whand)
                    print(seq,file=whand)
            else:
                print(header,file=whand)
                print(seq,file=whand)
            line = hand.readline()
            #if transcript not in fusion_genes:
        fusion_transcript_index = 0
        if fuse_strategy == "DISTRIBUTE":
            for f in fusions:
                pass
        elif fuse_strategy == "WILDTYPES":
            for fusion_index,f in enumerate(fusions):
                g1 = f[0]
                g2 = f[1]
                
                t1l = gene_to_transcript[g1]
                t2l = gene_to_transcript[g2]



                s1 = gene_to_strand[g1]
                s2 = gene_to_strand[g2]
                if f[8]:
                    max_depth_transcript1 = f[9]
                    seqtup = fused_seqs[f[9]]

                    sum_dep1 = 0
                    for t in t1l:
                        seqtup = fused_seqs[t]
                        depth = seqtup[2]
                        sum_dep1+=depth
                    max_depth1 = sum_dep1
                else:
                    max_depth_transcript1 = "-1"
                    max_depth1 = -1
                    for t in t1l:
                        seqtup = fused_seqs[t]
                        depth = seqtup[2]
                        if depth > max_depth1:
                            max_depth1 = depth
                            max_depth_transcript1 = t

                if f[8]:
                    max_depth_transcript2 = f[10]
                    seqtup = fused_seqs[f[10]]
                    max_depth2 = seqtup[2]

                else:

                    max_depth_transcript2 = "-1"
                    max_depth2 = -1
                    for t in t2l:
                        seqtup = fused_seqs[t]
                        depth = seqtup[2]
                        if depth > max_depth2:
                            max_depth2 = depth
                            max_depth_transcript2 = t
                    

                exons1 = sorted(transcript_to_exons[max_depth_transcript1],key=lambda x: x[1]) #Sort by start pos of exons since - strand exons will be reversely sorted
                exons2 = sorted(transcript_to_exons[max_depth_transcript2],key=lambda x: x[1])

                breakpoint1 = f[4]
                breakpoint2 = f[5]
                

                fusion_sequence = ""

                bases1 = 0
                bases2 = 0
                if s1 == s2: #deletion cases
                    bases1=get_base_count(exons1,breakpoint1,s1 == "+")
                    fusion_sequence += fused_seqs[max_depth_transcript1][1][:bases1]
                    
                    bases2=get_base_count(exons2,breakpoint2,s1 != "+")
                    fusion_sequence += fused_seqs[max_depth_transcript2][1][-bases2:]

                else: #inversion cases
                    if s1 == "+":

                        bases1=get_base_count(exons1,breakpoint1,True)
                        fusion_sequence += fused_seqs[max_depth_transcript1][1][:bases1]

                        bases2=get_base_count(exons2,breakpoint2,True)
                        fusion_sequence += fused_seqs[max_depth_transcript2][1][-bases2:]


                    else:
                        bases1=get_base_count(exons1,breakpoint1,False)
                        fusion_sequence += fused_seqs[max_depth_transcript2][1][:bases2]

                        bases2=get_base_count(exons2,breakpoint2,False)
                        fusion_sequence += fused_seqs[max_depth_transcript1][1][-bases1:]



                fusion_gene_id = "FUSG{0:011}".format(fusion_index)
                fusion_transcript_id = "FUST{0:011}".format(fusion_transcript_index)
                fusion_name = "{}::{}".format(gene_id_to_name[g1], gene_id_to_name[g2])
                print(fusion_name,bases1,bases2,sep="\t",file=sys.stderr)

                header = ">{0} depth={1} cdna {10} {11} source_gene1={2} source_gene2={3} breakpoints={4}-{5}:{6}-{7} source_transcript1={8} source_transcript2={9}" \
                            .format(fusion_transcript_id, max_depth1, g1, g2, chr1, chr2, pos1, pos2, \
                                max_depth_transcript1,max_depth_transcript2, fusion_gene_id, fusion_name)
                print(header, file=whand)
                print(fusion_sequence, file=whand)
                
                print(fusion_gene_id,fusion_transcript_id,fusion_name,sep="\t",file=ihand)
                fusion_transcript_index+=1
        else:
            pass

if __name__ == "__main__":
    exit(main(len(sys.argv),sys.argv))
