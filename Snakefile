configfile: 'config.yaml'

import os
import sys

sys.setrecursionlimit(100)
def get_abs_path(path):
    abs_path = os.popen('readlink -f {}'.format(path)).read()
    return abs_path.rstrip("\r\n")

def make_slurm():
    os.makedirs('slurm'.format(outpath), mode=0o777, exist_ok=True)

outpath = get_abs_path(config['outpath'])
make_slurm()
config['exec']['ltrsim'] = get_abs_path(config['exec']['ltrsim'])

map_d     = '{}/map'.format(outpath)
train_d   = '{}/train'.format(outpath)
reads_d   = '{}/reads'.format(outpath)
batches_d = '{}/batches'.format(reads_d)

localrules:
    all,
    batched_get_throughput
    
rule all:
    input:
        # expand('{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}.fasta'.format(reads_d),
        #        sample=config['samples'],
        #        degradation_level=config['degradation_level'],),
        # expand('{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}.batch-{{bid}}.fasta'.format(batches_d),
        #        sample=config['samples'],
        #        degradation_level=config['degradation_level'],
        #        bid=list(range(config['badread']['batches']))),
        # expand('{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}.batch-{{bid}}.reads.fastq'.format(batches_d),
        #        sample=config['samples'],
        #        degradation_level=config['degradation_level'],
        #        bid=list(range(config['badread']['batches']))),
        f1=expand('{}/{{sample}}.L-{{degradation_level}}.{{extension}}'.format(reads_d),
               sample=config['samples'],
               degradation_level=config['degradation_level'],
               extension=['fastq','tsv']),
        frp=expand('{}/{{sample}}.R-{{degradation_level}}.fastq'.format(reads_d),
               sample=config['samples'],
               degradation_level=config['degradation_level']),

rule git_badread:
    output:
        'extern/Badread/badread-runner.py'
    shell:
        'git clone --branch {tag} {url} extern/Badread'.format(tag=config['badread']['tag'], url=config['badread']['url'])

rule reads_to_bam:
    conda:
        'conda.env'
    input:
        reads  = lambda wildcards: config['samples'][wildcards.sample],
        target = lambda wildcards: config['references'][wildcards.target]
    output:
        paf = '{}/{{sample}}.{{target}}.paf'.format(map_d),
    params:
        mapping_settings = lambda wildcards: config['mapping_settings'][wildcards.target]
    conda:
        'conda.env'
    threads:
        32
    shell:
        'minimap2 -Y -x {params.mapping_settings} --eqx --MD -t {threads} {input.target} {input.reads} > {output.paf} '

rule expression_rate:
    input:
        paf = '{}/{{sample}}.cdna.paf'.format(map_d)
    output:
        exp = '{}/{{sample}}.expression_rate.tsv'.format(train_d)
    run:
        print('Capturing expression rate from: {}'.format(input.paf))
        tid_to_rcnt = dict()
        for idx,line in enumerate(open(input.paf)):
            if idx % 100000 == 0:
                if idx < 1000000:
                    print('Read {:4d}K alns'.format(idx//1000))
                else:
                    print('Read {:4.1f}M alns'.format(idx/1000000))
            if not 'tp:A:P' in line:
                continue
            line = line.rstrip().split('\t')
            # rid = l[0]
            tid = line[5].split('_')[-1].split('.')[0]
            tid_to_rcnt[tid] = tid_to_rcnt.get(tid, 0) + 1
        outfile = open(output.exp, 'w+')
        for tid,rcnt in tid_to_rcnt.items():
            print('{}\t{}'.format(tid, rcnt), file=outfile)
        outfile.close()

rule sim_transcriptome:
    input:
        cdna = config['references']['cdna'],
        exp  = '{}/{{sample}}.expression_rate.tsv'.format(train_d)
    output:
        cdna = '{}/{{sample}}.cdna.fasta'.format(reads_d)
    run:
        print('Simulating transcriptome from {}'.format(input))
        tid_to_rcnt = dict()
        genes = set(config['genes'])
        chroms = set([str(c) for c in config['chroms']])
        print(chroms)
        for line in open(input.exp):
            line = line.rstrip().split('\t')
            tid = line[0]
            cnt = int(line[1])
            tid_to_rcnt[tid] = cnt
        outfile = open(output.cdna, 'w+')
        seq=list()
        for line in open(input.cdna):
            line = line.rstrip()
            if line[0] == '>':
                if len(seq) > 0:
                    print(''.join(seq), file=outfile)
                    seq = list()
                gene_name_field_name = 'gene:'
                gene_name_field_sepa = '.'
                gene_name = line[line.find(gene_name_field_name):]
                gene_name = gene_name[len(gene_name_field_name) : gene_name.find(gene_name_field_sepa)]
                chrom_name_field_name = 'chromosome:GRCh38:'
                chrom_name_field_sepa = ':'
                chrom_name = line[line.find(chrom_name_field_name) + len(chrom_name_field_name):].split(chrom_name_field_sepa)[0]
                if not chrom_name in chroms or not ('All_genes' in genes or gene_name in genes):
                    flag = False
                    continue
                flag = True
                line = line[1:]
                line = line.split()
                tid = line[0].split('.')[0]
                cnt = tid_to_rcnt.get(tid, 0) + 1
                record = [tid, 'depth={}'.format(max(cnt,0))] + line[1:]
                print('>{}'.format(' '.join(record)), file=outfile)
            elif flag:
                seq.append(line)
        if len(seq) > 0:
            print(''.join(seq), file=outfile)
        outfile.close()

rule prep_rt_tsv:
    input:
        make_rt=config["exec"]["makert"],
        rt_file=config["fusionsim"]["rtfile"],
        gtf    = config['annotations']['gtf'],
    output:
        tsv  = '{}/{{sample}}.rt.tsv'.format(reads_d),
    shell:
        "PYTHONHASHSEED=0 python {input.make_rt} {input.rt_file} {input.gtf} {output.tsv}"
rule merge_fusion_and_rt_tsv:
    input:
        fus = config['fusionsim']['fusion_info'],
        rt = '{}/{{sample}}.rt.tsv'.format(reads_d),
    output:
        '{}/{{sample}}.fusion_and_rt.tsv'.format(reads_d),
    shell:
        "cat {input.fus} {input.rt} > {output}"
rule fuse_genes:
    input:
        make_fusions = config['exec']['fuser'],
        cdna  = '{}/{{sample}}.cdna.fasta'.format(reads_d),
        gtf    = config['annotations']['gtf'],
        fusion_info='{}/{{sample}}.fusion_and_rt.tsv'.format(reads_d),
    #    fusion_info    = config['fusionsim']['fusion_info'],
    output:
        cdna  = '{}/{{sample}}.cdna.fused.fasta'.format(reads_d),
        fus_index  = '{}/{{sample}}.cdna.fused.tsv'.format(reads_d),
    params:
        valid_chrs = ",".join([str(x) for x in config["chroms"]])
    shell:
        'PYTHONHASHSEED=0 python3 {input.make_fusions} {input.cdna}'
        ' {input.gtf} {input.fusion_info} {output.cdna} {output.fus_index} {params.valid_chrs}'
rule add_poly_A:
    conda:
        'conda.env'
    input:
        polyA = config['exec']['polyA'],
        cdna  = '{}/{{sample}}.cdna.fused.fasta'.format(reads_d),
    output:
        cdna_A = '{}/{{sample}}.cdna.polyA.fasta'.format(reads_d)
    params:
        mean   = config['polyA_config']['mean'],
        stddev = config['polyA_config']['stddev'],
        mincut = config['polyA_config']['mincut'],
        maxcut = config['polyA_config']['maxcut'],
        seed   = config['seed'],
    shell:
        'PYTHONHASHSEED=0 python {input.polyA} '
        '   {input.cdna} {output.cdna_A}'
        '   {params.mean} {params.stddev}'
        '   {params.mincut} {params.maxcut}'
        '   {params.seed}'


rule degrade:
    conda:
        'conda.env'
    input:
        degrade = config['exec']['degrade'],
        cdna_A = '{}/{{sample}}.cdna.polyA.fasta'.format(reads_d)
    output:
        cdna_A_degraded = '{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}.fasta'.format(reads_d)
    wildcard_constraints:
        sample='|'.join(config['samples']),
        degradation_level='|'.join(config['degradation_level'])
    params:
        intercept = lambda wildcards: config['degradation_level'][wildcards.degradation_level]['intercept'],
        slope     = lambda wildcards: config['degradation_level'][wildcards.degradation_level]['slope'],
        seed      = config['seed'],
    shell:
        'PYTHONHASHSEED=0 python {input.degrade}'
        '   {input.cdna_A} {output.cdna_A_degraded}'
        '   {params.slope} {params.intercept}'
        '   {params.seed}'

rule batch_fasta:
    input:
        split_file      = config['exec']['split_file'],
        cdna_A_degraded = '{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}.fasta'.format(reads_d)
    output:
        ['{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}.batch-{}.fasta'.format(batches_d,x) for x in range(int(config["badread"]["batches"]))]
    wildcard_constraints:
        sample='|'.join(config['samples']),
        degradation_level='|'.join(config['degradation_level'])
    params:
        lines_per_fasta = 2
    shell:
        '{input.split_file} -i {input.cdna_A_degraded} -o {output} -l {params.lines_per_fasta}'

rule batched_get_throughput:
    input:
        cdna_A_degraded = '{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}.batch-{{bid}}.fasta'.format(batches_d)
    output:
        cdna_A_degraded_cov = '{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}.batch-{{bid}}.cov.txt'.format(batches_d)
    wildcard_constraints:
        sample='|'.join(config['samples']),
        degradation_level='|'.join(config['degradation_level']),
        bid='|'.join([str(b) for b in range(config['badread']['batches'])])
    run:
        amplified_throughput = 0.0
        throughput = 0.0
        for line_num,line in enumerate(open(input.cdna_A_degraded)):
            line = line.rstrip()
            if line_num % 2 == 0:
                header_list = line.split()
            	for idx,field in enumerate(header_list):
                    if 'depth=' in field:
                        depth = int(field.split('=')[1])
            elif line_num % 2 == 1:
                amplified_throughput+=len(line)*depth
                throughput+=len(line)
                print(len(line),depth,amplified_throughput,throughput)
        print('{:.2f}'.format(amplified_throughput/throughput), end='', file=open(output.cdna_A_degraded_cov, 'w+'))

rule batched_generate_reads:
    conda:
        'conda.env'
    input:
        badread             = config['exec']['badread'],
        cdna_A_degraded     = temp('{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}.batch-{{bid}}.fasta'.format(batches_d)),
        cdna_A_degraded_cov = '{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}.batch-{{bid}}.cov.txt'.format(batches_d)
    output:
        fastq               = '{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}.batch-{{bid}}.reads.fastq'.format(batches_d),
    wildcard_constraints:
        sample='|'.join(config['samples']),
        degradation_level='|'.join(config['degradation_level']),
        bid='|'.join([str(b) for b in range(config['badread']['batches'])])
    params:
        # coverage = config['badread']['coverage'],
        seed     = config['seed'],
        other_params = " --force_strand "
    shell:
        'PYTHONHASHSEED=0 {input.badread} simulate {params.other_params} --seed {params.seed} --length 100000,0'
        '   --reference {input.cdna_A_degraded} --quantity $(cat {input.cdna_A_degraded_cov})x'
        '   --chimeras 0 --random_reads 0 --junk_reads 0 --glitches 0,0,0'
        '   > {output.fastq}'

rule randomly_pair:
    input:
        fastq  = '{}/{{sample}}.S-{{degradation_level}}.fastq'.format(reads_d),
        script = config["exec"]["randomly_pair"],
    output:
        fastq = '{}/{{sample}}.R-{{degradation_level}}.fastq'.format(reads_d),
        index = '{}/{{sample}}.R-{{degradation_level}}.fastq.rpindex'.format(reads_d),
    params:
        pair_ratio=config["fusionsim"]["random_pair_ratio"],
    shell:
        "python {input.script} {input.fastq} {params.pair_ratio} {output.fastq} {output.index}"
    
rule shuffle_reads:
    input:
       fastq = '{}/{{sample}}.L-{{degradation_level}}.fastq'.format(reads_d) 
    output:
       fastq = '{}/{{sample}}.S-{{degradation_level}}.fastq'.format(reads_d)
    threads:
        15
    shell:
        "cat {input.fastq} | seqkit shuffle -j15 > {output.fastq}"

rule reads_info:
    input:
        fastqs = ['{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}.batch-{}.reads.fastq'.format(batches_d,x) for x in range(int(config["badread"]["batches"]))],
        gtf    = config['annotations']['gtf'],
        fus_index  = '{}/{{sample}}.cdna.fused.tsv'.format(reads_d),
    output:
        tsv   = '{}/{{sample}}.L-{{degradation_level}}.tsv'.format(reads_d),
        fastq = '{}/{{sample}}.L-{{degradation_level}}.fastq'.format(reads_d),
    wildcard_constraints:
        sample='|'.join(config['samples']),
        degradation_level='|'.join(config['degradation_level']),
    run:
        print('Parsing read information from {}'.format(input))
        t_headers = [
            'gene_id',
            'transcript_id',
            'gene_name',
        ]
        tid_to_info = dict(
            non={f:'NA' for f in t_headers}
        )

        for line in open(input.fus_index):
            line = line.rstrip()
            fields = line.split('\t')
            tid_to_info[fields[1]] = { x:y for x,y in zip(t_headers,fields)}
        for line in open(input.gtf):
            if line[0] == '#':
                continue
            line = line.rstrip()
            line = line.split('\t')
            if line[2] != 'transcript':
                continue
            info = {x.split()[0] : x.split()[1].strip('"') for x in line[8].strip('; ').split(';')}
            tid_to_info[info['transcript_id']] = {f:info[f] for f in t_headers}
        print(len(tid_to_info))
        r_headers = [
            'read_id',
            'strand',
            'start',
            'end',
            'length',
            'error-free_length',
            'read_identity',
            'type',
            'chimera'
        ]
        out_tsv   = open(output.tsv, 'w+')
        out_fastq = open(output.fastq, 'w+')
        rid = 0
        print('\t'.join(t_headers + r_headers), file=out_tsv)
        for input_fastq in input.fastqs:
            print(input_fastq)
            for line in open(input_fastq):
                if line[0] == '@':
                    tid = line.split()[1].split(',')[0]
                    print('@{}_{:08d}'.format(tid, rid), file=out_fastq)
                    rid+=1
                else:
                    print(line, end='', file=out_fastq)
                if line[0] != '@':
                    continue
                line = line.rstrip()
                line = line.split()
                rid_info = {x:'NA' for x in r_headers}
                rid_info['type'] = 'normal'

                field = line[0]
                rid_info['read_id'] = field[1:]

                line = line[1:]
                field = line[0]
                if field in ['junk_seq', 'random_seq']:
                    tid = 'non'
                    rid_info['type'] = field
                else:
                    field = field.split(',')
                    tid                = field[0].split('_')[0]
                    rid_info['strand'] = field[1][0]
                    rid_info['start']  = field[2].split('-')[0]
                    rid_info['end']    = field[2].split('-')[0]
                while line[1] == 'chimera':
                    if rid_info['chimera'] == 'NA':
                        rid_info['chimera'] = ''
                    line = line[2:]
                    rid_info['chimera'] += '{};'.format(line[0])
                    rid_info['type']    += ';chimera'
                try:
                    assert(len(line) == 4)
                except Exception:
                    print(line)
                    assert(len(line) == 4)
                while len(line)>1:
                    line = line[1:]
                    field = line[0]
                    field = field.split('=')
                    rid_info[field[0]] = field[1]
                out_line = list()
                out_line.extend([tid_to_info[tid][h] for h in t_headers])
                out_line.extend([rid_info[h] for h in r_headers])
                print('\t'.join(out_line), file=out_tsv)
        out_fastq.close()
        out_tsv.close()
