configfile: 'config.yaml'

import os
import sys
import requests
import gzip

sys.setrecursionlimit(100)
def get_abs_path(path):
    abs_path = os.popen('readlink -f {}'.format(path)).read()
    return abs_path.rstrip("\r\n")

def make_slurm():
    os.makedirs('slurm'.format(outpath), mode=0o777, exist_ok=True)

# outpath = get_abs_path(config['outpath'])
outpath = config['outpath']
make_slurm()

map_d     = '{}/map'.format(outpath)
train_d   = '{}/train'.format(outpath)
reads_d   = '{}/reads'.format(outpath)
batches_d = '{}/batches'.format(reads_d)

localrules:
    all,
    batched_get_throughput

rule all:
    input:
        expand('refs/{species}/{species}.{ref_type}',
            species=config['refs'], 
            ref_type=['annot.gtf', 'cdna.fa', 'dna.fa', 'cdna.fa.fai', 'dna.fa.fai']),
        # expand('{}/{{sample}}.cdna.paf'.format(map_d),
        #        sample=config['samples'],),
        # expand('{}/{{sample}}.expression_rate.tsv'.format(train_d),
        #        sample=config['samples'],),
        # expand('{}/{{sample}}.cdna.fasta'.format(reads_d),
        #        sample=config['samples'],),
        # expand('{}/{{sample}}.cdna.polyA.fasta'.format(reads_d),
        #        sample=config['samples'],),
        # ['{dir}/{s}.cdna.polyA.degraded-{dl}.fasta'.format(
        #     dir=reads_d, s=s, dl=config['samples'][s]['degradation_level']) for s in config['samples']],
        # ['{dir}/{s}.cdna.polyA.degraded-{dl}/batch-{b}.fasta'.format(
        #     dir=batches_d, 
        #     s=s, 
        #     dl=config['samples'][s]['degradation_level'],
        #     b=b,
        #     )
        #     for s in config['samples'] for b in range(config['badread']['batches'])],
        # ['{dir}/{s}.cdna.polyA.degraded-{dl}/batch-{b}.cov.txt'.format(
        #     dir=batches_d, 
        #     s=s, 
        #     dl=config['samples'][s]['degradation_level'],
        #     b=b,
        #     )
        #     for s in config['samples'] for b in range(config['badread']['batches'])],
        # ['{dir}/{s}.cdna.polyA.degraded-{dl}/batch-{b}.reads.fastq'.format(
        #     dir=batches_d, 
        #     s=s, 
        #     dl=config['samples'][s]['degradation_level'],
        #     b=b,
        #     )
        #     for s in config['samples'] for b in range(config['badread']['batches'])],
        ['{dir}/{s}.L-{dl}.{ext}'.format(
            dir=reads_d, 
            s=s, 
            dl=config['samples'][s]['degradation_level'],
            ext=ext,
            )
            for s in config['samples'] for ext in ['fastq','tsv']],

rule git_badread:
    output:
        'extern/Badread/badread-runner.py'
    shell:
        'git clone --branch {tag} {url} extern/Badread'.format(tag=config['badread']['tag'], url=config['badread']['url'])

rule download_dna:
    output:
        outfile='refs/{species}/{species}.dna.fa'
    params:
        link = lambda wildcards: config['refs'][wildcards.species]['dna.fa']
    run:
        print('downloading:', params.link)
        os.system('wget {} -O {}.gz'.format(params.link, output))
        flag = True
        outfile = open(output.outfile, 'w+')
        print('writing to file:', output)
        for line in gzip.open('{}.gz'.format(output), 'rt'):
            if line[0] == '>':
                flag = not any((f in line) for f in config['refs'][wildcards.species]['dna_contig_filter'])
            if flag:
                outfile.write(line)
        outfile.close()
        os.remove('{}.gz'.format(output))

rule download_annot:
    input:
        index = 'refs/{species}/{species}.dna.fa.fai'
    output:
        outfile='refs/{species}/{species}.annot.gtf'
    params:
        link = lambda wildcards: config['refs'][wildcards.species]['annot.gtf']
    run:
        print('downloading:', params.link)
        os.system('wget {} -O {}.gz'.format(params.link, output))
        flag = True
        outfile = open(output.outfile, 'w+')
        contigs = {l.split()[0] for l in open(input.index)}
        print('writing to file:', output)
        for line in gzip.open('{}.gz'.format(output), 'rt'):
            if line[0]=='#' or line.split('\t')[0] in contigs:
                outfile.write(line)
        outfile.close()
        os.remove('{}.gz'.format(output))


rule download_cdna:
    input:
        annot = 'refs/{species}/{species}.annot.gtf'
    output:
        outfile='refs/{species}/{species}.cdna.fa'
    params:
        link = lambda wildcards: config['refs'][wildcards.species]['cdna.fa']
    run:
        print('downloading:', params.link)
        os.system('wget {} -O {}.gz'.format(params.link, output))
        tids = set()
        for l in open(input.annot):
            if l[0]=='#':
                continue
            l = l.rstrip().split('\t')
            if l[2]!='transcript':
                continue
            info = l[8]
            info = [x.strip().split(' ') for x in info.strip(';').split(';')]
            info = {x[0]:x[1].strip('"') for x in info}
            tids.add(info['transcript_id'])
        outfile = open(output.outfile, 'w+')
        print('writing to file:', output)
        for line in gzip.open('{}.gz'.format(output), 'rt'):
            if line[0]=='>':
                line = line.split()
                tid = line[0][1:].split('.')[0]
                line[0] = '>{}'.format(tid)
                line = ' '.join(line)+'\n'
                flag = tid in tids
            if flag:
                outfile.write(line)
        outfile.close()
        os.remove('{}.gz'.format(output))


rule index_ref:
    input:
        fasta = 'refs/{species}/{species}.{ref_type}'
    output:
        index = 'refs/{species}/{species}.{ref_type}.fai'
    wildcard_constraints:
        ref_type='cdna\.fa|dna\.fa'
    shell:
        'samtools faidx {input.fasta}'

rule sra_download:
    conda:
        'conda.env'
    output:
        reads  = 'samples/{sample}.fastq',
    params:
        SRA = lambda wildcards: config['samples'][wildcards.sample]['SRA']
    shell:
        "fastq-dump --split-files {params.SRA} --stdout > {output.reads}"

rule map_to_cdna:
    conda:
        'conda.env'
    input:
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        target = lambda wildcards: 'refs/{s}/{s}.cdna.fa'.format(s=config['samples'][wildcards.sample]['ref'])
    output:
        paf = '{}/{{sample}}.cdna.paf'.format(map_d),
    conda:
        'conda.env'
    threads:
        32
    shell:
        'minimap2 -Y -x map-ont -t {threads} {input.target} {input.reads} > {output.paf} '

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
            tid = line[5].split('_')[-1].split('.')[0]
            tid_to_rcnt[tid] = tid_to_rcnt.get(tid, 0) + 1
        outfile = open(output.exp, 'w+')
        for tid,rcnt in tid_to_rcnt.items():
            print('{}\t{}'.format(tid, rcnt), file=outfile)
        outfile.close()

rule sim_transcriptome:
    input:
        reads = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        cdna = lambda wildcards: 'refs/{s}/{s}.cdna.fa'.format(s=config['samples'][wildcards.sample]['ref']),
        gtf  = lambda wildcards: 'refs/{s}/{s}.annot.gtf'.format(s=config['samples'][wildcards.sample]['ref']),
        exp  = '{}/{{sample}}.expression_rate.tsv'.format(train_d)
    output:
        cdna = '{}/{{sample}}.cdna.fasta'.format(reads_d)
    run:
        print('Simulating transcriptome from \n{}'.format(input))
        
        genes = set([str(g) for g in config['samples'][wildcards.sample]['genes']])
        chroms = set([str(c) for c in config['samples'][wildcards.sample]['chroms']])
        any_chrom = 'All_chroms' in chroms
        any_gene = 'All_genes' in genes
        tids = set()
        print(chroms)
        for line in open(input.gtf):
            if line[0] == '#':
                continue
            line = line.rstrip()
            line = line.split('\t')
            if line[2] != 'transcript':
                continue
            chrom = line[0]
            if not (any_chrom or chrom in chroms):
                continue
            info = {x.split()[0] : x.split()[1].strip('"') for x in line[8].strip('; ').split(';')}
            if not (any_gene or info['gene_name'] in genes or info['gene_id'] in genes):
                continue
            tids.add(info['transcript_id'].split('.')[0])
        print('Looking at generating {} transcripts'.format(len(tids)))
        assert len(tids)>0
        tid_to_rcnt = dict()
        for line in open(input.exp):
            line = line.rstrip().split('\t')
            tid = line[0]
            if not tid in tids:
                continue
            cnt = int(line[1])
            tid_to_rcnt[tid] = cnt
        print('Only {} transcripts are experessed'.format(len(tid_to_rcnt)))
        assert len(tid_to_rcnt)>0
        outfile = open(output.cdna, 'w+')
        seq=list()
        flag = False
        for line in open(input.cdna):
            line = line.rstrip()
            if line[0] == '>':
                if len(seq) > 0:
                    print(''.join(seq), file=outfile)
                    seq = list()
                line = line[1:]
                line = line.split()
                tid = line[0].split('.')[0]
                if not tid in tid_to_rcnt:
                    flag = False
                    continue
                flag = True
                cnt = tid_to_rcnt[tid]
                record = [tid, 'depth={}'.format(cnt)] + line[1:]
                print('>{}'.format(' '.join(record)), file=outfile)
            elif flag:
                seq.append(line)
        if len(seq) > 0:
            print(''.join(seq), file=outfile)
        outfile.close()

rule add_poly_A:
    conda:
        'conda.env'
    input:
        polyA = config['exec']['polyA'],
        cdna  = '{}/{{sample}}.cdna.fasta'.format(reads_d),
    output:
        cdna_A = '{}/{{sample}}.cdna.polyA.fasta'.format(reads_d)
    params:
        mean   = config['polyA_config']['mean'],
        stddev = config['polyA_config']['stddev'],
        mincut = config['polyA_config']['mincut'],
        maxcut = config['polyA_config']['maxcut'],
        seed   = config['seed'],
    shell:
        'PYTHONHASHSEED=0 python3 {input.polyA} '
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
        'PYTHONHASHSEED=0 python3 {input.degrade}'
        '   {input.cdna_A} {output.cdna_A_degraded}'
        '   {params.slope} {params.intercept}'
        '   {params.seed}'

rule batch_fasta:
    input:
        split_file      = config['exec']['split_file'],
        cdna_A_degraded = '{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}.fasta'.format(reads_d)
    output:
        ['{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}/batch-{}.fasta'.format(batches_d,b) for b in range(int(config["badread"]["batches"]))]
    wildcard_constraints:
        sample='|'.join(config['samples']),
        degradation_level='|'.join(config['degradation_level'])
    params:
        lines_per_fasta = 2,
        batches = config["badread"]["batches"],
        out_pattern = lambda wildcards: '{}/{}.cdna.polyA.degraded-{}/batch-{{}}.fasta'.format(
            batches_d,
            wildcards.sample,
            config['samples'][wildcards.sample]['degradation_level']),
    shell:
        '{input.split_file} -i {input.cdna_A_degraded} -b {params.batches} -op {params.out_pattern} -l {params.lines_per_fasta}'

rule batched_get_throughput:
    input:
        cdna_A_degraded = '{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}/batch-{{bid}}.fasta'.format(batches_d)
    output:
        cdna_A_degraded_cov = '{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}/batch-{{bid}}.cov.txt'.format(batches_d)
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
        cdna_A_degraded     = '{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}/batch-{{bid}}.fasta'.format(batches_d),
        cdna_A_degraded_cov = '{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}/batch-{{bid}}.cov.txt'.format(batches_d)
    output:
        fastq               = '{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}/batch-{{bid}}.reads.fastq'.format(batches_d),
    wildcard_constraints:
        sample='|'.join(config['samples']),
        degradation_level='|'.join(config['degradation_level']),
        bid='|'.join([str(b) for b in range(config['badread']['batches'])])
    params:
        seed     = config['seed'],
    shell:
        'PYTHONHASHSEED=0 {input.badread} simulate --seed {params.seed} --length 100000,0'
        '   --reference {input.cdna_A_degraded} --quantity $(cat {input.cdna_A_degraded_cov})x'
        '   --chimeras 0 --random_reads 0 --junk_reads 0 --glitches 0,0,0'
        '   > {output.fastq}'

rule reads_info:
    input:
        fastqs = ['{}/{{sample}}.cdna.polyA.degraded-{{degradation_level}}/batch-{}.reads.fastq'.format(batches_d,x) for x in range(int(config["badread"]["batches"]))],
        gtf  = lambda wildcards: 'refs/{s}/{s}.annot.gtf'.format(s=config['samples'][wildcards.sample]['ref']),
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
