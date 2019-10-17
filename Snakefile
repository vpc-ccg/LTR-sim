configfile: 'config.yaml'

import os

def get_abs_path(path):
    abs_path = os.popen('readlink -f {}'.format(path)).read()
    return abs_path.rstrip("\r\n")

def make_slurm():
    os.makedirs('slurm'.format(outpath), mode=0o777, exist_ok=True)

outpath = get_abs_path(config['outpath'])
make_slurm()
config['exec']['ltrsim'] = get_abs_path(config['exec']['ltrsim'])

map_d   = '{}/map'.format(outpath)
train_d = '{}/train'.format(outpath)
reads_d = '{}/reads'.format(outpath)

rule all:
    input:
        expand('{}/{{sample}}.cdna.polyA.fasta'.format(reads_d), sample=config['samples']),
        expand('{}/{{sample}}.cdna.polyA.reads.fastq'.format(reads_d), sample=config['samples']),
        expand('{}/{{sample}}.cdna.polyA.reads.tsv'.format(reads_d), sample=config['samples']),


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
        tid_to_rcnt = dict()
        for idx,line in enumerate(open(input.paf)):
            if idx % 100000 == 0:
                print('Read {} reads'.format(idx))
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
        tid_to_rcnt = dict()
        genes = set(config['genes'])
        for line in open(input.exp):
            line = line.rstrip().split('\t')
            tid = line[0]
            cnt = int(line[1])
            tid_to_rcnt[tid] = cnt
        outfile = open(output.cdna, 'w+')
        for line in open(input.cdna):
            line = line.rstrip()
            if line[0] == '>':
                gene_name_field_name = 'gene:'
                gene_name_field_sepa = '.'
                gene_name = line[line.find(gene_name_field_name):]
                gene_name = gene_name[len(gene_name_field_name) : gene_name.find(gene_name_field_sepa)]
                if not gene_name in genes:
                    flag = False
                    continue
                flag = True
                print(gene_name)
                line = line[1:]
                line = line.split()
                tid = line[0].split('.')[0]
                cnt = tid_to_rcnt.get(tid, 0)
                if cnt > 0:
                    record = [tid, 'depth={}'.format(max(cnt,0))] + line[1:]
                    print('>{}'.format(' '.join(record)), file=outfile)
            elif flag:
                if cnt > 0:
                    print(line, file=outfile)
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
    shell:
        'python {input.polyA} {input.cdna} {output.cdna_A} {params.mean} {params.stddev} {params.mincut} {params.maxcut}'


rule generate_reads:
    conda:
        'conda.env'
    input:
        badread = config['exec']['badread'],
        cdna_A = '{}/{{sample}}.cdna.polyA.fasta'.format(reads_d)
    output:
        fastq_gz   = '{}/{{sample}}.cdna.polyA.reads.fastq.gz'.format(reads_d),
    params:
        coverage = config['badread']['coverage']
    shell:
        '{input.badread} simulate --seed 42 --length 100000,0 --reference {input.cdna_A} --quantity {params.coverage} | gzip > {output.fastq_gz}'

rule decompress_reads:
    input:
        fastq_gz = '{}/{{sample}}.cdna.polyA.reads.fastq.gz'.format(reads_d),
    output:
        fastq    = '{}/{{sample}}.cdna.polyA.reads.fastq'.format(reads_d),
    shell:
        'zcat {input.fastq_gz} > {output.fastq}'

rule reads_info:
    input:
        fastq = '{}/{{sample}}.cdna.polyA.reads.fastq'.format(reads_d),
        gtf   = config['annotations']['gtf'],
    output:
        tsv   = '{}/{{sample}}.cdna.polyA.reads.tsv'.format(reads_d),
    run:
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
        outfile = open(output.tsv, 'w+')
        print('\t'.join(t_headers + r_headers), file=outfile)
        for line in open(input.fastq):
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
                tid                = field[0]
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
            print('\t'.join(out_line), file=outfile)
        outfile.close()












#
