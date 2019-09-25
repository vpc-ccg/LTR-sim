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
        expand('{}/{{sample}}.cdna.fasta'.format(reads_d), sample=config['samples']),
        expand('{}/{{sample}}.cdna.reads.fastq'.format(reads_d), sample=config['samples']),


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
        for line in open(input.exp):
            line = line.rstrip().split('\t')
            tid = line[0]
            cnt = line[1]
            tid_to_rcnt[tid] = cnt
        # print(tid_to_rcnt)
        outfile = open(output.cdna, 'w+')
        for line in open(input.cdna):
            line = line.rstrip()
            if line[0] == '>':
                line = line[1:]
                line = line.split()
                tid = line[0].split('.')[0]
                cnt = tid_to_rcnt.get(tid, 0)
                record = [tid, 'depth={}'.format(cnt)] + line[1:]
                print('>{}'.format(' '.join(record)), file=outfile)
            else:
                print(line, file=outfile)
        outfile.close()

rule generate_reads:
    conda:
        'conda.env'
    input:
        badread = 'extern/Badread/badread-runner.py',
        cdna    = '{}/{{sample}}.cdna.fasta'.format(reads_d),
    output:
        fastq   = '{}/{{sample}}.cdna.reads.fastq'.format(reads_d),
    shell:
        '{input.badread} simulate --seed 42 --length 100000,0 --reference {input.cdna} --quantity 1x | gzip > {output.fastq}'
















#
