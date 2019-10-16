# LTR-sim
Long transcriptomic reads simulator

## Design of simulator:

### Input:

- Raw LTRs
- Reference transcriptome
- Reference transcriptome annotation

### Preprocessing

#### Mapping to transcriptome and genome

- Characterize rates of:
  - In vivo
    - [ ] poly-A tail length --> use mean of 60, stdev of 20, and min cutoff of 20
    - expression per transcript --> ok for now
    - gene fusion --> ignore for now
  - In vitro
    - [ ] degradation --> linear function with `m` slope `b` for offset to cut from the end of the reads
    - sequence-borne truncation (PCR with mid-transcript poly-A) --> duplicate ref contig per polyA > 8nt, cut after poly-A, and distrubute `depth` between them
    - palindromic chimerism (PCR strand switching)
    - concatenation by ligation chimerism (library)
    - reverse strand start (?) --> Badread
  - In silico
    - concatenation by processing chimerism
    - splitting by processing
- Badread quantification error-rate from genome alignment

#### Build transcriptome reference

- Poly-A tails
- Expression
- Gene fusion transcripts
- Degraded transcripts
- Palindromic transcripts
- Sequence-borne truncation
- Ligation-borne concatenation

### Read generation

#### Badread simulator

- Generate full length transcript 
  - Reverse strand start rate
  - Concatenation by processing chimerism rate
  - Badread sequencing error rate

### Postprocessing

- Signal-borne splitting

