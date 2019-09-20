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
    - poly-A tail length
    - expression per transcript
    - gene fusion
  - In vitro
    - degradation
    - palindromic chimerism (PCR strand switching)
    - sequence-borne truncation (PCR with mid-transcript poly-A)
    - concatenation by ligation chimerism (library)
    - reverse strand start (?)
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

