# TeloVision

TeloVision is a python package which determines the presence of telomeres in genome assembly scaffolds. Subsequently the scaffolds are visualised along with their GC content. Telomeres will be visualised as black boxes at the ends of telomeres, while the absence of telomeres will be indicated by a white box at the en of telomeres. 

## Installation

TeloVision can be installed by first cloning this GitHub repository. Subsequently, go into the location where the repository was saved and use pip to install the package.
```bash
# Clone repository
git clone https://github.com/TimVerschuren/TeloVision.git

# Go into repository location
cd .../TeloVision/

# Install package with pip
pip install .
```

## Usage

TeloVision has 2 required inputs and 3 optional ones. An input and output file are required, while the k-mer size, minimum repeat length and sequence size can be altered if the user wishes so. If the optional inputs are not selected, default values will be used. 
```bash
usage: telovision [-h] -i INPUT -o OUTPUT [-k KMER_SIZE] [-r MIN_REPEAT_LENGTH] [-s SEQUENCE_SIZE]

TeloVision - Determine and visualise the presence of telomeres. For more information see: https://github.com/TimVerschuren/TeloVision

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Location of fasta file
  -o OUTPUT, --output OUTPUT
                        Name of output file
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        Size of k-mer sliding window. (default: 5000)
  -r MIN_REPEAT_LENGTH, --min_repeat_length MIN_REPEAT_LENGTH
                        Minimum length for a repetitive region to be counted as a telomere. (default: 30)
  -s SEQUENCE_SIZE, --sequence_size SEQUENCE_SIZE
                        Size of sequence taken from the 5' and 3' ends of the sequence for analysis. (default: 200)
```
