MashMap
========================================================================


## Installation
Follow [`INSTALL.txt`](INSTALL.txt) to compile and install the code

## Usage

* Map set of long reads against a reference genome:
  ```sh
  mashmap -s reference.fna -q query.fa -o output.txt
  ```
  The output is space-delimited with each line consisting of query name, length,
  0-based start, end, strand, target name, length, start, end, mapping nucleotide
  identity, count of shared sketch elements and the sketch size.

* Map set of long reads against a list of reference genomes:
  ```sh
  mashmap --sl referenceList.txt -q query.fa -o output.txt
  ```
  File containing the list of reference genomes should contain path to reference 
  genomes, one per line.
