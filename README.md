MashMap
========================================================================

MashMap is a fast and approximate long read (ONT/PacBio) mapper. It maps a read against a reference region if and only if it's estimated alignment identity is above a specified threshold. It doesn't compute the alignments explicitly, but rather computes the [Jaccard similarity](https://en.wikipedia.org/wiki/Jaccard_index) measure and uses a sequence error model to estimate the identity. Jaccard similarity is estimated using an algorithm based on [MinHash](https://en.wikipedia.org/wiki/MinHash) and [Winnowing](http://www.cs.princeton.edu/courses/archive/spr05/cos598E/bib/p76-schleimer.pdf) techniques. We have tested its speed, scalability and accuracy by mapping ONT as well as PacBio data to complete RefSeq database. Unlike BWA or Bowtie, MashMap does not include the pairwise alignment steps. Using Jaccard similarity estimates, it only outputs the mapping coordinate of the read as well as the nucleotide identity. It is about as fast as [minimap](https://github.com/lh3/minimap), and >200x faster than BWA-MEM. 

Computing Jaccard using MinHash formula requires deciding an appropriate sketch size. MashMap autocomputes the sampling rate to ensure certain statistical significance (p-value) of the mapping results. Refer to [paper](http://biorxiv.org) to learn more about the algorithm and results.

## Installation
Follow [`INSTALL.txt`](INSTALL.txt) to compile and install MashMap.

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
  File 'referenceList.txt' containing the list of reference genomes should contain path to the reference genomes, one per line.
  
## Parameters

For most of the use cases, default values should be appropriate. However, different parameters and their purpose can be checked using the help page `mashmap -h`. Some of the important ones are mentioned below:

* Identity threshold (--perc_identity, --pi) : By default, its set to 85, implying read mappings with 85% identity should be reported. It can be set to 80% to account for more noisy read datasets.

* Minimum read length (-m, --minReadLen) :  Defult is 5,000 bp. This is set to 5K as the current average read lengths for both ONT and PacBio are >10K. Reads below this length are ignored. 

* Protein sequences (-a, --protein) : Use this parameter when mapping protein sequences. MashMap adjusts alphabet and k-mer size accordingly.

* Report all mappings, not just the best ones (--all)
