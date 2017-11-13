MashMap
========================================================================

MashMap is a fast and approximate long read (PacBio/ONT) mapper. It maps a read against a reference region if and only if its estimated alignment identity is above a specified threshold. It does not compute the alignments explicitly, but rather estimates a *k*-mer based [Jaccard similarity](https://en.wikipedia.org/wiki/Jaccard_index) using a combination of [Winnowing](http://www.cs.princeton.edu/courses/archive/spr05/cos598E/bib/p76-schleimer.pdf) and [MinHash](https://en.wikipedia.org/wiki/MinHash). This is then converted to an estimate of sequence identity using the [Mash](http://mash.readthedocs.org) distance. An appropriate *k*-mer sampling rate is automatically determined given minimum local alignment length and identity thresholds. The efficiency of the algorithm improves as both of these thresholds are increased.

Unlike traditional read mappers, MashMap does not compute gapped pairwise alignments. Instead it approximates mapping positions and identities using only *k*-mers. As a result, MashMap is both extremely fast and memory efficient, enabling rapid long-read mapping to large reference databases like NCBI RefSeq. We describe the full algorithm and report on speed, scalability, and accuracy of the software here: ["A fast approximate algorithm for mapping long reads to large reference databases"](http://biorxiv.org/content/early/2017/01/27/103812).

MashMap is in early development and currently reports only full-length mappings of reads to references. For split-read mapping, see Heng Li's [minimap](https://github.com/lh3/minimap), which is based on a similar idea, but does not provide identity estimates for the mapping targets it reports. We plan to add support for split-read mapping in future versions of MashMap.

## Installation
Follow [`INSTALL.txt`](INSTALL.txt) to compile and install MashMap.

## Usage

* Map set of long reads against a reference genome:
  ```sh
  mashmap -s reference.fna -q query.fa
  ```
  The output is space-delimited with each line consisting of query name, length,
  0-based start, end, strand, target name, length, start, end and mapping nucleotide
  identity.

* Map set of long reads against a list of reference genomes:
  ```sh
  mashmap --sl referenceList.txt -q query.fa
  ```
  File 'referenceList.txt' containing the list of reference genomes should contain path to the reference genomes, one per line.

## Parameters

For most of the use cases, default values should be appropriate. However, different parameters and their purpose can be checked using the help page `mashmap -h`. Important ones are mentioned below:

* Identity threshold (--perc_identity, --pi) : By default, its set to 85, implying read mappings with 85% identity should be reported. It can be set to 80% to account for more noisy read datasets.

* Minimum segment length (-m, --segLength) :  Default is 5,000 bp. Reads below this length are ignored. Mashmap provides guarantees on reporting local alignments twice the segment length.

## Release

Use the [latest release](https://github.com/marbl/MashMap/releases) for a stable version. In case your goal is to reproduce the results in the [Mashmap paper](http://biorxiv.org/content/early/2017/01/27/103812), you should use the very  [first version](https://github.com/marbl/MashMap/releases/tag/v1.0). 
