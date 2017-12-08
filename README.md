MashMap
========================================================================

MashMap is a fast and approximate software for mapping long reads (PacBio/ONT) or assembly to reference genome(s). It maps a query sequence against a reference region if and only if its estimated alignment identity is above a specified threshold. It does not compute the alignments explicitly, but rather estimates a *k*-mer based [Jaccard similarity](https://en.wikipedia.org/wiki/Jaccard_index) using a combination of [Winnowing](http://www.cs.princeton.edu/courses/archive/spr05/cos598E/bib/p76-schleimer.pdf) and [MinHash](https://en.wikipedia.org/wiki/MinHash). This is then converted to an estimate of sequence identity using the [Mash](http://mash.readthedocs.org) distance. An appropriate *k*-mer sampling rate is automatically determined given minimum local alignment length and identity thresholds. The efficiency of the algorithm improves as both of these thresholds are increased.

Unlike traditional mappers, MashMap does not compute exact sequence alignments. Instead it approximates mapping positions and identities using only *k*-mers. As a result, MashMap is both extremely fast and memory efficient, enabling rapid query mapping to large reference databases like NCBI RefSeq. We describe the full algorithm (first version)  and report on speed, scalability, and accuracy of the software here: ["A fast approximate algorithm for mapping long reads to large reference databases"](http://biorxiv.org/content/early/2017/01/27/103812).

We have extended this software to compute split-read mappings as well. Based on the specified minimum length requirements for local alignments, it segments the query sequence accordingly to guarantee reporting of the requested local alignment boundaries, with high probability. In future, we plan to add an optional alignment support to generate base-to-base alignments.

## Installation
Follow [`INSTALL.txt`](INSTALL.txt) to compile and install MashMap. We also provide dependency-free linux and OSX binaries with each release for user convenience.

## Usage

* Map set of query sequences against a reference genome:
  ```sh
  mashmap -r reference.fna -q query.fa
  ```
  The output is space-delimited with each line consisting of query name, length,
  0-based start, end, strand, target name, length, start, end and mapping nucleotide
  identity.

* Map set of query seqences against a list of reference genomes:
  ```sh
  mashmap --rl referenceList.txt -q query.fa
  ```
  File 'referenceList.txt' containing the list of reference genomes should contain path to the reference genomes, one per line.

## Parameters

For most of the use cases, default values should be appropriate. However, different parameters and their purpose can be checked using the help page `mashmap -h`. Important ones are mentioned below:

* Identity threshold (--perc_identity, --pi) : By default, it is set to 85, implying mappings with 85% or more identity should be reported. For example, it can be set to 80% to account for more noisy long-read datasets or 95% for mapping human genome assembly to human reference.

* Minimum segment length (-s, --segLength) :  Default is 5,000 bp. Sequences below this length are ignored. Mashmap provides guarantees on reporting local alignments of length twice this value.

* Filtering options (-f, --filter_mode) : Similar to [delta-filter](http://mummer.sourceforge.net/manual/#filter) in nucmer, different filtering options are provided that are suitable for long read or assembly mapping. `-f map` is suitable for reporting the best mappings for long reads, whereas `-f one-to-one` is suitable for reporting orthologous mappings among all computed assembly to genome mappings.   

## Visualize

We provide a perl [script](scripts) for generating dot-plots to visualize mappings. It takes Mashmap's mapping output as its input. This script requires availability of gnuplot. Below is an example demonstrating mapping of [canu NA12878 human genome assembly](https://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md) (y-axis) to hg38 reference (x-axis).

<p align="center">
<img src="https://alurulab.cc.gatech.edu/sites/all/images/mashmap/dotplot_canu_hg38.png" height="500"/>
</p>

## Release

Use the [latest release](https://github.com/marbl/MashMap/releases) for a stable version. 
