MashMap
========================================================================
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/mashmap.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/mashmap)
[![GitHub Downloads](https://img.shields.io/github/downloads/marbl/MashMap/total.svg?style=social&logo=github&label=Download)](https://github.com/marbl/MashMap/releases)

MashMap implements a fast and approximate algorithm for computing local alignment boundaries between long DNA sequences. It can be useful for mapping genome assembly or long reads (PacBio/ONT) to reference genome(s). Given a minimum alignment length and an identity threshold for the desired local alignments, Mashmap computes alignment boundaries and identity estimates using *k*-mers. It does not compute the alignments explicitly, but rather estimates an unbiased *k*-mer based [Jaccard similarity](https://en.wikipedia.org/wiki/Jaccard_index) using a combination of [minmers](https://www.biorxiv.org/content/10.1101/2023.05.16.540882v1) (a novel winnowing scheme) and [MinHash](https://en.wikipedia.org/wiki/MinHash). This is then converted to an estimate of sequence identity using the [Mash](http://mash.readthedocs.org) distance. An appropriate *k*-mer sampling rate is automatically determined using the given minimum local alignment length and identity thresholds. 

As an example, Mashmap can map a human genome assembly to the human reference genome in about one minute total execution time and < 4 GB memory using just 8 CPU threads, achieving more than an order of magnitude improvement in both runtime and memory over alternative methods. We describe the algorithms associated with Mashmap, and report on speed, scalability, and accuracy of the software in the publications listed [below](#publications). Unlike traditional mappers, MashMap does not compute exact sequence alignments. In future, we plan to add an optional alignment support to generate base-to-base alignments.

## MashMap3 important changes
- The automatic sampling rate has increased slightly relative to MashMap2, resulting in more accurate mappings and ANI prediction at the cost of more RAM. For a more dense sampling rate, users can provide the `--dense` flag. 
- The default output file format is now [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md), but the MashMap2 output format can be obtained with the `--legacy` flag. The `id` tag in the PAF output corresponds to the predicted ANI. 
- More details of the MashMap3 changes can be found in the [v3.0.1 changelog](https://github.com/marbl/MashMap/releases/tag/v3.0.1).

## Installation
Follow [`INSTALL.txt`](INSTALL.txt) to compile and install MashMap. We also provide dependency-free linux and OSX binaries for user convenience through the [latest release](https://github.com/marbl/MashMap/releases).

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

* Identity threshold (`--perc_identity`, `--pi`) : By default, it is set to 85, implying mappings with 85% or more identity should be reported. For example, it can be set to 80% to account for more noisy long-read datasets or 95% for mapping human genome assembly to human reference.

* Minimum segment length (`-s`, `--segLength`) :  Default is 5,000 bp. Sequences below this length are ignored. Mashmap provides guarantees on reporting local alignments of length twice this value.

* Filtering options (`-f`, `--filter_mode`) : Mashmap implements a [plane-sweep](https://en.wikipedia.org/wiki/Sweep_line_algorithm) based algorithm to perform the alignment filtering. Similar to [delta-filter](http://mummer.sourceforge.net/manual/#filter) in nucmer, different filtering options are provided that are suitable for long read or assembly mapping. Option `-f map` is suitable for reporting the best mappings for long reads, whereas `-f one-to-one` is suitable for reporting orthologous mappings among all computed assembly to genome mappings.   

* Sketch size (`-J`, `--sketchSize`) : This parameter sets the seed density of the winnowing scheme, gauranteeing that the minhash will be calculated from a sample of `sketchSize` k-mers for each segment. It is set automatically based on `--pi` but can be manually set as well. 

* Dense sketching (`--dense`) : This flag will increase the seed density substantially, resulting in a density of roughly `0.02 * (1 + (1 - pi) / .05)` where `pi` is the `perc_identity` threshold. This leads to longer runtimes and higher RAM usage, but significantly more accurate estimates of ANI. 

## Visualize

We provide a perl [script](scripts) for generating dot-plots to visualize mappings. It takes Mashmap's mapping output as its input. This script requires availability of gnuplot. Below is an example demonstrating mapping of [canu NA12878 human genome assembly](https://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md) (y-axis) to hg38 reference (x-axis).

<p align="center">
<img src="https://i.postimg.cc/HskJNzNg/readme-mashmap.jpg" height="300"/>
</p>

## Release

Use the [latest release](https://github.com/marbl/MashMap/releases) for a stable version. 

## <a name=“publications”></a>Publications

- **Bryce Kille, Erik Garrison, Todd Treangen, and Adam M Phillippy** ["Minmers are a generalization of minimizers that enable unbiased local Jaccard estimation"](https://doi.org/10.1101/2023.05.16.540882). *bioRxiv* (2023).

- **Chirag Jain, Sergey Koren, Alexander Dilthey, Adam M. Phillippy, and Srinivas Aluru**. ["A Fast Adaptive Algorithm for Computing Whole-Genome Homology Maps"](https://doi.org/10.1093/bioinformatics/bty597). *Bioinformatics (ECCB issue)*, 2018.

- **Chirag Jain, Alexander Dilthey, Sergey Koren, Srinivas Aluru, and Adam M. Phillippy**. ["A fast approximate algorithm for mapping long reads to large reference databases."](https://link.springer.com/chapter/10.1007/978-3-319-56970-3_5) In *International Conference on Research in Computational Molecular Biology*, Springer, Cham, 2017.
