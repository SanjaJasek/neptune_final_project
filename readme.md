# Neptune Computational Biology - Final Project

## Guidelines - you can delete this section before submission

This repository is a stub for your final project. Fork it, develop your project, and submit it as a pull request. Edit/ delete the text in this readme as needed.

Some guidelines and tips:

- Use the stubs below to write up your final project.

- For information on formatting text files with markdown, see https://guides.github.com/features/mastering-markdown/ . You can use markdown to include images in this document by linking to files in the repository, eg `![Figure 1](./Figure1.png?raw=true)`.

- The project must be entirely reproducible. In addition to the results, the repository must include all the data (or links to data) and code needed to reproduce the results.

- If you are working with unpublished data that you would prefer not to publicly share at this time, please contact me to discuss options. In most cases, the data can be anonymized in a way that putting them in a public repo does not compromise your other goals.

- Paste references (including urls) into the reference section, and cite them with the general format (Smith at al. 2003).

- Commit and push often as you work.

OK, here we go.

# k-mer mapper

## Introduction and Goals

The goal of my project is to answer the question, What is...?

The goal of my project is to write a k-mer mapper from scratch using python.


Procedure:
reference - I will download a transcriptome from the internet (fasta)
reads - generate with sim_reads (fastq)

1. k-merize reads
2. generate hash table = dictionary where: 
	* key=k-mer sequence
	* value = number of times k-mer appears in fastq
3. search reference with k-mers, each match gets k-mer value added to the sequence



The methods I will use to do this are...


The data I will use are (my own data/ data publicly available at YYY/ simulations)

## Methods

I downloaded cow transcriptome from NCBI. I'm using only the first 4 sequences for practice.
I simuleted reads with sim_reads with default settings:

sim_reads cow/cow_mRNA_only_4seq{,reads}.fasta

From that I use only first 21 reads for practice.

The tools I used were... See analysis files at (links to analysis files).

## Results

![Figure 1](./Figure1.png?raw=true)

In Figure 1...

## Discussion

These results indicate...

The biggest difficulty in implementing these analyses was...

If I did these analyses again, I would...

## References

Put references here

