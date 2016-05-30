# EpiCSA DNA Methylation Workshop

June 2016

Presenters: Jimmy Breen, Ben Mayne  

DNA methylation hands-on tutorial for the Epigenetics Consortium of South Australia's (EpiCSA) Bioinformatics workshop

## What you're need  

- Laptop (preferrably a Mac or Linux)

For this tutorial we'll need the following programs:

- NCBI SRA tookit
- bowtie2
- bismark (http://www.bioinformatics.babraham.ac.uk/projects/download.html#bismark)
- samtools/htslib/bcftools

### MacOSX

To download some tools that we will need on MacOSX, we recommend homebrew:

	/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

After installing homebrew, we can easily install our requirements by running:

	brew install homebrew/science/sambamba
	brew install homebrew/science/htslib
	brew install homebrew/science/bcftools
	brew install homebrew/science/bowtie2
	brew install homebrew/science/sratoolkit
	brew install wget
