#rad-scripts

We foster the openness, integrity, and reproducibility of scientific research.

Various scripts and tools to process RAD-sequences

##Associated publication

> under-submission

##How to use this repository?

This repository host, various scripts and tools to process RAD-sequences. Feel free to adapt the scripts and tools, but remember to cite their authors!

To look at our scripts, **browse** through this repository. If you want to use some of the scripts, you will need to **clone** this repository. If you want to use our scripts for our own research, **fork** this repository and **cite** the authors.

##Available scripts (so far)

* `ddRAD_enzyme.pl`: A very simple bioperl script predicting the distribution of the fragment length after a ddRAD.

```
#Download a genome (e.g. Myotis lucifugus)
wget ftp://ftp.ensembl.org/pub/release-84/fasta/myotis_lucifugus/dna/Myotis_lucifugus.Myoluc2.0.dna.toplevel.fa.gz
gunzip Myotis_lucifugus.Myoluc2.0.dna.toplevel.fa.gz

#List of all the possible fragments produced by the association of SbfI and SphI
./ddRAD_enzyme.pl -i Myotis_lucifugus.Myoluc2.0.dna.toplevel.fa -2 SbfI -1 SphI > cuts.SbfI_SphI.txt

#List of all the possible fragments available for a ddRAD - (SbfI------SphI)
./ddRAD_enzyme.pl -i Myotis_lucifugus.Myoluc2.0.dna.toplevel.fa -2 SbfI -1 SphI -rad > rad.SbfI_SphI.txt
```

* `find_pattern.pl`: Test for diagnostic alleles or patterns between populations using raw STACKS RAD outputs

```
#extented help (PLEASE READ IT)
./find_pattern.pl

#Simple test
./find_pattern.pl --haplotypes batch_1.haplotypes.tsv --population pop.txt -v --group 2 -d

```

##Issues

If you have any problems with or questions about the scripts, please contact us through a [GitHub issue](https://github.com/pseudogene/rad-scripts/issues).


##Contributing

You are invited to contribute new features, fixes, or updates, large or small; we are always thrilled to receive pull requests, and do our best to process them as fast as we can.


##License and distribution


This code is distributed under the GNU GPL license v3. The documentation, raw data and work are licensed under a Creative Commons Attribution-ShareAlike 4.0 International License.​

