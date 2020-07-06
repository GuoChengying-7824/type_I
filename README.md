# irepeat.upstream.pl

This software enables the user to determine the presence of inverted repeats in a region of DNA downstream of the gene of interest.

User Information: The Perl script irepeat.upstream.pl has five following input parameters: 

-geno: for the location of the genome file directory which contains fasta sequence files.

-gene: the gene sequence file (in fasta format) of type I downloaded from the REBASE database

-anno: annotation from the REBASE database (example: annotation.strand.xls)

-outpre: the filename of the output file.

-outdir: the directory name of the output file.

example: perl irepeat.upstream.pl -geno example/genome -gene example/gene/gene.txt -anno example/annotation/annotation.strand.xls -outpre test -outdir ./example/outdir
