#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my($genomedir,$genefile,$annotationfile,$outprefix,$outdir);

GetOptions(
	"help|h|?" => \&USAGE,
	"geno:s" => \$genomedir,
	"gene:s" => \$genefile,
	"anno:s" => \$annotationfile,
	"outpre:s" => \$outprefix,
	"outdir:s" => \$outdir,
	) or &USAGE;

&USAGE unless ($genomedir and $genefile and $annotationfile and $outprefix and $outdir);
my($name,%hash,%hash2,@irs_l);
my @list = glob("$genomedir/*.fasta");
foreach  $name(@list){
	if($name =~/$genomedir\/(\S+).fasta/){
		my $title = $1;
		print  "$title\n";
	}
local $/ = ">";
open (GENOME, "$name") || die "Can't open GENOME! \n" ;
while(<GENOME>){
	chomp;
	my ($titleline, $sequenceGenome) = split(/\n/,$_,2);
	next unless ($sequenceGenome && $titleline);
	$titleline = (split(/[\s]+/,$titleline))[0];
	my $ncbiID = (split(/\./,$titleline))[0];
	$sequenceGenome =~ s/\s//g;
	$sequenceGenome =~ s/\n//g;
	study($sequenceGenome);
	$hash{$ncbiID} = $sequenceGenome;
}

open (GENE,$genefile) || die "Can't open GENE! \n" ;
while(<GENE>){
	chomp;
	my ($lineID, $sequenceGene) = split(/\n/,$_,2);
	next unless ($lineID && $sequenceGene);
	my $geneID = (split(/\s/,$lineID))[0];
	$sequenceGene =~ s/\s//g;
	$sequenceGene =~ s/\n//g;
	study($sequenceGene);
	$hash2{$geneID} = $sequenceGene;
}
}	
@irs_l = 20..500;
open (GFF, $annotationfile) || die "Can't open GFF! \n" ;
open (OUT, ">$outdir/$outprefix.xls") || die "Can't creat OUTFILE! \n" ;
open (OUT2, ">$outdir/$outprefix.defeated.gff");
print OUT "#GeneName\tgeneStart\tgeneEnd\tOrganism\tGenBank\tgeneLength\tInverted_Repeats_number\tInverted_Repeats_Length\tGene_seq(5'-3')\tInverted_Repeats_gene_Start\tInverted_Repeats_gene_End\tInverted_Repeats_upstream_seq(5'-3')\tInverted_Repeats_upstream_Start\tInverted_Repeats_upstream_End\n";

while(<GFF>){
	chomp;
	my ($GeneName,$geneStart,$geneEnd,$Organism,$GenBank,$Strand) = split(/\t/,$_);
	if($geneStart ne "-" && $geneEnd ne "-"){
		my $geneLength = $geneEnd - $geneStart + 1;
		if(exists $hash2{$GeneName} && exists $hash{$GenBank}){
			my $seqGene = $hash2{$GeneName};
			my $seqGenome = $hash{$GenBank};	
			my $len = length($seqGenome);
			my $reseqGene = reverse($seqGene);
			$reseqGene =~ tr/ATGC/TACG/;
			if($geneStart < $len){
				if($geneStart >= 30000){
					my $upstream = substr($seqGenome,$geneStart-30000,30000);
				}
				else{
					my $upstream = substr($seqGenome,0,$geneStart);
				}
				my $ir_number = 1;
				my $repSeq = "";
				my %locations;
				for(my $i=0; $i<scalar(@irs_l); $i++){
					my $minreps = $irs_l[$i];
					for(my $l=0; $l<$geneLength-$minreps; $l++){
						if($Strand eq "-"){
					       		my $repSeq = substr($seqGene,$l,$minreps);
							if(my $upstream =~ /($repSeq)/ig){
								my $ir = uc($1);
								my $irLength = length($ir);
								my $irStart = index($seqGenome,$ir);
								my $irEnd = $irStart + $irLength;
								my $irTran = reverse($ir);
								$irTran =~ tr/ATGC/TACG/;
								my $irgeneStart = index($seqGene,$ir);
								my $irgeneEnd = $irgeneStart + $irLength;
							#	if (&novel($irgeneStart, \%locations)){
									print OUT join("\t", $GeneName, $geneStart, $geneEnd, $Organism, $GenBank,
									$geneLength, $ir_number++, $irLength, $ir,$irgeneStart, $irgeneEnd,
									 $irTran, $irStart,$irEnd), "\n";
							}
						}
						else{ ##if($Strand eq "+")
							my $repSeq = substr($reseqGene,$l,$minreps);
							if(my $upstream =~ /($repSeq)/ig){
								my $ir = uc($1);
								my $irLength = length($ir);
								my $irStart = index($seqGenome,$ir);
								my $irEnd = $irStart + $irLength;
								my $irTran = reverse($ir);
								$irTran =~ tr/ATGC/TACG/;
								my $irgeneStart = index($seqGene,$irTran);
								my $irgeneEnd = $irgeneStart + $irLength;
								#       if (&novel($irgeneStart, \%locations)){
									print OUT join("\t", $GeneName, $geneStart, $geneEnd, $Organism, $GenBank,
									$geneLength, $ir_number++, $irLength, $ir,$irgeneStart, $irgeneEnd,
									$irTran, $irStart,$irEnd), "\n";
							}
						}

					}
				}
			}

			else{
				print OUT2 "$_\t$geneStart\t$len\n";
			}
		}
	}
}
#######-----------------------------------------------------------########
# sub function
sub novel {
		my($position, $locationsref) = @_;
		if(defined $locationsref->{$position}){
				return undef;
		}
		else{
				$locationsref->{$position} = 1;
				return 1;
		}
}
sub USAGE {
	my $usage=<<"USAGE";
	Program:     $0
	Description: repeat upstream;
	Usage:
		-geno    Genome file directory  must be given
		-gene    Gene file              must be given
		-anno    Annotation file        must be given
		-outpre  Output file prefix     must be given
		-outdir  Output file directory  must be given
		-h       Help document
USAGE
	print $usage;
	exit;
}
