#!/usr/bin/perl -w
use strict;

unless(scalar(@ARGV) == 2){
  die "usage: perl $0 <chr> <width>\n";
}
my $CHR = $ARGV[0];
my $WIDTH = $ARGV[1];
print "chr=$CHR; width=$WIDTH\n";
open(IN,"/hernandez/netapp/rhernandez/GEUVADIS/matrix.eqtl/expressions.genes.meta.matrix.eqtl.txt") or die "Cannot read /hernandez/netapp/rhernandez/GEUVADIS/matrix.eqtl/expressions.genes.meta.matrix.eqtl.txt\n";
my %GENE = ();
my %GNAME = ();
my $nGENE = 0;
while(<IN>){
  chomp;
  my @sp = split('\t',$_);
  if($sp[0] =~ / /){ #skip genes with spaces in names
    next; 
  }
  if($sp[1] eq "chr$CHR"){
    if(1 && (-e "/hernandez/netapp/rhernandez/GEUVADIS/GENOW/$WIDTH/$sp[0].$WIDTH.txt.012.gz") 
       && (-e "/hernandez/netapp/rhernandez/GEUVADIS/GENOW/$WIDTH/$sp[0].$WIDTH.txt.pos.gz")){
      next;
    }
    my $start = $sp[2] - $WIDTH;
    $start = 1 if($start <= 0);
    my $stop = $sp[3] + $WIDTH;
    if(exists($GENE{$start})){
      if($GENE{$start}[0] > $stop){
	next;
      }
    }
    $GENE{$start}[0] = $stop;
    $GENE{$start}[1] = $sp[0];
    $GNAME{$sp[0]}++;
    $nGENE++;
  }
}
close(IN);
print "nGENE = $nGENE\n";
if($nGENE == 0){
  exit;
}
my @SGENE = ();
my $nSGENE = 0;
foreach my $s (sort {$a<=>$b} keys %GENE){
  $SGENE[$nSGENE][0] = $s;
  $SGENE[$nSGENE][1] = $GENE{$s}[0];
  $SGENE[$nSGENE][2] = $GENE{$s}[1];
  $nSGENE++;
}

my %SAMP = ();
my $nSAMP = 0;
open(IN,"/hernandez/netapp/rhernandez/GEUVADIS/GENO/indiv.txt") or die "cannot read /hernandez/netapp/rhernandez/GEUVADIS/GENO/indiv.txt\n";
while(<IN>){
  chomp;
  $SAMP{$_}++;
  $nSAMP++;
}
close(IN);
print "nSAMP=$nSAMP\n";

my @SCOL = ();
my $nSCOL = 0;
open(IN,"/hernandez/netapp/share/tgp/phase3/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz") or die "cannot read /hernandez/netapp/share/tgp/phase3/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz\n";
close(IN);
open(IN,"zcat /hernandez/netapp/share/tgp/phase3/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz |") or die "cannot zcat /hernandez/netapp/share/tgp/phase3/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz\n";
my %FOUNDSAMP = ();
my $nlines = 0;
my $sindex = 0;
my %GENOS = ();
my %SITES = ();
my %PRINTED = ();
while(<IN>){
  $nlines++;
  if($_ =~ /\#\#/){
    next;
  }
  elsif($_ =~ /\#CHROM/){
    chomp;
    my @sp = split(' ', $_);
    for(my $i=0; $i<scalar(@sp); $i++){
      if(exists($SAMP{$sp[$i]})){
	$SCOL[$nSCOL] = $i;
	$FOUNDSAMP{$sp[$i]} = $nSCOL;
	$nSCOL++;
      }
    }
  }
  elsif($_ =~ /$CHR\t(\d+)\t.*\t[ACGT]\t[ACGT]\t100/){
    my $site = $1;
    if($nlines % 10000 == 0){
      print "read line $nlines, site=$site [$sindex]:$SGENE[$sindex][0]-$SGENE[$sindex][1]\n";
    }
    if($site >= $SGENE[$sindex][0]){
      if($site > $SGENE[$sindex][1]){ #print and move on
	if(exists($GENOS{$SGENE[$sindex][2]})){
	  print "site $site; finished $sindex/$nSGENE: $SGENE[$sindex][2]; ";
	  print "$GENOS{$SGENE[$sindex][2]}[0] / $SITES{$SGENE[$sindex][2]}[0] SNPs\n";
	  my $outfile="/hernandez/netapp/rhernandez/GEUVADIS/GENOW/$WIDTH/$SGENE[$sindex][2].$WIDTH.txt.012.gz";
	  open(OUT,">$outfile") or die "cannot write to $SGENE[$sindex][2].txt.012.gz";
	  close(OUT);
	  open(OUT,"| gzip - > $outfile") or die "cannot write to $SGENE[$sindex][2].txt.012.gz";
	  for(my $j=0; $j<$nSCOL; $j++){
	    for(my $i=1; $i<=$GENOS{$SGENE[$sindex][2]}[0]; $i++){
	      print OUT "$GENOS{$SGENE[$sindex][2]}[$i][$j]\t";
	    }
	    print OUT "\n";
	  }
	  close(OUT);
	  $outfile = "/hernandez/netapp/rhernandez/GEUVADIS/GENOW/$WIDTH/$SGENE[$sindex][2].$WIDTH.txt.pos.gz";
	  open(OUT,">$outfile") or die "cannot write to $SGENE[$sindex][2].txt.pos.gz";
	  close(OUT);
	  open(OUT,"| gzip - > $outfile") or die "cannot write to $SGENE[$sindex][2].txt.pos.gz";
	  for(my $i=1; $i<=$SITES{$SGENE[$sindex][2]}[0]; $i++){
	    print OUT "$SITES{$SGENE[$sindex][2]}[$i]\n";
	  }
	  close(OUT);
	}
	$PRINTED{$SGENE[$sindex][2]}++;
	if(exists($GENOS{$SGENE[$sindex][2]})){
	  delete $GENOS{$SGENE[$sindex][2]};
	  delete $SITES{$SGENE[$sindex][2]};
	}
	$sindex++;
	if($sindex >= $nSGENE){
	  last;
	}
	redo;
      }
      else{
	chomp;
	my @sp = split(' ',$_);
	my @genos = ();
	my $freq = 0;
	for(my $i=0; $i<$nSCOL; $i++){
	  my @sp2 = split('',$sp[$SCOL[$i]]);
	  $genos[$i] = $sp2[0]+$sp2[2];
	  $freq += $sp2[0]+$sp2[2];
	}
	if($freq == 0 || $freq == 2*$nSCOL){
	  next;
	}
	for(my $i=$sindex; $i<$nSGENE; $i++){
	  if($site >= $SGENE[$i][0] && $site <= $SGENE[$i][1]){
	    $GENOS{$SGENE[$i][2]}[0]++;
	    $SITES{$SGENE[$i][2]}[0]++;
	    $SITES{$SGENE[$i][2]}[$SITES{$SGENE[$i][2]}[0]] = $site;
	    for(my $j=0; $j<$nSCOL; $j++){
	      $GENOS{$SGENE[$i][2]}[$GENOS{$SGENE[$i][2]}[0]][$j] = $genos[$j];
	    }
	  }
	}
      }
    }
  }
}
close(IN);
foreach my $g (keys %GNAME){
  if(!exists($PRINTED{$g})){
    open(OUT,">/hernandez/netapp/rhernandez/GEUVADIS/GENOW/$WIDTH/$SGENE[$sindex][2].$WIDTH.txt.012.gz") or die "cannot write to $SGENE[$sindex][2].txt.012.gz";
    close(OUT);
    open(OUT,"| gzip - > /hernandez/netapp/rhernandez/GEUVADIS/GENOW/$WIDTH/$SGENE[$sindex][2].$WIDTH.txt.012.gz") or die "cannot write to $SGENE[$sindex][2].txt.012.gz";
    for(my $j=0; $j<$nSCOL; $j++){
      for(my $i=1; $i<=$GENOS{$g}[0]; $i++){
	print OUT "$GENOS{$g}[$i][$j]\t";
      }
      print OUT "\n";
    }
    close(OUT);
    open(OUT,">/hernandez/netapp/rhernandez/GEUVADIS/GENOW/$WIDTH/$g.$WIDTH.txt.pos") or die "cannot write to $g.$WIDTH.txt.pos";
    for(my $i=1; $i<=$SITES{$g}[0]; $i++){
      print OUT "$SITES{$g}[$i]\n";
    }
    close(OUT);
  }
}

