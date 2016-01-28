#!/usr/bin/perl 

###########################################################################################
# a script to annotate your assemblies via a "weak" recip blast method & to identify and  #
# break up chimeric contigs                                                               #
# external dependencies: blastall/makeblastdb (>2.2.17), exonerate, cd-hit-est            #
# This script DOES NOT require bioperl. 												  #
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 29 January 2012                #
# Added command flags by Ke Bi. July 27 2013                                              #
# Minor bugs fixed. Now use makeblastdb instead of formatdb                               # 
# Some minor modifications in subroutine orfFinesse                                       #
# Blasting is one directional instead of reciprocal                                       #                                                #
###########################################################################################

use warnings;
use strict;
use List::Util qw[min max];
use File::Basename;
use Getopt::Std;

die(qq/
8annotationTranscriptome.pl [options] 
external dependencies: blastall\/makeblastdb (>2.2.17), exonerate, cd-hit-est
options:
-a  FILE  folder with all trinity assemblies (named as libray_name.fasta). MUST PROVIDE A FULL PATH HERE!
-b  FILE  reference protein database (e.g. Xenopus_tropicalis.fa)
-c  INT   evalue for blast search [1e-20]
-e  INT   how many processors can you use for blasting [4]
-f  FILE  a file that contains all wikigene names\/descriptions associated with reference protein database (using Ensemble BioMart tool) 
examples for gene name file (-f)
Ensembl Protein ID, Associated Gene Name, WikiGene Description
ENSACAP00000010213,FAM204A,uncharacterized protein C10orf84-like
ENSACAP00000020968,,
ENSACAP00000009946,ENO4,
ENSACAP00000007831,,importin-4-like
\n\n/) unless (@ARGV);

my %opts = (a=>undef, b=>undef, c=>1e-20, d=>undef, e=>4,f=>undef);
getopts('a:b:c:d:e:f:', \%opts);

my $dir;
 
if ($opts{a} =~ m/\/$/ ){
	$dir = $opts{a}; 
	}
else {
	$dir = $opts{a} . "/";
	}

my $np = $opts{e};
my @assemblies = <$dir*fasta>;
my $dbP = $opts{b};
my $evalue = $opts{c};

my $name=$opts{f};

my %protein;
open(IN, "<$dbP");
my $id;
while (<IN>){ 
	chomp ($_);
	if ($_ =~ m/^>(\S+)/){
		$id = $1;
		}
	else {
		$protein{$id}{'seq'} .= $_;
		}
	}
close(IN);


########################
# run the subroutines! #
########################

#formats the protein database unless it already has been done
unless (-f $dbP . '.pin') {
	my $call = system("makeblastdb -in $dbP -dbtype prot");
	}

foreach my $assembly (@assemblies) {
    my $as = $1 if basename($assembly) =~ m/(\S+).fasta/; 
    my $resDir =$dir.$as;  
    mkdir $resDir unless -d $resDir;
    my $seqout = $assembly . ".annotated";
    print "Doing assembly $assembly now!\n";
    unless (-f $seqout) {
		my $outfile1 = blastProteins($assembly); 
	        print "Done blasting 1 assembly $assembly now!\n";
		my $outfile2 = chimericTest($assembly,$outfile1,\%protein);
		my $outfile3 = blastProteins($outfile2);

		print "Done blasting 2 assembly $assembly now!\n";
		my $recip = recipBlast($outfile2);
		print "Done recip blasting assembly $assembly now!\n";
	       	
		my ($seqref,$annoref) = makeHash($outfile2,$dbP,$name); 
		my $seq = annotateProt($outfile3, $seqref, $annoref, $recip, \%protein);
		print "Done annotating assembly $assembly now!\n";
		$seq = orfFinesse($seq,\%protein);
	        

########################
# report the results!  #
########################

		my $seqout = $assembly . ".annotated";
		open (SEQOUT, ">$seqout");
		my %seq = %$seq;
		foreach my $id (sort {$a cmp $b} keys %seq) {
		    print SEQOUT ">", $id, "\t";
		    if ($seq{$id}{'info'}) {
				print SEQOUT $seq{$id}{'info'}, "\t" if $seq{$id}{'info'};
	      
				print SEQOUT $seq{$id}{'match'}[0]{'gene'}, "\t" if $seq{$id}{'match'};

				print SEQOUT $seq{$id}{'abbr'}, "\t" if $seq{$id}{'abbr'};
		
				print SEQOUT $seq{$id}{'desc'}, "\t" if $seq{$id}{'desc'};

				print SEQOUT $seq{$id}{'match'}[0]{'eval'}, "\t" if $seq{$id}{'match'};

			    }
		    print SEQOUT "\n";
		    print SEQOUT $seq{$id}{'seq'}, "\n";
			}
		#my $call = system("rm error.log");
		}
	system ("mv $dir$as'.fasta.'*   $resDir ");
	system ("mv $dir/'interruptedORFs.fa' $resDir");
	}

###########################
# behold the subroutines  #
###########################

sub recipBlast {
	my ($seq) = @_;	
	my $recip = $seq . '.recipBlast.out';
	my $call1 = system("makeblastdb -in $seq -dbtype nucl");
	my $call2 = system("tblastn -db $seq -query $dbP -num_threads $np -evalue $evalue -outfmt 6 -out $recip -max_target_seqs 10");

	my %r; my $score;
	my $recip1 = $seq . '.recipBlast_sorted.out';
	system ("sort -k 1,1 -k 12n,12r -k 11n,11r $recip > $recip1");
	open(IN, "<$recip1");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		unless ($r{$d[0]}) {
			$r{$d[0]}{$d[1]} = $d[10];
			}
		}

	close(IN);
	#unlink($recip);
	#unlink($recip1);
	my $call3 = system("rm $seq.n*");
	return(\%r);
	}

sub makeHash {
	my ($assembly,$database, $name) = @_;
	my (%seq, $id);
	
	open(IN, "<$assembly");
	while(<IN>){
		chomp(my $line = $_);
		if ($line =~ m/>(.*)/){
			$id = $1;
			}
		else {
			$seq{$id}{'seq'} .= $line;
			}
		}
	close(IN);

	my %anno;
	open(DATA, "<",$database);
	while(<DATA>){
	  chomp(my $line = $_);
	  if ($line =~ m/^>(\S+)/) {
	    my $ID = $1;
	    $anno{$ID} = {'gene'=>'NA', 'info' => 'NA'};	    
	  }
	}
	close DATA;
	
	open(NAME, "<$name");
	# get rid of the header
	chomp(my $line = <NAME>);
	while(<NAME>){
		chomp(my $line = $_);	  
		my @d = split('\t', $line);	
	      	$anno{$d[0]} = {'gene' => $d[1], 'info'=> $d[2]};
	     # print $id, "\t", $anno{$id}{'gene'},"\t",$anno{$id}{'info'},"\n";
	      
	    }	
	close(NAME);	 
	return(\%seq,\%anno)
      }

sub blastProteins {
  my ($assembly) = @_;
  
  my $masterout = $assembly . '.blast.out';
  my $call = system("blastx -db $dbP -query $assembly -num_threads $np -evalue $evalue -outfmt 6 -out $masterout -max_target_seqs 20");
  return($masterout);	
}

sub parseBlast {
  my ($seq, $out) = @_;
  
  my %seq = %$seq;
  
	my $tracker;
  open(IN, "<$out");
  while(<IN>) {
    chomp(my $line = $_);
    my @d = split(/\t/,$line);
    #if the line is for the same contig as before or not
    if ($seq{$d[0]}{'match'}) {
      #want to consider other high scoring matches
      unless ($d[1] eq $seq{$d[0]}{'match'}[$tracker]{'gene'}) {				
	#only consider addt'l matches if evalue is good and if there isn't a sharp decline in quality of match
	$d[10] = 1e-200 if $d[10] =~ m/^0.0$/;
	if ($d[10]/$seq{$d[0]}{'match'}[0]{'eval'} <= 10000) {
	  if (scalar(@{$seq{$d[0]}{'match'}}) < 10){
	    push @{$seq{$d[0]}{'match'}}, {'gene' => $d[1], 'eval' => $d[10]};
	    $tracker++;
	  }
	}	
      }			
    }
    #this is a new contig	
    else {
      #will be dividing by this later, so it cannot equal 0!
      $d[10] = 1e-200 if $d[10] =~ m/^0.0$/; 
      push @{$seq{$d[0]}{'match'}}, {'gene' => $d[1], 'eval' => $d[10]};	
      $tracker = 0;
    }
		}											
  close(IN);
  
  return(\%seq);
}

sub annotateProt {
	my ($out, $seq, $anno,$recip, $protein) = @_;
	
	my %pro =  %{$protein};
	$seq = parseBlast($seq,$out);	
	
	my %seq = %$seq; my %anno = %$anno; my %recip = %$recip; 
	
	#need to use exonerate to define utr etc; call to external program
	#my $db = Bio::DB::Fasta->new($dbP);
		
	foreach my $id (keys %seq) {
		if ($seq{$id}{'match'}) {
			if ($recip{$seq{$id}{'match'}[0]{'gene'}}{$id}) {
				#exonerate first
				my $query = "query.fa";
				my $target = "target.fa";
				open(QUERY, ">$query");
				open(TARGET, ">$target");
				print QUERY ">$id\n$seq{$id}{'seq'}\n";
				my $proid =$seq{$id}{'match'}[0]{'gene'};
				my $seq = $pro{$seq{$id}{'match'}[0]{'gene'}}{'seq'};
			       
				my $protlength = length($seq);
				my $contiglength = length($seq{$id}{'seq'});			
				print TARGET ">$proid\n$seq\n";
		       
				close(QUERY); close(TARGET);					
				my @call3 = `exonerate $target $query -m protein2genome --showalignment no --showcigar 0`;
				unlink($query); unlink($target);
				my $info;
			
				my @match;
				my @prot5;
				my @prot3;
				foreach (@call3) {
					chomp (my @line = split /\s+/, $_);
					if ($line[0] =~m /vulgar/) {
						push @match, $line[6]+1;
						push @match, $line[7];
						push @prot3, $line[3];
						push @prot5, $line[2];
						}
					}
				my $start = min (@match);
				my $end = max (@match);
				my $fiveu = min (@prot5);
				my $threeu = max (@prot3);
				my @d = split(/\s+/,$call3[2]);
				if ($d[8]){
					if ($d[8] eq '+' || $d[8] eq '-') {
						if ($d[8] eq '+') {
							#this is in 5->3
							#identify gene start and stop
							$info = 'gs' . $start . '_ge' . $end;
							if ($threeu/$protlength > 0.9 || $protlength - $threeu < 11) {
								#yes, i am going to call it a  3' utr
								my $utr3 = $end+1;
								$info = $info . '_3u' . $utr3;
								}
							if ($fiveu/$protlength < 0.1 || $fiveu < 11) {
								#yes, i am going to call it a 5' utr
								my $utr5 = $start - 1;
								$info = '5u' . $utr5 . "_" . $info;
								}
							}		
						else {
							#this is in 3->5;
							$seq{$id}{'seq'} = reverse($seq{$id}{'seq'});
							$seq{$id}{'seq'} =~ tr/ATGC/TACG/;	
						
							my $gs = $contiglength - $end+2;
							my $ge = $contiglength - $start;
							$info = 'gs' . $gs . '_ge' . $ge;
							if ($threeu/$protlength > 0.9 || $protlength - $threeu < 11) {
								#yes, i am going to call it a  3' utr
								my $utr3 = $ge+1;
								$info = $info . '_3u' . $utr3;
								}
							if ($fiveu/$protlength < 0.1 || $fiveu < 11) {
								#yes, i am going to call it a 5' utr
								my $utr5 = $gs - 1;
								$info = '5u' . $utr5 . '_' . $info;
						    	}				
							}
				    	}
					}
				#a blast match but no exonerate match?	
				if (!$d[8]) {
					print "Huh, this is odd. BLAST hit but not exonerate hit for $id?\n";
					}

				$seq{$id}{'info'} = $info;
				$seq{$id}{'desc'} = $anno{$seq{$id}{'match'}[0]{'gene'}}{'info'} if ($anno{$seq{$id}{'match'}[0]{'gene'}});
				$seq{$id}{'abbr'} = $anno{$seq{$id}{'match'}[0]{'gene'}}{'gene'} if ($anno{$seq{$id}{'match'}[0]{'gene'}}) ;

				delete $recip{$seq{$id}{'match'}[0]{'gene'}}{$id};
				}
			}				
		}
	
	foreach my $id (keys %recip) {
	  foreach my $contig (keys %{$recip{$id}}) {
	    if ($recip{$id}{$contig} < $evalue) {
	    	unless ($seq{$contig}{'info'}) {	
		
				my $query = "query.fa";
				my $target = "target.fa";
				open(QUERY, ">$query");
				open(TARGET, ">$target");
				print QUERY ">$id\n$seq{$contig}{'seq'}\n";
				#my $proid =$seq{$id}{'match'}[0]{'gene'};
				my $seq = $pro{$id}{'seq'};
				my $protlength = length($seq);
				my $contiglength = length($seq{$contig}{'seq'});			
				print TARGET ">$id\n$seq\n";
		
				close(QUERY); close(TARGET);					
				my @call3 = `exonerate $target $query -m protein2genome --showalignment no --showcigar 0`;
				unlink($query); unlink($target);
				my $info;
		
				my @match;
				my @prot5;
				my @prot3;
				foreach (@call3) {
					chomp (my @line = split /\s+/, $_);
					if ($line[0] =~m /vulgar/) {
		    			push @match, $line[6]+1;
		    			push @match, $line[7];
		    			push @prot3, $line[3];
		   				push @prot5, $line[2];
		  				}
					}

				my $start = min (@match);
				my $end = max (@match);
				my $fiveu = min (@prot5);
				my $threeu = max (@prot3);
				my @d = split(/\s+/,$call3[2]);
				if ($d[8]){
					if ($d[8] eq '+' || $d[8] eq '-') {
		    			if ($d[8] eq '+') {
		    				#this is in 5->3
							#identify gene start and stop
		    				$info = 'gs' . $start . '_ge' . $end;
		   					if ($threeu/$protlength > 0.9 || $protlength - $threeu < 11) {
								#yes, i am going to call it a  3' utr
								my $utr3 = $end+1;
								$info = $info . '_3u' . $utr3;
		    					}
		    				if ($fiveu/$protlength < 0.1 || $fiveu < 11) {
								#yes, i am going to call it a 5' utr
								my $utr5 = $start - 1;
								$info = '5u' . $utr5 . "_" . $info;
		    					}
		   	 				}		
		    			else {
		    				#this is in 3->5;
		    				$seq{$contig}{'seq'} = reverse($seq{$contig}{'seq'});
		    				$seq{$contig}{'seq'} =~ tr/ATGC/TACG/;	
		      
		    				my $gs = $contiglength - $end+2;
		      				my $ge = $contiglength - $start;
		      				$info = 'gs' . $gs . '_ge' . $ge;
		    				if ($threeu/$protlength > 0.9 || $protlength - $threeu < 11) {
								#yes, i am going to call it a  3' utr
								my $utr3 = $ge+1;
								$info = $info . '_3u' . $utr3;
		      					}
		    				if ($fiveu/$protlength < 0.1 || $fiveu < 11) {
								#yes, i am going to call it a 5' utr
								my $utr5 = $gs - 1;
								$info = '5u' . $utr5 . '_' . $info;
		    					}				
		    				}
		  				}
					$seq{$contig}{'info'} = 'one_way_blast_' . $info;
					$seq{$contig}{'desc'} = $anno{$id}{'info'} if ($anno{$id});
					$seq{$contig}{'abbr'} = $anno{$id}{'gene'} if ($anno{$id});
					$seq{$contig}{'match'}[0]{'gene'} = $id;
					$seq{$contig}{'match'}[0]{'eval'} = $recip{$id}{$contig};
					}
				#a blast match but no exonerate match?	
				if (!$d[8]) {
		  			print "Huh, this is odd. BLAST hit but not exonerate hit for $id?\n";
					}
	    		}
	  		}		
		}
	}
	return(\%seq)	
	}	
	
sub chimericTest {
	my ($seq,$out,$protein) = @_;
	my %pro =  %{$protein};
	my $che = $seq . ".chimeric_contigs.fasta";
	open (CHE,">", $che);
	my %seq; my $id;
	open(IN, "<$seq");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			$id = $1;
			}
		else {
			$seq{$id} .= $line;
			}
		}
	close(IN);

	#parses the blast output
	open(OUT, "<$out");
	my %match;
	my %chimera;
	my $t;
	while(<OUT>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
	
		#it is a gene that i already have	
		if ($match{$d[0]}) {	
			#it is a match to that gene i already have
			if ($match{$d[0]}[$t]{'gene'} eq $d[1]) {
				push @{$match{$d[0]}[$t]{'start'}}, $d[6];
				#ends
				push @{$match{$d[0]}[$t]{'end'}}, $d[7];
				}
			#it is a new match to that gene
			else {
				$t++;
				push @{$match{$d[0]}[$t]{'start'}}, $d[6];
				#ends
				push @{$match{$d[0]}[$t]{'end'}}, $d[7];
				$match{$d[0]}[$t]{'gene'} = $d[1];
				}
			}
		#it is a new gene
		else {
			$t = 0;
			#starts
			push @{$match{$d[0]}[$t]{'start'}}, $d[6];
			#ends
			push @{$match{$d[0]}[$t]{'end'}}, $d[7];
			$match{$d[0]}[$t]{'gene'} = $d[1];
			}
		}
	close(OUT);	
	
	#figures out which contigs appear to be chimeric	
	foreach my $c (sort {$a cmp $b} keys %match) {
		my %l;
		my @m = @{$match{$c}};
	
		my $start = min(@{$m[0]{'start'}},@{$m[0]{'end'}});
		my $end = max(@{$m[0]{'start'}},@{$m[0]{'end'}});
		for (my $j = $start; $j <= $end; $j++) {
			$l{$c}{$j}++;
			}
		
		for (my $i = 1; $i < scalar(@m); $i++) {	
			my $start = min(@{$m[$i]{'start'}},@{$m[$i]{'end'}});
			my $end = max(@{$m[$i]{'start'}},@{$m[$i]{'end'}});
			my $chimera = 1;
			for (my $j = $start; $j <= $end; $j++) {
				$chimera = 0 if $l{$c}{$j};	
				}
			if ($chimera) {	
				$chimera{$c}{$m[0]{'gene'}}++;
		     
				$chimera{$c}{$m[$i]{'gene'}}++;
				for (my $j = $start; $j <= $end; $j++) {
					$l{$c}{$j}++;
					}
				}
			}
		}	
	
	#determines if these chimeras are real
	#my $db = Bio::DB::Fasta->new($dbP);
	my %newseq; my $tracker = 1;
	foreach my $c (keys %chimera) {	
		my $query = "query.fa";
		my $target = "target.fa";
		open(QUERY, ">$query");
		open(TARGET, ">$target");
		print QUERY ">$c\n$seq{$c}\n";
		foreach my $c2 (keys %{$chimera{$c}}) {
			my $seq = $pro{$c2}{'seq'};
			print TARGET ">$c2\n$seq\n";
			}
		close(QUERY); close(TARGET);					
	
		my @call = `exonerate $target $query -m protein2genome --showalignment no --showcigar 0`;
		my @s; my @e;
		if (@call) {
			my %m;
			for (my $i = 2; $i < scalar(@call) - 1; $i++) {
				my @d = split(/\s+/,$call[$i]);
				unless ($m{$d[1]}) {
					if ($d[8] eq '+') {
						push @s, $d[6];
						push @e, $d[7];
						}
					else {
						push @s, $d[7];
						push @e, $d[6];
						}
					$m{$d[1]}++;
					}
				}		
			}	
		#a blast match but no exonerate match?	
		else {
			print "Huh, this is odd. BLAST hit but no exonerate hit for $id?\n";	
			}
	
		#ensure that the two annotated parts do not overlap
		my $overlap = 0;
		for (my $i = 0; $i < scalar(@s); $i++) {
			for (my $j = 0; $j < scalar(@e); $j++) {
				unless ($i == $j) {
					if (($s[$i] < $s[$j] &&   $s[$j] < $e[$i]) || ($s[$i] < $e[$j] &&  $e[$j] < $e[$i] )) {
						$overlap++;
						}
					}
				}
			}			
		if ($overlap) {
			my $max = 0;
			my $s; my $e;
			for (my $i = 0; $i < scalar(@s); $i++) {
				if ($e[$i] - $s[$i] > $max) {
					$s = $s[$i];
					$e = $e[$i];
					}
				@s = ($s); @e = ($e);
				}			
			}

		#separate the chimeric contig out to separate contigs
		for (my $i = 0; $i < scalar(@s); $i++) { 
			my $ln = $e[$i] - $s[$i] + 1;
			$newseq{$tracker} = substr $seq{$c}, $s[$i], $ln; 
			$tracker++;
			}
		print CHE ">", $c, "\n";
                print CHE $seq{$c}, "\n";
		delete $seq{$c} if $seq{$c};	
		unlink($target); unlink($query);
		}
	close CHE;
	my $outfile = $seq . ".noChimera";
	open(OUT, ">$outfile");
	foreach my $c (keys %newseq) {
		print OUT ">contig", $c, "\n", $newseq{$c}, "\n";
		}
	foreach my $c (keys %seq) {
		print OUT ">contig", $tracker, "\n", $seq{$c}, "\n";
		$tracker++;
		}
	close(OUT);	
       
	my $call = system("cd-hit-est -i $outfile -M 0 -o temp.fa -c 1.00");
	
	open(IN, "<temp.fa");
	open(OUT2, ">temp2.fa"); 
	my $tracker2 = 1;
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>/) {
			print OUT2 ">contig", $tracker2, "\n";
			$tracker2++;
			}
		else {	
			print OUT2 $line, "\n";
			}
		}
	close(IN); close(OUT); 
	my $call2 = system("rm temp.fa*");
	my $call3 = system("mv temp2.fa $outfile");
	
	return($outfile);
	}
	
sub orfFinesse {
	my ($seq, $protein) = @_;
	my %pro =  %{$protein};
	my %seq = %$seq;
	my $frameProt = $dir . 'interruptedORFs.fa';
	
	open(OUT, ">$frameProt");
	foreach my $c (keys %seq) {
		if ($seq{$c}{'info'}) {
			my $info = $seq{$c}{'info'};
			my $gs = $1 if $info =~ m/gs(\d+)/;
			my $ge = $1 if $info =~ m/ge(\d+)/;
			my $gs0 = $gs - 1;
			my $length = $ge - $gs + 1;
		
			my $s = $seq{$c}{'seq'};
		
			my $def_orf = substr $s, $gs0, $length;
			my $pot_orf = substr $s, $gs0;
		
			my $def_aa = translate($def_orf);
			my $pot_aa = translate($pot_orf);
		
			if ($pot_aa =~ m/\*/) {
				$pot_aa = $1 if $pot_aa =~ m/^([A-Z]+)\*/;
				}
			if ($pot_aa =~ m/^\*/) {
			    $pot_aa = 1;
			}
		
			if (length($pot_aa) >  length($def_aa)) {
			    if (length($pot_aa) - length($def_aa) < 20 || ((length($pot_aa) - length($def_aa))/length($def_aa)) < 0.2) {
				my $newlength = 3 * length($pot_aa);
				my $newend = $gs + $newlength - 1;
				$info =~ s/ge$ge/ge$newend/;
			
				if ($info =~ m/_3u(\d+)/) {
					my $u3 = $1;
					my $u3new = $newend + 1;
					$info =~ s/_3u$u3/_3u$u3new/;
					}
				$seq{$c}{'info'} = $info;		
				}
			}
			elsif (length($pot_aa) < length($def_aa)) {
				#short enough to chop
				if (length($def_aa) - length($pot_aa) < 10 || (length($def_aa) - length($pot_aa))/length($def_aa) < 0.1) {
					my $newlength = 3 * length($pot_aa);
					my $newend = $gs + $newlength - 1;
					$info =~ s/ge$ge/ge$newend/;
			
					if ($info =~ m/_3u(\d+)/) {
						my $u3 = $1;
						my $u3new = $newend + 1;
						$info =~ s/_3u$u3/_3u$u3new/;
						}
					$seq{$c}{'info'} = $info;	
					}
				else {	
				    print OUT ">", $c, "\n", $seq{$c}{'seq'}, "\n";
					}
				}	
			}	
		}
	close(OUT);
	
	#my $db = Bio::DB::Fasta->new($dbP);
	#need to check ORF again
	foreach my $c (keys %seq) {
		if ($seq{$c}{'info'}) {
			my $info = $seq{$c}{'info'};
			my $gs = $1 if $info =~ m/gs(\d+)/;
			my $ge = $1 if $info =~ m/ge(\d+)/;
			my $gs0 = $gs - 1;
			my $length = $ge - $gs + 1;
		
			my $s = $seq{$c}{'seq'};
		
			my $def_orf = substr $s, $gs0, $length;
			my $pot_orf = substr $s, $gs0;
		
			my $def_aa = translate($def_orf);
			my $pot_aa = translate($pot_orf);
		
			if ($pot_aa =~ m/\*/) {
				$pot_aa = $1 if $pot_aa =~ m/^([A-Z]+)\*/;
				}
			if ($pot_aa =~ m/^\*/) {
			    $pot_aa = 1;
				}
		
			if (length($pot_aa) < length($def_aa)) {
				#short enough to chop
				if (length($def_aa) - length($pot_aa) < 10 || (length($def_aa) - length($pot_aa))/length($def_aa) < 0.1) {
					my $newlength = 3 * length($pot_aa);
					my $newend = $gs + $newlength - 1;
					$info =~ s/ge$ge/ge$newend/;
			
					if ($info =~ m/_3u(\d+)/) {
						my $u3 = $1;
						my $u3new = $newend + 1;
						$info =~ s/_3u$u3/_3u$u3new/;
						}
					$seq{$c}{'info'} = $info;	
					}
				else {	
					#re-exonerate
					my $query = "query.fa";
					my $target = "target.fa";
					print "re exonerating for contig $c\n";
					open(QUERY, ">$query");
					open(TARGET, ">$target");
					print QUERY ">$c\n$seq{$c}{'seq'}\n";
					my $seq = $pro{$seq{$c}{'match'}[0]{'gene'}}{'seq'};
					my $protlength = length($seq);
					my $proid = $seq{$c}{'match'}[0]{'gene'};
					my $contiglength = length($seq{$c}{'seq'});			
					print TARGET ">$proid\n$seq\n";
					close(QUERY); close(TARGET);					
					my @call3 = `exonerate $target $query -m protein2genome --showalignment no --showcigar 0`;
					unlink($query); unlink($target);
					my $info;
					
					my @match;
					my @prot5;
					my @prot3;
					foreach (@call3) {
						chomp (my @line = split /\s+/, $_);
						if ($line[0] =~m /vulgar/) {
							push @match, $line[6]+1;
							push @match, $line[7];
							push @prot3, $line[3];
							push @prot5, $line[2];
							}
						}
					my $start = min (@match);
					my $end = max (@match);
					my $fiveu = min (@prot5);
					my $threeu = max (@prot3);
					
					my @d = split(/\s+/,$call3[2]);
					if ($d[8] eq '+' || $d[8] eq '-') {
						if ($d[8] eq '+') {
							#this is in 5->3
							#identify gene start and stop
							$info = 'gs' . $start . '_ge' . $end;
							if ($threeu/$protlength > 0.9 || $protlength - $threeu < 11) {
								#yes, i am going to call it a  3' utr
								my $utr3 = $end+1;
								$info = $info . '_3u' . $utr3;
								}
							if ($fiveu/$protlength < 0.1 || $fiveu < 11) {
								#yes, i am going to call it a 5' utr
								my $utr5 = $start - 1;
								$info = '5u' . $utr5 . "_" . $info;
								}
							$seq{$c}{'info'} = $info;
							}		
						}
					#a blast match but no exonerate match?	
					else {
						print "Huh, this is odd. BLAST hit but not exonerate hit for $c?\n";
						delete($seq{$c}{'info'});
					}
				}
			}
		}
	}
	return(\%seq);	
}
	
	
sub translate {
	my $string = shift;
	$string = uc($string);
	my @codons = $string =~ m/(\S\S\S)/g;
	my %codons = (	'ATG'=>'M','ACG'=>'T','CTG'=>'L','CCG'=>'P','GTG'=>'V','GCG'=>'A','TTG'=>'L','TCG'=>'S',
					'ATA'=>'I','ACA'=>'T','CTA'=>'L','CCA'=>'P','GTA'=>'V','GCA'=>'A','TTA'=>'L','TCA'=>'S',
					'ATC'=>'I','ACC'=>'T','CTC'=>'L','CCC'=>'P','GTC'=>'V','GCC'=>'A','TTC'=>'F','TCC'=>'S',
					'ATT'=>'I','ACT'=>'T','CTT'=>'L','CCT'=>'P','GTT'=>'V','GCT'=>'A','TTT'=>'F','TCT'=>'S',
					'AGG'=>'R','AAG'=>'K','CGG'=>'R','CAG'=>'Q','GGG'=>'G','GAG'=>'E','TGG'=>'W','TAG'=>'*',
					'AGA'=>'R','AAA'=>'K','CGA'=>'R','CAA'=>'Q','GGA'=>'G','GAA'=>'E','TGA'=>'*','TAA'=>'*',
					'AGC'=>'S','AAC'=>'N','CGC'=>'R','CAC'=>'H','GGC'=>'G','GAC'=>'D','TGC'=>'C','TAC'=>'Y',
					'AGT'=>'S','AAT'=>'N','CGT'=>'R','CAT'=>'H','GGT'=>'G','GAT'=>'D','TGT'=>'C','TAT'=>'Y');
	my $translate;
	foreach(@codons) {
		if ($codons{$_}) {
			$translate = $translate . $codons{$_};
			}
		else {
#			print "ERROR: ILLEGAL PASS TO CODON TRANSLATION: $_ is not a codon!\n";
			$translate = $translate . 'X';
			}
		}
	return($translate);
	}