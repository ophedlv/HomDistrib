#!/usr/bin/perl

=head1 NAME

pipeline.pl

=head1 AUTHOR

Ophélie Da Silva MSc student (UPSUD)

=head1 DESCRIPTION

A Perl pipeline to filter blast tabular to graph files.
Blast results are filtered by several metrics : 
- % of identity
- maximum number of mismatches
- evalue
- alignment coverage

    date : 06-2017

=head1 USAGE

    homdistrib.pl \
    --database <file.gz> \
    --fasta <file.fasta> \
    --evalue <scientific notation> (i.e. 1e-20) (default : = 1e-20) \
    --identity <integer> (default : = 80.0) \
    --cov <integer> (default : = 80.0) \
    --max_seq <integer> (default : = 0 means no cut) \
    --max_all <integer> (default : = 0 means no cut) \
    --nb_thread <integer> (default : = 20) \
    
    Note : FASTA input is needed to retrieve sequences' length

=head1 NOTE

    Only works on BLAST results formated as follow :

    -outfmt 8 (-m 8)
    
    Column identifiers :
    
    1 Query	
    2 Subject
    3 % id
    4 alignment length
    5 mistmatches
    6 gap openings
    7 q.start	
    8 q.end
    9 s.start
    10 s.end
    11 e-value
    12 bit score

    Source : http://betascience.blogspot.no/2010/01/filtering-blast-hits-by-coverage-using.html

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Bio::SeqIO;
use Bio::SearchIO;
use Cwd qw(abs_path);
use IO::Uncompress::Gunzip;


# FUNCTIONS ------------------------------------------------------------
sub db_creation {
	my ($dir_db,$db_path,$db_file, $db_name) = @_;
	
	unless(-e $dir_db) {
		mkdir($dir_db) or warn ("Error : impossible to create db folder \n"); 
	}

	unless(-e $db_path && glob("$db_path/*")) { # if specific db folder not exists or empty
		unless(-e $db_path) {
			mkdir($db_path) or warn ("Error : impossible to create $db_name folder \n"); 
		}
		system("zcat $db_file > $db_path/tmp");
		system("makeblastdb -in $db_path/tmp -dbtype nucl -out $db_path/$db_name");
		system("rm $db_path/tmp");
	}
}



sub blast {
	my ($dir_bl, $db_path, $db_name, $seq_file, $seq_name, $nb_thread) = @_;
	
	unless(-e $dir_bl) {
		mkdir($dir_bl) or warn ("Error : impossible to create bl folder \n"); 
	}

	unless(-e "$dir_bl$seq_name\_tblastn_res.bl") {
		system ("tblastn -query $seq_file -db $db_path/$db_name -num_threads $nb_thread -outfmt '7 qseqid sseqid qlen slen length pident evalue' -out $dir_bl$seq_name\_tblastn_res.bl");
	}
}

sub filtering {
	my ($dir_res,$seq_name,$dir_bl,$evalue_opt,$identity_opt,$cov_opt) = @_;
	
	unless(-e $dir_res) {
		mkdir($dir_res) or warn ("Error : impossible to create res folder \n"); 
	}

	unless(-e "$dir_res$seq_name") {
		mkdir($dir_res.$seq_name) or warn ("Error : impossible to create res/$seq_name folder \n");
	}

	open IN, "$dir_bl$seq_name\_tblastn_res.bl";
		my $min;
		my $cov;
		my $qcov;
		my $scov;
		
		my $resultats_detailles = "#query\ttarget\tevalue\tidentity\talign\tqcov\tscov\n";
		my $resultats_synth = "#query\ttarget\n";
		while(<IN>) {
			unless (/^#/) {
				chomp;
				my ($qid,$sid,$qlen,$slen,$alen,$pid,$evalue)=split(/\t/,$_);
				$evalue =~ s/e/1e/ if ($evalue =~ /^e/);
				
				if ($qlen < $slen) {
					$min = $qlen;
				}
				else {
					$min = $slen;
				}
				$cov = $cov_opt*$min/100;
				$qcov = $alen/$qlen*100;
				$scov = $alen/$slen*100;
				
				if (($evalue <=  $evalue_opt) && (($pid >= $identity_opt) && ($alen >= $cov))) {
					$resultats_detailles.="$qid\t$sid\t$evalue\t$pid\t$alen\t$qcov\t$scov\n";
					$resultats_synth.="$qid\t$sid\n";
				}
			}
		}
	close IN;
	
	open (RES, ">$dir_res$seq_name/$seq_name\_$evalue_opt\_id$identity_opt\_cov$cov_opt\_res_detail.txt") || die ("Vous ne pouvez pas créer le fichier");
	print RES $resultats_detailles;
	close RES;
	
	open (RES, ">$dir_res$seq_name/$seq_name\_$evalue_opt\_id$identity_opt\_cov$cov_opt\_res_synth.txt") || die ("Vous ne pouvez pas créer le fichier");
	print RES $resultats_synth;
	close RES;
}


sub selection {
	
	my ($dir_res,$seq_name,$evalue_opt,$identity_opt,$cov_opt,$max_seq_opt,$max_all_opt) = @_;
	
	if($max_seq_opt > 0) {
		my %seq = ();
		my $res_all = "#query\ttarget\tevalue\tid\talen\cov\n";
		my $res_synth = "#query\ttarget\n";
			
		open(IN, "$dir_res$seq_name/$seq_name\_$evalue_opt\_id$identity_opt\_cov$cov_opt\_res_detail.txt");
		while(<IN>) {
			unless(/^#/) {
				chomp;
				my ($query,$target,$evalue,$id,$alen,$cov)=split(/\t/,$_);
				if(exists $seq{$query}) {
					$seq{$query}+=1;
				}
				else {
					$seq{$query} = 1;
				}
				
				if($seq{$query} <= $max_seq_opt) {
					$res_all .= "$query\t$target\t$evalue\t$id\t$alen\t$cov\n";
					$res_synth .= "$query\t$target\n";
				}
			}
		}
		close(IN);
		
		open (RES, ">$dir_res$seq_name/$seq_name\_$evalue_opt\_id$identity_opt\_cov$cov_opt\_res$max_seq_opt\_detail.txt") || die ("Vous ne pouvez pas créer le fichier");
		print RES $res_all;
		close RES;
		open (RES, ">$dir_res$seq_name/$seq_name\_$evalue_opt\_id$identity_opt\_cov$cov_opt\_res$max_seq_opt\_synth.txt") || die ("Vous ne pouvez pas créer le fichier");
		print RES $res_synth;
		close RES;
	}
		
	if($max_all_opt > 0) {
		my @all = ();	
		my $res_all = "#query\ttarget\tevalue\tid\talen\tcov\n";
		my $res_all_synth = "#query\ttarget\n";
		
		open(RES, "$dir_res$seq_name/$seq_name\_$evalue_opt\_id$identity_opt\_cov$cov_opt\_res_detail.txt") or die ("Can't open file");
		my $i=0;
		while(<RES>) {
			unless(/^#/) {
				chomp;
				my @tmp=split(/\t/,$_); # la colonne 2 contient la e-value
				$all[$i][0] = $tmp[2];
				$all[$i][1] = $_;
				$i+=1;
			}
		}
		@all = sort {$a->[0] <=> $b->[0]} @all;
		
		for ($i = 0; $i < $max_all_opt; $i++) {
			my ($query,$target,$evalue,$id,$alen,$cov) = split(/\t/,$all[$i][1]);
			
			$res_all .= "$query\t$target\t$evalue\t$id\t$alen\t$cov\n";
			$res_all_synth .= "$query\t$target\n";
		}
		close(RES);
		
		open (RES, ">$dir_res$seq_name/$seq_name\_$evalue_opt\_id$identity_opt\_cov$cov_opt\_res_all$max_all_opt\_detail.txt") || die ("Vous ne pouvez pas créer le fichier");
		print RES $res_all;
		close RES;
		open (RES, ">$dir_res$seq_name/$seq_name\_$evalue_opt\_id$identity_opt\_cov$cov_opt\_res_all$max_all_opt\_synth.txt") || die ("Vous ne pouvez pas créer le fichier");
		print RES $res_all_synth;
		close RES;
	}
}


sub abund_matrix {
	my ($dir_res,$dir_data,$seq_name,$evalue_opt,$identity_opt,$cov_opt,$max_seq_opt,$max_all_opt) = @_;	
		# Analysis (results -> genes)
	open (IN, "$dir_res$seq_name/$seq_name\_$evalue_opt\_id$identity_opt\_cov$cov_opt\_res\_detail.txt");
	my %genes;
	my $homologue;
	while(<IN>) {
		unless(/^#/) {
			$homologue = (split(/\t/,$_))[1];
			if (exists $genes{$homologue}) {
				$genes{$homologue} += 1;
			}
			else {
				$genes{$homologue} = 1;
			}
		}
	}
	close(IN);

		# MetaG matrix (genes -> sampleID)
	my %sample;
	my @tmp;
	my @abund;
	my @sum;
	my $sub_matrix = "";
	my $n;

	my $IN = IO::Uncompress::Gunzip->new("$dir_data/metaG.gz") or die "IO::Uncompress::Gunzip failed: Error\n";
	while (<$IN>) {
		if(/^#/) {
			my $str = $_; 
			$str =~ s/^#//;
			$sub_matrix .= $str;
			@tmp = split(/\t/);
			$n = @tmp;
			@sum = 0 x $n;
		}
		else {
			@abund = split(/\t/);
			if (exists $genes{$abund[0]}) {
				$sub_matrix .= $_;
				for (my $i = 1; $i < $n; $i++) {
					unless($abund[$i] == 0) {
						$sum[$i]+=$abund[$i];
					}
				}
				$sub_matrix .= "\n";
			}
		}
	}
	
		# Abundance sub-matrix saving 
	open(RES, ">$dir_res$seq_name/$seq_name\_$evalue_opt\_id$identity_opt\_cov$cov_opt\_sub_matrix_abund.txt");
	print RES $sub_matrix;
	close(RES);
}




# MAIN -----------------------------------------------------------------

# Retrieving arguments
my $h = 0;
my $help = 0;

# Default option values
my $evalue_opt = 1e-20;
my $identity_opt = 80.0;
my $cov_opt = 80.0;
my $max_seq_opt = 0;
my $max_all_opt = 0;
my $nb_thread = 20;

# File variable
my $db_file;
my $seq_file;



GetOptions('help|?' => \$help, 
            man => \$h,
            'database=s' => \$db_file,
            'fasta=s' => \$seq_file,
            'evalue=f' => \$evalue_opt,
            'identity=i' => \$identity_opt,
            'cov=f' => \$cov_opt,
            'max_seq=i' => \$max_seq_opt,
            'max_all=i' => \$max_all_opt,
            'nb_thread=i' => \$nb_thread) or pod2usage(2);
            
pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $h;


# Directories variables
my $dir = abs_path($0);
my @tmp = ((split(/\//,$dir))[-1],(split(/\//,$dir))[-2]);
$dir =~ s/$tmp[0]//;
$dir =~ s/$tmp[1]\///;

my $dir_db = $dir."db/";
my $dir_bl = $dir."bl/";
my $dir_res = $dir."res/";
my $dir_data = $dir."data/";

my $seq_name = $seq_file;
$seq_name =~ s/(.*)\///;
$seq_name =~ s/\.fasta//;

my $db_name = (split(/\./,(split(/\//,$db_file))[-1]))[0]; # nom de la base de données
my $db_path = $dir_db.$db_name;


# DataBase creation
db_creation($dir_db,$db_path,$db_file, $db_name);

# Blast
blast($dir_bl, $db_path, $db_name, $seq_file, $seq_name, $nb_thread);

# Data Filtering
filtering($dir_res,$seq_name,$dir_bl,$evalue_opt,$identity_opt,$cov_opt);

# Data Selection
selection($dir_res,$seq_name,$evalue_opt,$identity_opt,$cov_opt,$max_seq_opt,$max_all_opt);

# Sub-matrix of abundance creation
abund_matrix($dir_res,$dir_data,$seq_name,$evalue_opt,$identity_opt,$cov_opt,$max_seq_opt,$max_all_opt);

# Plot
system("Rscript ".$dir."pipeline/src/plot_genes.R $dir_res$seq_name/$seq_name\_$evalue_opt\_id$identity_opt\_cov$cov_opt\_sub_matrix_abund.txt")
