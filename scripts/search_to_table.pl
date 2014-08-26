#!env perl
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

my $gff_dir = 'GFF';
my $results_dir = 'results';
my $gff_ext = 'gff3';
my $ref ='GAGcluster.fas';
GetOptions(
    'g|gff:s'     => \$gff_dir,
    'r|results:s' => \$results_dir,
    'ref:s'       => \$ref,
    );

open(my $r => $ref ) || die "cannot open ref file $ref: $!";
my %refs;
my $order = 0;
while(<$r>) {
    if( /^>(\S+)/ ) {
	$refs{$1}->{order} = $order++;
    }
}

print join("\t", qw(CLUSTER_GENE TARGET_GENE PERCENT_ID EVALUE 
TARGET_CHR TARGET_START TARGET_END TARGET_STRAND)), "\n";

opendir(my $resdir => $results_dir) || die "cannot open $results_dir: $!";
for my $file ( readdir($resdir) ) {
    next unless $file =~ /-vs-(\S+)\.(FASTA|BLAST).tab/;
    my $sp = $1;
    my $tab;
    if( $file =~ /\.gz/ ) {
	open($tab => "zcat $results_dir/$file |" ) || die "cannot open $results_dir/$file: $!";
    } else {
	open($tab => "$results_dir/$file" ) || die "cannot open $results_dir/$file: $!";
    }
    my %cluster;
    while(<$tab>) {
	my ($query,$target,$pid,@row) = split;
	my $score = pop @row;
	my $evalue = pop @row;
	push @{$cluster{$query}}, [$target, $pid, $evalue, @row ];
    }
    #warn("$sp\n");
    my $gff = File::Spec->catfile($gff_dir,"$sp.$gff_ext");
    my $cmd;
    if( ! -f $gff ) {
	my $first = $gff;
	if( ! -f $gff.".gz" ) {
	    warn("cannot find $gff or $gff.gz for $sp\n");
	    next;
	}  
	$gff .= ".gz";	
	$cmd = "zcat $gff |";
    } else {
	$cmd = "< $gff";
    }
    
    open(my $fh => $cmd) || die "cannot open $gff: $!";
    my %genes;
    while(<$fh>) {
	next if /^\#/;
	chomp;
	my @row = split(/\t/,$_);
	if( lc($row[2]) eq 'mrna' ) {
	    my %lastcol;
	    eval { 
		%lastcol = map { my ($k,$v) = split(/=/,$_);
		($k,$v) } split(/;/,pop @row);
	    };
	    if( $@ ) {
		warn($@,"\n");
	    }
	    if( my $id = ($lastcol{'Name'} || $lastcol{'ID'} ) ) {
		$id =~ s/\.tr//;		
		if( $sp =~ /Bgran3/ ) {
		    $id =~ s/\.\d+$//;
		}
		if( $id =~ /\|/ ) {
		    my ($p,$stripid) = split(/\|/,$id);
		    if( $stripid ) {
			$id = $stripid;
		    } else {
			warn("no second value for $id ($sp)\n");
		    }
		}
	       
		$genes{$id} = [ $row[0], $row[3], $row[4], $row[6] ];
	    } else {
		warn("cannot find an ID for ",
		     join(";", map { sprintf("%s=%s",$_, 
					     $lastcol{$_}) } keys %lastcol),
		     "\n");
	    }
	}
    }
    print "#$sp\t", (scalar keys %genes), " genes\n";
    for my $gene ( sort { $refs{$a}->{order} <=> $refs{$b}->{order} } 
		   keys %refs ) {
	
	for my $n ( @{$cluster{$gene}} ) {
	    my ($spp_n,$hit_name) = split(/\|/,$n->[0]);
	    print join("\t", 
		       $gene,$n->[0],$n->[1],$n->[2],
		       @{$genes{$hit_name} || []}), "\n";
	    
	}
    }
}
