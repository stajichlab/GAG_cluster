#!env perl
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

my $gff_dir = 'GFF';
my $results_dir = 'results';
my $gff_ext = 'gff3';
GetOptions(
    'g|gff:s' => \$gff_dir,
    'r|results:s' => \$results_dir,
    );

opendir(my $resdir => $results_dir) || die "cannot open $results_dir: $!";
for my $file ( readdir($resdir) ) {
    next unless $file =~ /-vs-(\S+)\.(FASTA|BLAST).tab/;
    my $sp = $1;

    #warn("$sp\n");
    my $gff = File::Spec->catfile($gff_dir,"$sp.$gff_ext");
    if( ! -f $gff && ! -f $gff.".gz") {
	warn("cannot find $gff for $sp\n");
	next;
   }
    
}
