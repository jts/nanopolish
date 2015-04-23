#! /usr/bin/env perl
#
# Extract fast5 file names from a poretools fasta file and prepend an index to uniquify names
#
use strict;
use Bio::Perl;
use Getopt::Long;

my $inFile  = Bio::SeqIO->new(-file => $ARGV[0], '-format' => 'Fasta');
my $outFile = Bio::SeqIO->new(-fh => \*STDOUT, '-format' => 'Fasta');
open(FM, ">$ARGV[0].fast5.fofn");

my $id = 0;

while(my $seq = $inFile->next_seq())
{
    my $in_name = $seq->display_id();
    my $out_name = "$id.$in_name";
    my $desc = $seq->desc();
    my ($path) = $desc =~ /([^:\s]+\/\S+)/;
    if (not defined $path) { $path = $desc; }
    
    print FM "$out_name\t$path\n";

    $seq->display_id($out_name);
    $outFile->write_seq($seq);
    $id++;
}
