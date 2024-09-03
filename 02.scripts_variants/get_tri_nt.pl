#!/usr/bin/perl

use strict;
use warnings;

# Read the reference genome sequence
my $genome_file = '../00.data/B_uniformis.fasta';
#my $genome_file = shift @ARGV;
my $reference_genome_sequence = read_fasta($genome_file);

# Open the input mutation file
my $mutation_file = shift @ARGV;
#my $mutation_file = 'mutation_table.txt';
open my $fh, '<', $mutation_file or die "Couldn't open $mutation_file: $!";

# Print the updated header
print "Position\tREF\tTrinucleotide\tRC_trinucleotide\n";

# Process each line of the input mutation file
while (my $line = <$fh>) {
    chomp $line;
    my ($position, $ref, $chromosome) = split("\t", $line);

    my $trinucleotide = substr($reference_genome_sequence, $position - 2, 3);

    my $revcomp = reverse $trinucleotide;
    $revcomp =~ tr/ATGCatgc/TACGtacg/;

    print "$position\t$ref\t$trinucleotide\t$revcomp\n";
}

# Close the input file
close $fh;

# Subroutine to read a FASTA file into a hash
sub read_fasta {
    my ($file) = @_;
    open my $fh, '<', $file or die "Couldn't open $file: $!";
    my $sequence = '';
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^>/;  # Skip header lines
        $sequence .= $line;
    }
    close $fh;
    return $sequence;
}
