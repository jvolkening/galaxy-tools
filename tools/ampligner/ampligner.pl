#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use BioX::Seq::Stream;
use Getopt::Long;
use File::Temp;

my $fn_in;
my $fn_out;
my $tgt_depth = 100;
my $min_len   = 1;

our $VERSION = 0.001;

GetOptions(
    'in=s'         => \$fn_in,
    'out=s'        => \$fn_out,
    'depth=i'      => \$tgt_depth,
    'min_length=i' => \$min_len,
    'version'      => sub{ say $VERSION; exit; },
);

die "Failed to find or read input file\n"
    if (! -r $fn_in);

die "Output path failed taint check\n"
    if ($fn_out =~ /( \.\. | \| | ; | > | \& )/);

my $p = BioX::Seq::Stream->new($fn_in);
my $n_seq;
while (my $seq = $p->next_seq) {
    next if (length($seq) < $min_len);
    ++$n_seq;
}

my $aligned = File::Temp->new(CLEANUP => 1);

open my $mafft, '|-',
    "mafft --localpair --maxiterate 1000 --adjustdirection --op 0.1 --lop -0.1 --quiet - > $aligned";


$p = BioX::Seq::Stream->new($fn_in);

while (my $seq = $p->next_seq) {
    next if (length($seq) < $min_len);
    next if (rand() > $tgt_depth/$n_seq);
    print {$mafft} $seq->as_fasta;
}

close $mafft;

# calculate consensus

my @seqs;
$p = BioX::Seq::Stream->new($aligned);

while (my $s = $p->next_seq) {
    push @seqs, $s;
}

my $l = length($seqs[0]);

my $cons;

for my $i (0..$l-1) {

    my @b = map {substr($_, $i, 1)} @seqs;
    my %counts;
    for (@b) {
        ++$counts{$_};
    }
    my @srt = sort {$counts{$b} <=> $counts{$a}} keys %counts;

    my $c = scalar(@srt) == 1 ? $srt[0]
          : $counts{$srt[0]} == $counts{$srt[1]} ? 'N'
          : $srt[0];

    $c = '' if ($c eq '-');
    $cons .= uc $c;

}

open my $out, '>', $fn_out;
my $seq = BioX::Seq->new( $cons, 'consensus' );
print {$out} $seq->as_fasta;

exit;
