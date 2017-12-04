#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use YAML::XS qw/LoadFile/;
use autodie;

my ($fn_in, $fn_out) = @ARGV;

die "Can't find or read input file: $!\n"
    if (! -r $fn_in);

# set output filehandle based on arguments
my $fh = \*STDOUT;
if (defined $fn_out) {
    open $fh, '>', $fn_out;
}

my $yaml = LoadFile($ARGV[0]);

convert($yaml);

sub convert {

    my ($yaml) = @_;

    print {$fh} header();

    say {$fh} "    <h3>Summary statistics</h3>";

    for my $grp (sort keys %$yaml) {

        my $ref = $yaml->{$grp};

        next if (! ref $ref);
        next if (! defined $ref->{'total.gigabases'});
        
        print {$fh} <<"CONTENT"
    <table>
        <caption>$grp</caption>
        <tr>
            <td>Total Yield (Gb)</td>
            <td>$ref->{'total.gigabases'}</td>
        </tr>
        <tr>
            <td>Total Reads</td>
            <td>$ref->{'total.reads'}</td>
        </tr>
        <tr>
            <td>Mean Length</td>
            <td>$ref->{'mean.length'}</td>
        </tr>
        <tr>
            <td>Median Length</td>
            <td>$ref->{'median.length'}</td>
        </tr>
        <tr>
            <td>Max Length</td>
            <td>$ref->{'max.length'}</td>
        </tr>
        <tr>
            <td>Mean Q</td>
            <td>$ref->{'mean.q'}</td>
        </tr>
        <tr>
            <td>Median Q</td>
            <td>$ref->{'median.q'}</td>
        </tr>
    </table>
CONTENT

    }

    my %figs = (
        'length_histogram'  => "Read length distribution",
        'q_histogram'       => "Mean quality score distribution",
        'reads_per_hour'    => "Yield over time",
        'cumulative_yield'  => "Cumulative yield over time",
        'yield_summary'     => "Yield by read length cutoff",
        'flowcell_overview' => "Median read quality per channel",
        'length_by_hour'    => "Read length over time",
        'q_by_hour'         => "Read quality over time",
        'length_vs_q'       => "Read length vs. quality",
    );

    my @order = qw/
        length_histogram
        q_histogram
        reads_per_hour
        cumulative_yield
        yield_summary
        flowcell_overview
        length_by_hour
        q_by_hour
        length_vs_q
    /;


    say {$fh} "    <h3>QC plots</h3>";
    say {$fh} "    <p>(Click on plot for hi-resolution version)</p>";

    for my $base (@order) {

        my $caption = $figs{$base} // die "No caption found for $base";

        print {$fh} <<"CONTENT"
    <a href="$base.png">
        <figure>
            <img src="$base.screen.png" alt="$base" />
            <figcaption>$caption</figcaption>
        </figure>
    </a>
CONTENT

    }

    print {$fh} footer();

}



sub header {

    return <<'HEADER';
<?xml version="1.0" encoding="utf-8"?>
<!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
    <title></title>
    <meta http-equiv="content-type" content="application/xhtml+xml; charset=utf-8" />
    <style>
        h2 {
            padding: 0.3em;
            background-color: #000000;
            color: #ffffff;
            margin: 1em 0 2em 0;
        }
        h3 {
            padding: 0em 0.2em 0em 0.2em;
            color: #555555;
            border: solid 1px black;
            border-width: 0px 0px 1px 0px;
            margin: 2em 0 0.4em 0;
        }
        tr {
            margin: 0;
        }
        tr:nth-child(even) {
            background-color: #bbbbbb;
        }
        tr:nth-child(odd) {
            background-color: #eeeeee;
        }
        caption {
            text-align: left;
            font-weight: bold;
            background-color: #550000;
            color: #ffffff;
            padding: 0.1em 0.2em;
        }
        table {
            margin: 1em;
            padding: 0.3em;
        }
        td {
            margin: 0;
            padding: 0 0.4em;
        }
        tr td:nth-child(1) {
            color: #550000;
        }
        figure {
            display: table;
            margin: 2em 0;
        }
        figcaption {
            display: table-caption;
            caption-side: top;
            font-size: 1.1em;
            text-decoration: none;
            text-align: center;
            font-weight: bold;
            background-color: #550000;
            color: #ffffff;
            padding: 0.1em 0.2em;
            margin: 2em 0 0.7em 0;
        }
            
    </style>
</head>
<body>

    <h2>NanoporeQC Report</h2>
HEADER

}

sub footer {

    return <<'FOOTER';

</body>
</html>
FOOTER

}
