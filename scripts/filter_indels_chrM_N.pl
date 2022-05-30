#!/usr/bin/env perl

# This script takes VCF with header on the STDIN and filters out:
# - chrM
# - indels (Selene can't parse them currently)
# - non-AGTC alleles (probably these are just Ns after conversion to hg19)
# Filtered stuff goes to STDERR, good stuff to STDOUT, both headered

while (<>) {
    if (not $x) {
    # header
        if (/^#/) {
            print $_;
            print STDERR $_;
            next
        } else {
            $x = 1
        }
    }
    @f = split;
    if ($f[0] eq "chrM" or
        length($f[3]) > 1 or length($f[4]) > 1 or 
        $f[3] !~ /^[AGTC]$/ or $f[4] !~ /^[AGTC]$/) {
        print STDERR $_
    } else {
        print $_
    }
}
