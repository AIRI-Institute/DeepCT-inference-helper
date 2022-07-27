#!/usr/bin/env perl

# This script just takes VCF without header from STDIN
# and puts same VCF with REF and ALT swapped to STDOUT

while (<>) {
    s/\tFail\(REF\=\=ALT\)//;
    @f = split;
    print join("\t", @f[0..2], @f[4,3], @f[5..9])."\n"
}
