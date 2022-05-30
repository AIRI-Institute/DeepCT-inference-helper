#!/usr/bin/env perl

# This script takes a VCF on STDIN
# and outputs a TSV for Selene to STDOUT

while (<>) {
    @f = split;
    print join("\t", @f[0,1,3,4])."\n"
}
