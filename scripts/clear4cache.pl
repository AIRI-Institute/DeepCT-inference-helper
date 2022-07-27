#!/usr/bin/env perl

# This script takes a VCF on STDIN
# tries to remove FILTER, QUAL, and INFO/* except for INFO/DEEPCT
# and outputs it back to STDOUT

while (<>) {
    if (not $x) {
        if (/^#/) {
            # We're in header
            if (/^#CHROM/) {
                # Last header line, remove FORMAT and samples
                s/\sFORMAT.+$//;
                $x = 1;
            }
            print;
            next;
        } else {
            # There was no header
            $x = 1
        }
    }
    # Regular VCF line
    if (m/DEEPCT/) {
        s/^(
            \S+\t # CHROM
            \S+\t # POS
            \S+\t # ID
            \S+\t # REF
            \S+\t # ALT
           )
            \S+\t # QUAL
            \S+\t # FILTER
           (\S+?;)? # Start of INFO
           (
            DEEPCT[^\t]+
           )
          (\t.+)?$  # FORMAT and samples
        /$1.\t.\t$3/x;
    } else {
        s/^(
            \S+\t # CHROM
            \S+\t # POS
            \S+\t # ID
            \S+\t # REF
            \S+\t # ALT
           )
            \S+\t # QUAL
            \S+\t # FILTER
            \S+   # INFO
           (\t.+)?$ # FORMAT and samples
        /$1.\t.\t./x;
    }
    print;
}
