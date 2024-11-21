#!/bin/awk -f
BEGIN {
    OFS='\t'
    FS='\t'
}
{
    if (/^[^@]/) {
        gsub("/2","",$1); gsub("/1","",$1); print 
    }
    else {
        print $0
    }
}
