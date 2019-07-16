#!/usr/bin/awk -f

BEGIN {
    OFS="\t"
}

$3 == "cds" {
    $2 = "p2g"
    $3 = "CDSpart"
    $9 = "source=P"
    print $0
}
