#! /bin/bash

./abb2gd.pl evidence.txt evidence.gd
gdtools MUTATIONS -r REL606.6.GENES.1400-flanking.gff3 -r REL606.6.GENES.MASKED.gff3 -o mutations.gd evidence.gd
./gd2abb.pl mutations.gd mutations.txt
