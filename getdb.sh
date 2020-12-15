#!/bin/bash

echo "This script needs wget and diamond"
echo "Downloading databases from http://klif.uu.nl/download/plasmid_db/"
echo "Please be patient. This is a large download"
echo
echo "Download plasmid genes databases"
wget -O plasmiddb_cge.faa http://klif.uu.nl/download/plasmid_db/plasmiddb_cge.faa
wget -O plasmiddb_total.faa http://klif.uu.nl/download/plasmid_db/plasmiddb_total.faa
echo "Indexing files"
diamond makedb --in plasmiddb_cge.faa -d plasmiddb_cge
diamond makedb --in plasmiddb_total.faa -d plasmiddb_total
echo "Done"
echo ""
exit 1

