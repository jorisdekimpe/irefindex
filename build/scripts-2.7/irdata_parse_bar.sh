#!/bin/sh

# Copyright (C) 2012 Ian Donaldson <ian.donaldson@biotek.uio.no>
# Original author: Paul Boddie <paul.boddie@biotek.uio.no>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

if [ -e "irdata-config" ]; then
    . "$PWD/irdata-config"
elif [ -e "scripts/irdata-config" ]; then
    . 'scripts/irdata-config'
else
    . 'irdata-config'
fi

if [ "$1" = '--help' ]; then
    cat 1>&2 <<EOF
Usage: $PROGNAME <output data directory> <filename>...

Process the given MITAB files, producing import data in the given output
directory.
EOF
    exit 1
fi

DATADIR=$1
shift 1
FILENAMES=$*

if [ ! "$DATADIR" ] || [ ! "$FILENAMES" ]; then
    echo "$PROGNAME: A data directory and input filenames must be specified." 1>&2
    exit 1
fi


#it is assumed that all records in the virushost file refer to protein-protein interactions
#it is not possible to check for these (as done in the below commented code) because virushost
#does not follow regular mitab conventions and interactor type is not provided
#should this change in future, the following code could be uncommented
#
#preprocess input files to remove non protein-protein interaction records
#save a copy of the original
for THISFILE in $FILENAMES; do
    ORIGINAL="$THISFILE".original
   ## awk 'BEGIN {FS="\t";}{if ($21 == "psi-mi:\"MI:0326\"(protein)"  && $22 == "psi-mi:\"MI:0326\"(protein)") print $0;}' < "$THISFILE" > tmp
    awk 'BEGIN {FS="\t";}{if ($12 == "psi-mi:\"MI:1047\"(protein protein)") print $0;}' < "$THISFILE" > tmp

    #catch errors - $? means return value of last function
    if [ $? != 0 ]; then
        echo "$PROGNAME: pre-processing of $THISFILE failed." 1>&2
        exit 1
    fi
    mv $THISFILE $ORIGINAL
    mv tmp $THISFILE  
done

"$TOOLS/irdata_parse_mitab.py" 'BAR' "$DATADIR" $FILENAMES
