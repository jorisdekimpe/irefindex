#!/bin/sh

# Copyright (C) 2011, 2012 Ian Donaldson <ian.donaldson@biotek.uio.no>
# Copyright (C) 2013 Paul Boddie <paul@boddie.org.uk>
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

Process the RefSeq files (typically in the "feature table" text format),
producing data suitable for iRefIndex in the output data directory.
EOF
    exit 1
fi

DATADIR=$1

if [ ! "$DATADIR" ]; then
    echo "$PROGNAME: A data directory must be specified." 1>&2
    exit 1
fi

shift 1
FILENAMES=$*

if [ ! "$FILENAMES" ]; then
    echo "$PROGNAME: Input filenames must be specified." 1>&2
    exit 1
fi

# Parse the data files.

  "$SCRIPTS/argument-per-line" $FILENAMES \
| "$SCRIPTS/irparallel" "\"$TOOLS/irdata_parse_refseq.py\" \"$DATADIR\" {}"

# Concatenate the output data.
echo "$DATADIR cat step" 1>&2 
cat "$DATADIR"/*_proteins > "$DATADIR/refseq_proteins.txt"
cat "$DATADIR"/*_identifiers > "$DATADIR/refseq_identifiers.txt"

# Process the sequence data.

"$TOOLS/irdata_process_signatures.sh" "$DATADIR" --append --append-length
exit $?
