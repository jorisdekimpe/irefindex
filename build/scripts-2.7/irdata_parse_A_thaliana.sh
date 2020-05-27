#!/bin/sh

# Copyright (C) 2011, 2012 Ian Donaldson <ian.donaldson@biotek.uio.no>
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
Usage: $PROGNAME <output data directory> <filename>

Process the UniProt FlyBase or Yeast file (typically fly.txt or yeast.txt),
producing data suitable for iRefIndex in a file with the same name in the output
data directory.
EOF
    exit 1
fi

DATADIR=$1
FILENAME=$2
echo "$FILENAME" 1>&2
if [ ! "$DATADIR" ] || [ ! "$FILENAME" ]; then
    echo "$PROGNAME: A data directory and an input filename or - must be specified." 1>&2
    exit 1
fi

OUTFILE=`basename "$FILENAME"`
FILETYPE=${OUTFILE%%.*}
echo "$FILETYPE"
# Find the last table header marker.
#this is taken from fly but not necessart since there is not comment, it is just header
#followed by data
#START=`grep -ne '^__' "$FILENAME" | tail -n 1 | cut -d ':' -f 1`

# Find the footer.

#END=`grep -ne '^----' "$FILENAME" | tail -n 2 | head -n 1 | cut -d ':' -f 1`

# cat the file and present it to the parser.
echo "$DATADIR/$OUTFILE" 1>&2 
  cat "$FILENAME" \
| "$TOOLS/irdata_parse_A_thaliana2.py" "$FILETYPE" --discard-ill-formed \
> "$DATADIR/$OUTFILE"
