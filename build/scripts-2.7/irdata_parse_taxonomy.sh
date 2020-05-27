#!/bin/sh

# Copyright (C) 2011 Ian Donaldson <ian.donaldson@biotek.uio.no>
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

OUTFILE='names.txt'

if [ "$1" = '--help' ]; then
    cat 1>&2 <<EOF
Usage: $PROGNAME <output data directory> <filename>

Process the Taxonomy names file (typically names.dmp), producing data suitable
for iRefIndex in a file called $OUTFILE in the output data directory.
EOF
    exit 1
fi

DATADIR=$1
FILENAME=$2

if [ ! "$DATADIR" ] || [ ! $FILENAME ]; then
    echo "$PROGNAME: A data directory and an input filename must be specified." 1>&2
    exit 1
fi

# Convert the bizarre format to a plain tab-separated format.
# Tab character used in the sed command.
TAB=`printf '\t'`

sed "s/${TAB}|${TAB}/${TAB}/g;s/${TAB}|$//;" "$FILENAME" > "$DATADIR/tmp1.txt"
  
# In mny cases, column 3 is NULL 
# In many cases, there may be two or more rows that have identical values
# for columns 1, 2 and 4 and only one row will have a non-NULL value in 
# column 3.
# This causes a problem for the import stage which enforces a unique compound key
# for the taxonomy_names table composed of columns 1,2 and 3.
# So, remove column 3 entirely (just leaving a TAB) and unique the resulting table
# such that each row will have a unique compound key

# replace column 3 with an empty value
awk 'BEGIN{FS="\t"; OFS="";}{ print $1,"\t",$2,"\t","\t",$4}' < $DATADIR/tmp1.txt > $DATADIR/tmp2.txt
# ensure that each row is unique
sort -u $DATADIR/tmp2.txt > $DATADIR/$OUTFILE
#report
echo "$PROGNAME: Redundant entries in $OUTFILE were removed." 1>&2
echo "$PROGNAME: Lines in original $OUTFILE." 1>&2
#wc -l $DATADIR/tmp2.txt
echo "$PROGNAME: Lines in final $OUTFILE." 1>&2
#wc -l $DATADIR/$OUTFILE

exit $?
