#!/usr/bin/python

"""
Parse FASTA format text files, writing protein flat files to the output data
diimport reacrectory.

--------

Copyright (C) 2012 Ian Donaldson <ian.donaldson@biotek.uio.no>
Original author: Paul Boddie <paul.boddie@biotek.uio.no>

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import re

try:
    set
except NameError:
    from sets import Set as set

class Parser:

    "A parser for FASTA format text files."
    import re
    NULL = r"\N"

    def __init__(self, source, f, f_main, identifier_types, output_identifier_types):

        """
        Initialise the parser with the given input file object 'f', an output
        file object 'f_main', a list of 'identifier_types' to be extracted, and
        a list of 'output_identifier_types' to be produced (potentially
        including types that are not provided in the data and will therefore be
        null in the output).
        """
        self.source = source
        self.f = f
        self.f_main = f_main
        self.identifier_types = identifier_types
        self.output_identifier_types = output_identifier_types
        self.lineno = None

    def close(self):
        if self.f is not None:
            self.f.close()
            self.f = None
        if self.f_main is not None:
            self.f_main.close()
            self.f_main = None

    def parse_header(self, line):
        records = []

        # Records are separated by ^A (ASCII) bytes.

        for header_record in line.rstrip("\n").lstrip(">").split("\x01"):
            identifiers = {}
            header_elements = []

            # Record data is separated from descriptions by white-space.

            # depending on source e.g. Uniprot different iterables
            if self.source in ["PDB"]:
                # Values are separated by space character and chain is split by _
                header_elements = header_record.split(" ")[0].split("_")
            elif self.source in ["UNIPROT", "IPI"]:
                # Values are separated by pipe/bar characters.
                header_elements = header_record.split("|")
            elif self.source in ["GENPEPT"]:
                # Values are formatted >ACCESSION.VERSION text with spaces [organism]
                # >ERL63622.1 hypothetical protein L248_2681 [Lactobacillus shenzhenensis LY-73]
                #TO DO:If there is no organism known, the line and sequence will be skipped
                #
                bracketmatch = re.search(r'^(.+\.[0-9])\s(.+)\(\[(.+)\]$', header_record)
                if bracketmatch:
                    bracketmatch = re.split(r'^(.+\.[0-9])\s(.+)\(\[(.+)\]$', header_record)
                    if len(bracketmatch) == 5:
                        header_elements = bracketmatch[1:-1]
                    else:
                        raise ValueError, "Header %s is not formatted as accession description [organism]  and should be corrected." % (header_record)
                match = re.search(r'^(.+\.[0-9])\s(.+)\s\[(.+)\]$', header_record)
                if match:
                    match = re.split(r'^(.+\.[0-9])\s(.+)\s\[(.+)\]$', header_record)
                    if len(match) == 5:
                        header_elements = match[1:-1]
                    else:
                        raise ValueError, "Header %s is not formatted as accession description [organism]  and should be corrected." % (header_record)
                #print >> sys.stderr, "brmatch: %s" % (bracketmatch)
                #print >> sys.stderr, "match: %s" % (match)
                # Note: raise will stop the run!
                #patent = re.search(r'^(.+)\.([0-9])\s(.+)[0-9]$', header_record)
                if bracketmatch is None and match is None:
                    header_elements = None 
                    print >>sys.stderr, "Header %s is not formatted correctly (no organism)  and should be corrected" % (header_record)
                

                # Values are formatted >ACCESSION.VERSION text with spaces([organism]
                # Values are formatted >ACCESSION.VERSION very long text with spaces
                # only get the description without organism 
                #if (False and  match or matchall):
                #    matchall = re.split(r'^(.+)\.([0-9])\s(.+)$', header_record)
                #    header_elements = matchall[-1]
                #else:
                #    raise ValueError, "Header %s is not formatted as expect and should be corrected." % (header_record) 
            else:
                raise ValueError, "Source type is not known. Parser has to be adapted."

            # Record data is separated from descriptions by white-space.
            if header_elements is not None:
                for identifier_type, field in map(None, self.identifier_types, header_elements):
                    if field is not None:
                        identifiers[identifier_type] = field
                    else:
                        identifiers[identifier_type] = ''

                records.append(identifiers)

        # Convert the records to lists.

        converted_records = []
        for identifiers in records:
            converted_record = []

            # Check non-output identifier types to see if they match the label
            # used in the file.
            # NOTE: This is a primitive format test.

            for identifier_type in self.identifier_types:
                if identifier_type not in self.output_identifier_types:
                    if identifiers.get(identifier_type) != identifier_type:
                        raise ValueError, "Identifier type %r was given as %r at line %d." % (identifier_type, identifiers.get(identifier_type), self.lineno + 1)

            # Produce a record from the output identifiers.

            for identifier_type in self.output_identifier_types:
                converted_record.append(identifiers.get(identifier_type, self.NULL))
            converted_records.append(converted_record)
            #print >>sys.stderr, "Parsing identifier %s" % (converted_record)

        return converted_records

    def parse(self):
        records = []
        sequence = []
        self.lineno = 1
        line = self.f.readline()
        while line:
            if line.startswith(">"):
                if records:
                    for record in records:
                        record.append("".join(sequence))
                        self.write_record(record)
                records = self.parse_header(line)
                sequence = []
            else:
                if records != []:
                    sequence.append(line.rstrip("\n"))
            self.lineno += 1
            line = self.f.readline()
        else:
            if records:
                for record in records:
                    record.append("".join(sequence))
                    self.write_record(record)

    def write_record(self, record):
        self.f_main.write("\t".join(record) + "\n")

if __name__ == "__main__":
    from irdata.cmd import get_progname
    from os.path import extsep, join, split, splitext
    import sys, gzip

    progname = get_progname()

    try:
        i = 1
        source = sys.argv[i]
        data_directory = sys.argv[i+1]
        identifier_types = sys.argv[i+2].split(",")
        output_identifier_types = sys.argv[i+3].split(",")
        filenames = sys.argv[i+4:]
    except IndexError:
        print >>sys.stderr, "Usage: %s <source> <output data directory> <identifier types> <output identifier types> <data file>..." % progname
        sys.exit(1)

    filename = None # used for exceptions

    try:
        for filename in filenames:
            leafname = split(filename)[-1]
            basename, ext = splitext(leafname)
            print >>sys.stderr, "%s: Parsing %s" % (progname, leafname)

            if ext.endswith("gz"):
                opener = gzip.open
                basename, ext = splitext(basename) # remove any remaining extension
            else:
                opener = open

            f_out = open(join(data_directory, "%s_proteins.txt" % basename), "w")
            try:
                parser = Parser(source, opener(filename), f_out, identifier_types, output_identifier_types)
                try:
                    parser.parse()
                finally:
                    parser.close()
            finally:
                f_out.close()

    except Exception, exc:
        print >>sys.stderr, "%s: Parsing failed for file %s with exception: %s" % (progname, filename, exc)
        sys.exit(1)

# vim: tabstop=4 expandtab shiftwidth=4
