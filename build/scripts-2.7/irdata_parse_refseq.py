#!/usr/bin/python

"""
Parse RefSeq text files.

Various NCBI resources use the "feature table" format documented here:

ftp://ftp.ncbi.nlm.nih.gov/genbank/docs/FTv9_0.html

--------

Copyright (C) 2012 Ian Donaldson <ian.donaldson@biotek.uio.no>
Copyright (C) 2013 Paul Boddie <paul@boddie.org.uk>
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

def parse_line(line):
    if not line:
        raise EOFError
    parts = re.split(r"(\s+)", line)

    if parts:

        # First field is not empty.

        if parts[0]:
            code = parts[0]
            rest = "".join(parts[2:])

        # First field is empty.

        else:
            code = None
            rest = "".join(parts[2:])

        return code, rest.rstrip("\n")
    else:
        return None, None

class Parser:

    "A parser for UniProt text files."

    taxid_regexp = re.compile('/db_xref="taxon:(.*?)"')

    # Optional fields do not typically apply to the published RefSeq archives.
    # They are usually only relevant to old eUtils records.

    optional = "ACCESSION", "VERSION", "TAXID"

    def __init__(self, f, f_main, f_identifiers):
        self.f = f
        self.f_main = f_main
        self.f_identifiers = f_identifiers
        self.line = None
        self.return_last = 0

    def close(self):
        if self.f is not None:
            self.f.close()
            self.f = None
        if self.f_main is not None:
            self.f_main.close()
            self.f_main = None
        if self.f_identifiers is not None:
            self.f_identifiers.close()
            self.f_identifiers = None

    def next_line(self):
        if not self.return_last:
            self.line = parse_line(self.f.readline())
        else:
            self.return_last = 0
        return self.line

    def save_line(self):
        self.return_last = 1

    def parse_locus(self, line, record):
        code, rest = line

        # Filter out DNA. This is useful if having to deal with eUtils.

        if not ' DNA ' in rest:
            record["LOCUS"] = rest.split()[0]

    def parse_accession(self, line, record):
        code, rest = line

        # Handle missing accessions in old eUtils records.

        if rest:
            record[code] = rest.split()[0] # ignore trailing accessions (presumably expired)

    def parse_dbsource(self, line, record):
        code, rest = line
        parts = rest.split()
        if parts[0] == "REFSEQ:":
            record["REFSEQ"] = parts[-1]

    def parse_features(self, line, record):
        code, rest = self.next_line()
        while not code:
            match = self.taxid_regexp.match(rest)
            if match:
                record["TAXID"] = match.group(1)
                return
            code, rest = self.next_line()
        else:
            self.save_line()

    def parse_origin(self, line, record):
        sequence = []
        code, rest = self.next_line()
        while not code:
            rest = "".join(rest.strip().split()[1:])
            sequence.append(rest)
            code, rest = self.next_line()
        else:
            self.save_line()
        record["SEQUENCE"] = "".join(sequence).upper()

    def parse_reference(self, line, record):
        pmids = []
        code, rest = self.next_line()
        while not code:
            fields = rest.split()
            if fields[0] == "PUBMED":
                pmids.append(fields[1])
            code, rest = self.next_line()
        else:
            self.save_line()
        record["PUBMED"] = pmids

    def parse_version(self, line, record):
        code, rest = line
        version_plus_gi = rest.split()

        # Handle missing versions in old eUtils records.

        if len(version_plus_gi) > 1:
            version, gi = version_plus_gi[:2]
            record["VERSION"] = version
            if ':' in gi:
                record["GI"] = gi.split(":")[1] # strip "GI:" from the identifier
        else:
            version = version_plus_gi[0]
        # in the current refseq files no GI in line VERSION
        #record["GI"] = gi.split(":")[1] # strip "GI:" from the identifier
        record["VERSION"]= version
    handlers = {
        "LOCUS"     : parse_locus,
        "ACCESSION" : parse_accession,
        "DBSOURCE"  : parse_dbsource,
        "FEATURES"  : parse_features,
        "ORIGIN"    : parse_origin,
        "REFERENCE" : parse_reference,
        "VERSION"   : parse_version,
        }

    def parse(self):
        record = {}
        try:
            while 1:
                line = self.next_line()
                code, rest = line
                if code == "//":
                    self.write_record(record)
                    record = {}
                handler = self.handlers.get(code)
                if handler is not None:
                    handler(self, line, record)
        except EOFError:
            pass

    def write_record(self, record):
        if record.has_key("LOCUS"):

            # Handle missing fields in old eUtils records.
            for key in self.optional:
                if not record.has_key(key):
                    record[key] = r"\N"

            self.f_main.write("%(ACCESSION)s\t%(VERSION)s\t%(TAXID)s\t%(SEQUENCE)s\n" % record)
            if record.has_key("PUBMED"):
                for pos, pmid in enumerate(record["PUBMED"]):
                    self.f_identifiers.write("%s\tPubMed\t%s\t%s\n" % (record["ACCESSION"], pmid, pos))

if __name__ == "__main__":
    from irdata.cmd import get_progname
    from os.path import extsep, join, split, splitext
    import sys, gzip

    progname = get_progname()

    try:
        i = 1
        data_directory = sys.argv[i]
        filenames = sys.argv[i+1:]
    except IndexError:
        print >>sys.stderr, "Usage: %s <output data directory> <data file>..." % progname
        sys.exit(1)

    filename = None # used for exceptions

    try:
        for filename in filenames:
            leafname = split(filename)[-1]
            basename, ext = splitext(leafname)
            print >>sys.stderr, "%s: Parsing %s" % (progname, leafname)

            if ext.endswith("gz"):
                opener = gzip.open
            else:
                opener = open

            f_main = open(join(data_directory, basename + "_proteins"), "w")
            f_identifiers = open(join(data_directory, basename + "_identifiers"), "w")
            parser = Parser(opener(filename), f_main, f_identifiers)
            try:
                parser.parse()
            finally:
                parser.close()

    except Exception, exc:
        print >>sys.stderr, "%s: Parsing failed for file %s with exception: %s" % (progname, filename, exc)
        sys.exit(1)

# vim: tabstop=4 expandtab shiftwidth=4
