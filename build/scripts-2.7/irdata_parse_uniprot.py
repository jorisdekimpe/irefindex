#!/usr/bin/python

"""
Parse UniProt text files.

See: http://web.expasy.org/docs/userman.html

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

def parse_line(line):
    if not line:
        raise EOFError
    code, space, rest = line[:2], line[2:5], line[5:]
    if code != "//" and space != "   ":
        raise ValueError, "Line was not of the form 'code   data': %r" % line
    return code, rest.rstrip("\n")

class Parser:

    "A parser for UniProt text files."

    null = r"\N"
    date_regexp = re.compile(r"(\d{2})-([A-Z]{3})-(\d{4})")
    pubmed_regexp = re.compile(r"PubMed=(\d+);")
    gene_name_regexp = re.compile(r"Name=(.+?);")

    months = {}
    for i, name in enumerate(("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")):
        months[name] = i + 1

    def __init__(self, f, f_main, f_accessions, f_identifiers, f_gene_names):
        self.f = f
        self.f_main = f_main
        self.f_accessions = f_accessions
        self.f_identifiers = f_identifiers
        self.f_gene_names = f_gene_names
        self.line = None
        self.return_last = 0

    def close(self):
        if self.f_main is not None:
            self.f_main.close()
            self.f_main = None
        if self.f_accessions is not None:
            self.f_accessions.close()
            self.f_accessions = None
        if self.f_identifiers is not None:
            self.f_identifiers.close()
            self.f_identifiers = None
        if self.f_gene_names is not None:
            self.f_gene_names.close()
            self.f_gene_names = None

    def next_line(self):
        if not self.return_last:
            self.line = parse_line(self.f.readline())
        else:
            self.return_last = 0
        return self.line

    def save_line(self):
        self.return_last = 1

    def parse_identifier(self, line):
        code, rest = line
        return rest.split()[0]

    def parse_accessions(self, line):
        code, rest = line
        l = self._get_accessions(rest)
        code, rest = self.next_line()
        while code == "AC":
            l += self._get_accessions(rest)
            code, rest = self.next_line()
        else:
            self.save_line()
        return l

    def _get_accessions(self, s):
        return [i.strip() for i in s.split(";") if i.strip()]

    def parse_dates(self, line):

        "See: http://web.expasy.org/docs/userman.html#DT_line"

        code, rest = line
        creation = self._get_date(rest) # creation date
        sequence = None
        code, rest = self.next_line()
        while code == "DT":
            if not sequence:
                sequence = self._get_date(rest)
            code, rest = self.next_line()
        else:
            self.save_line()
        return sequence or creation or self.null

    def _get_date(self, s):
        match = self.date_regexp.match(s)
        if match:
            day, month, year = match.groups()
            return "%s%02d%s" % (year, self.months[month], day)
        else:
            return None

    def parse_taxonomy(self, line):

        "See: http://web.expasy.org/docs/userman.html#OX_line"

        code, rest = line
        key, value = rest.rstrip(";").split(" ")[0].split("=")
        if key == "NCBI_TaxID":
            #taxid = value.split(" ") #separate possible extra annotation
            return value 
        else:
            return None

    def parse_sequence(self, line):

        "See: http://web.expasy.org/docs/userman.html#SQ_line"

        # Get the molecular weight from the SQ line.

        code, rest = line
        stats = rest.split(";")
        length_part, mw_part, crc_part, _empty = stats
        mw = mw_part.strip().split()[0]

        # Get the sequence from subsequent lines.

        sequence = []
        code, rest = self.next_line()
        while code == "  ":
            sequence.append(rest)
            code, rest = self.next_line()
        else:
            self.save_line()

        # Return the molecular weight and sequence.

        return mw, "".join(sequence).replace(" ", "")

    def parse_pubmed(self, line):

        "See: http://web.expasy.org/docs/userman.html#RX_line"

        code, rest = line
        pmids = []
        while code == "RX":
            pmid = self._get_pubmed(rest)
            if pmid:
                pmids.append(pmid)
            code, rest = self.next_line()
        else:
            self.save_line()
        return pmids

    def _get_pubmed(self, s):
        match = self.pubmed_regexp.match(s)
        if match:
            return match.groups()[0]
        else:
            return None

    def parse_identifiers(self, line):

        "See: http://web.expasy.org/docs/userman.html#DR_line"

        code, rest = line
        identifiers = set()
        while code == "DR":
            parts = rest.rstrip(".").split(";")

            # Only the first identifier is recorded.

            if len(parts) > 1:
                type = parts[0].strip()
                value = parts[1].strip()
                identifiers.add((type, value))
            code, rest = self.next_line()
        else:
            self.save_line()
        return identifiers

    def parse_gene_names(self, line):

        "See: http://web.expasy.org/docs/userman.html#GN_line"

        code, rest = line
        names = set()
        while code == "GN":
            gene_name = self._get_gene_name(rest)
            if gene_name:
                names.add(gene_name)
            code, rest = self.next_line()
        else:
            self.save_line()
        return names

    def _get_gene_name(self, s):
        match = self.gene_name_regexp.search(s)
        if match:
            return match.groups()[0]
        else:
            return None

    handlers = {
        "ID" : parse_identifier,
        "AC" : parse_accessions,
        "DT" : parse_dates,
        "OX" : parse_taxonomy,
        "SQ" : parse_sequence,
        "RX" : parse_pubmed,
        "DR" : parse_identifiers,
        "GN" : parse_gene_names,
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
                    record[code] = handler(self, line)
        except EOFError:
            pass

    def write_record(self, record):
        record["AC1"] = record["AC"][0]
        record["MW"] = record["SQ"][0]
        record["SQ"] = record["SQ"][1]
        self.f_main.write("%(ID)s\t%(AC1)s\t%(DT)s\t%(OX)s\t%(MW)s\t%(SQ)s\n" % record)

        # Write accessions for each identifier.

        for accession in record["AC"]:
            self.f_accessions.write("%s\t%s\n" % (record["ID"], accession))

        # Write PubMed references.

        if record.has_key("RX"):
            for pos, pmid in enumerate(record["RX"]):
                self.f_identifiers.write("%s\tPubMed\t%s\t%s\n" % (record["ID"], pmid, pos))

        # Write cross-references.

        if record.has_key("DR"):
            for type, identifier in record["DR"]:
                self.f_identifiers.write("%s\t%s\t%s\t0\n" % (record["ID"], type, identifier))

        # Write gene names.

        if record.has_key("GN"):
            for pos, name in enumerate(record["GN"]):
                self.f_gene_names.write("%s\t%s\t%s\n" % (record["ID"], name, pos))

if __name__ == "__main__":
    from irdata.cmd import get_progname
    from os.path import extsep, join, split, splitext
    import sys, gzip, re

    progname = get_progname()

    try:
        i = 1
        data_directory = sys.argv[i]
        filename = sys.argv[i+1]
        if len(sys.argv) > i+2:
            format = sys.argv[i+2]
        else:
            format = None
    except IndexError:
        print >>sys.stderr, "Usage: %s <output data directory> <data file> [ <format> ]" % progname
        sys.exit(1)

    try:
        leafname = split(filename)[-1]
        basename, ext = splitext(leafname)

        # Use a default format like "uniprot_sprot_%s.txt".

        format = format or ("%s" % basename) + "_%s" + ("%stxt" % extsep)
        mainfile = join(data_directory, format % "proteins")
        accessionsfile = join(data_directory, format % "accessions")
        identifiersfile = join(data_directory, format % "identifiers")
        genenamesfile = join(data_directory, format % "gene_names")

        if filename == "-":
            print >>sys.stderr, "%s: Parsing standard input" % progname
            f = sys.stdin
        else:
            print >>sys.stderr, "%s: Parsing %s" % (progname, leafname)
            if ext.endswith("gz"):
                opener = gzip.open
            else:
                opener = open
            f = opener(filename)

        parser = Parser(f, open(mainfile, "w"), open(accessionsfile, "w"), open(identifiersfile, "w"), open(genenamesfile, "w"))
        try:
            parser.parse()
        finally:
            parser.close()
            if filename != "-":
                f.close()

    except Exception, exc:
        print >>sys.stderr, "%s: Parsing failed for file %s with exception: %s" % (progname, filename, exc)
        sys.exit(1)

# vim: tabstop=4 expandtab shiftwidth=4
