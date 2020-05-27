#!/usr/bin/python

"""
Open and parse MITAB files and write to tables in preparation for import into 
iRefIndex build.  File format is auto-detected based on number of tab-delimited columns
in first non-commented line.  Recognized file-types are one of
mitab25, mitab26, mitab27, irefindex, or unknown.  Unknown file types will throw an error 
in main.  

irefindex files are recognized but parsing is not yet supported
especially wrt bipartite represented complexes. Output of parser only supports fields 
found in mitab2.5 right now (see global list_fields for example).

Lines containing non protein-protein interactions are not supported and are skipped.

Lines lacking taxon id info for either A or B are skipped

Lines encoding a spoke of a complex record are supported by recognizing repeated 
interaction record ids in the file (the complex representation column present in 2.6 onwards
is not relied on for this).

Reactome distributes files with - instead of distinct interaction record ids : this is supported
by recognizing blank or - entries and replacing them with the line number in the file -
it is assumed that such files will not contain multi-line records (complexes). 

MPIDB distributes files (from their ftp site) with 17 columns that are still supported 
(by recognizing the 17 columns width)
Fixes related to MPIDB are left in mostly to provide as example code such as
-fixing of controlled-vocabulary terms
-support for lines that contain multiple experimental evidences for a single interaction
(see corresponding_fields and associated notes) although this is no longer used since
standard MPIDB mitab files are retrieved using PSICQUIC services. 

detailed notes on the output of this parser are described in the document
irefindex code review
--------

Copyright (C) 2012-2014 Ian Donaldson <ian@donaldsonresearch.com>
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

#for debugging - copy this line to the desired breakpoint
#import pdb; pdb.set_trace()

###################
# Global variables
###################

from os.path import extsep, join, split, splitext
import os
import re
import gzip

#list fields specific to each file format type
mitab25 = (
    #==> applies to MITAB 2.5, 2.6 and 2.7
    "uidA",                     #1.  ID(s) interactor A	
    "uidB",                     #2.  ID(s) interactor B	
    "altA",                     #3.  Alt. ID(s) interactor A	
    "altB",                     #4.  Alt. ID(s) interactor B	
    "aliasA",                   #5.  Alias(es) interactor A	
    "aliasB",                   #6.  Alias(es) interactor B	
    "method",                   #7.  Interaction detection method(s)	
    "authors",                  #8.  Publication 1st author(s)	
    "pmids",                    #9.  Publication Identifier(s)	
    "taxA",                     #10. Taxid interactor A	
    "taxB",                     #11. Taxid interactor B	
    "interactionType",          #12. Interaction type(s)	
    "sourcedb",                 #13. Source database(s)	
    "interactionIdentifiers",   #14. Interaction identifier(s)	
    "confidence"                #15. Confidence value(s)	
    )


mitab26 = (
    #==> only in MITAB 2.6 and 2.7
    "expansion",                #16. Expansion method(s)	
    "bioRoleA",                 #17. Biological role(s) interactor A	
    "bioRoleB",                 #18. Biological role(s) interactor B	
    "expRoleA",                 #19. Experimental role(s) interactor A	
    "expRoleB",                 #20. Experimental role(s) interactor B	
    "typeA",                    #21. Type(s) interactor A	
    "typeB",                    #22. Type(s) interactor B	
    "xrefA",                    #23. Xref(s) interactor A	
    "xrefB",                    #24. Xref(s) interactor B	
    "xrefInt",                  #25. Interaction Xref(s)	
    "annotA",                   #26. Annotation(s) interactor A	
    "annotB",                   #27. Annotation(s) interactor B	
    "annotInt",                 #28. Interaction annotation(s)	
    "hostOrg",                  #29. Host organism(s)	
    "intParam",                 #30. Interaction parameter(s)	
    "created",                  #31. Creation date	
    "updated",                  #32. Update date	
    "checksumA",                #33. Checksum(s) interactor A	
    "checksumB",                #34. Checksum(s) interactor B	
    "checksumInt",              #35. Interaction Checksum(s)	
    "negative"                  #36. Negative
    )	


mitab27 = (
    #==>only in MITAB 2.7
    "featA",                    #37. Feature(s) interactor A	
    "featB",                    #38. Feature(s) interactor B	
    "stoichA",                  #39. Stoichiometry(s) interactor A	
    "stoichB",                  #40. Stoichiometry(s) interactor B	
    "idMethA",                  #41. Identification method participant A	
    "idMethB"                   #42. Identification method participant B
    )


ireff = (
    #==only in iRefIndex MITAB as of release 13
    "origA",                    #37. Original Reference A
    "origB",                    #38. Original Refere3nce B
    "finalA",                   #39. Final Reference A
    "finalB",                   #40. Final Reference B
    "mapA",                     #41. Mapping Score A
    "mapB",                     #42. Mapping Score B
    "irogidA",                  #43. irogida
    "irogidB",                  #44. irogidb
    "irigid",                   #45. irigid
    "crogidA",                  #46. crogida
    "crogidB",                  #47. crogidb
    "crigidB",                  #48. crigid
    "icrogidA",                 #49. icrogida
    "icrogidB",                 #50. icrogidb
    "icrigid",                  #51. icrigid
    "imexid",                   #52. imexid
    "edgetype"                  #53. edgetype
    )

#define fields for each file format type
mitab25_fields = mitab25
mitab26_fields = mitab25 + mitab26
mitab27_fields = mitab25 + mitab26 + mitab27
irefindex_fields =    mitab25 + mitab26 + ireff


#custom mpidb file format - from mpidb ftp site
custom_mpidb_fields = mitab25_fields + (
    "evidence",
    "interaction"
    )

#the defaults below can be changed in main
file_format = "mitab25"
all_fields = mitab25_fields

#see get_experiment_data for an explanation
corresponding_fields = ()

#fields that have pipe-separated lists in mitab
#AND that will be processed
list_fields = ("altA", "altB", "aliasA", "aliasB", "method", "authors", "pmids", 
    "interactionType", "sourcedb", "interactionIdentifiers", "confidence")

mpidb_term_regexp       = re.compile(r'(?P<prefix>.*?)"(?P<term>.*?)"\((?P<description>.*?)\)')
standard_term_regexp    = re.compile(r'(?P<term>.*?)\((?P<description>.*?)\)')
taxid_regexp            = re.compile(r'taxid:(?P<taxid>[^(]+)(\((?P<description>.*?)\))?')

mpidb_taxid_regexp      = re.compile(r'taxid:(?P<taxid>[^(]+)')
uniprotkb_protein_regexp  = re.compile(r'uniprotkb:\"(.+|.+)\"')
term_regexps = [standard_term_regexp]

##################
# Classes
##################

class Parser:

    "A parser which produces output in different forms."

    def __init__(self, writer):
        self.writer = writer
        self.file_line_no = 0
        #keep track of multi-line interction records
        self.int_line_no = 0
        self.last_int_line_no = {}
        #self.last_interaction = None
        

    def close(self):
        if self.writer is not None:
            self.writer.close()
            self.writer = None

    def parse(self, filename):

        """
        Parse the file with the given 'filename', writing to the output stream.
        """

        leafname = split(filename)[-1]
        basename, ext = splitext(leafname)

        if ext.endswith("gz"):
            opener = gzip.open
        else:
            opener = open

        f = opener(filename)

        try:
            self.writer.start(filename)
            
            line = f.readline()
            if file_format == "custom_mpidb":
                line = f.readline() #the first line is an uncommented header

            while line:
                #skip over comment lines
                if line.startswith("#"):
                    line = f.readline()
                    self.file_line_no += 1
                    continue
                #parse all other lines
                 
                self.parse_line(line)
                line = f.readline()
                self.file_line_no += 1

        finally:
            f.close()

    def parse_line(self, line):

        "Parse the given 'line', appending output to the writer."

        data = dict(zip(all_fields, line.strip().split("\t")))
        
        print >>sys.stderr, "line: %s" % line
        # Convert all pipe-delimited list values into lists.
        for key in list_fields:
            data[key] = pipe_2_list(data[key], return_na = True)
        
        print >>sys.stderr, "data: %s" % data
        print >>sys.stderr, "typeA: %s" % data["typeA"] 
        print >>sys.stderr, "typeB: %s" % data["typeB"] 
        # Check that line is suitable for parsing otherwise skip to next line
        # Omit non protein-protein interactions and other records where A or B or taxid are ill-defined
        if data["uidA"] == "-" or data["uidB"] == "-":
            return
        if data["taxA"] == "taxid:-" or  data["taxB"] == "taxid:-":
            return
        if data["taxA"] == "-" or  data["taxB"] == "-":
            return
        if file_format == "mitab26" or file_format == "mitab27":
            if data["typeA"] != 'psi-mi:"MI:0326"(protein)' or  data["typeB"] != 'psi-mi:"MI:0326"(protein)':
                return
        
        # Fix aliases.
        for key in ("aliasA", "aliasB"):
            data[key] = map(fix_alias, data[key])

        # Fix controlled vocabulary fields.
        for key in ("method", "interactionType", "sourcedb"):
            print >>sys.stderr, "key: %s" % key 
            print >>sys.stderr, "data: %s" % data[key] 
            data[key] = map(fix_vocabulary_term, data[key])
            print >>sys.stderr, "data full: %s" % data 
        '''
        Detect multi-line interactions representing protein complexes.
        Distinct interaction records are distinguished based on there being
        a unique interaction record identifier present (see get_interaction_id).
        Some databases lack these in which case get_interaction_id will provide 
        the line number in the mitab file being parsed as a surrogate record indetifier.
        In this later case, it is assumed that each line in the mitab file represents a separate record 
        - multi-line interactions will not be possible for these databases.

        The number of lines in a multi-line intxn record is kept track of using a hash table where
        the unique interaction record id is used as a key (last_int_line_no[interaction_id])      
        
        This will allow multiple lines for the same interaction (complex) record
        to occur non-contiguously in the MITAB file.  The line number for each line of a multi-line 
        interaction record is added to the current data line being processed as data["int_line_number"]
        "last_int_line_no" counts the number of times that an interaction identifier
        has been seen so far during the parse of this MITAB file.
        '''
        data["file_line_no"] = self.file_line_no
        interaction_id = get_interaction_id(data)

        if interaction_id not in self.last_int_line_no:
            self.last_int_line_no[interaction_id] = 0
        else:
            self.last_int_line_no[interaction_id] += 1

        data["int_line_no"] = self.last_int_line_no[interaction_id]


        #finally write the data out to a series of text, tab-delimited tables
        self.writer.append(data)
        
class Writer:

    "Support for writing to files."

    def __init__(self, source, directory):
        self.input_source = source
        self.directory = directory
        self.filename = None
        self.init()

    def start(self, filename):
        self.filename = filename
        self.output_line = 1

    def write_line(self, out, values):
        print >>out, "\t".join(map(str, values))

    def get_experiment_data(self, data):

        """
        Observe correspondences between multivalued fields in 'data'. 
        At present this is only implelmented for the MPIDB source.

        For example...the following three tab-delimited fields occurs in one mitab line
        about an interaction between A and B

        A  B .....method1|method2   author1|author2    pmid1|pmid2.....

        in order to capture two experimental evidences for one interaction 

        the line could be rewritten as two lines where each of the three 
        "corresponding fields" could be written with just one value so the 
        new lines would be

        A  B.....method1    author1    pmid1....
        A  B.....method2    author2    pmid2....
        
        
        This is only applicable at present for custom MPIDB files provided on their ftp site
        and probably does notapply to MPIDB psicquic files
        """
        
        
        # Where no correspondences are being recorded, return the data as the
        # only experiment entry, and with only a single additional entry
        # indicating a unique output line number.

        if not corresponding_fields:
            data["line"] = self.output_line
            self.output_line += 1
            return [data]

        # Obtain the values for each of the fields.

        fields = []
        length = None
        for key in corresponding_fields:
            values = data[key]

            # Ensure a consistent length for all fields.

            if length is None:
                length = len(values)
            elif length != len(values):
                raise ValueError, "Field %s has %d values but preceding fields have %d values." % (key, len(values), length)

            fields.append(values)

        # Get values for the fields for each position in the correspondence.

        experiment_data = []

        for values in zip(*fields):

            # Each value will be on its own in the list of values for the field.

            new_data = {"line" : self.output_line}
            for key, value in zip(corresponding_fields, values):
                new_data[key] = [value]

            # Write the unpacked correspondences.

            experiment_data.append(new_data)
            self.output_line += 1

        return experiment_data

class MITABWriter(Writer):

    '''
    A standard MITAB format file writer.
    Presently not used
    '''

    def init(self):
        self.out = None

    def close(self):
        if self.out is not None:
            self.out.close()
            self.out = None

    def get_filename(self):
        #imd - inspect this later - hard-coded 
        return join(self.directory, "mpidb_mitab.txt")

    def start(self, filename):
        Writer.start(self, filename)

        if self.out is not None:
            return

        if not os.path.exists(self.directory):
            os.mkdir(self.directory)

        self.out = open(self.get_filename(), "w")
        print >>self.out, "#" + "\t".join(mitab25_fields)

    def append(self, data):

        "Write tidied MITAB from the given 'data'."

        for exp_data in self.get_experiment_data(data):
            new_data = {}
            new_data.update(data)
            new_data.update(exp_data)

            # Convert values back into strings.

            for key in list_fields:
                new_data[key] = list_2_pipe(new_data[key])

            self.write_line(self.out, [new_data[key] for key in mitab25_fields])

class iRefIndexWriter(Writer):

    "A writer for iRefIndex-compatible data."

    filenames = (
        "uid", "alias", # collecting more than one column each
        "alternatives", "method", "authors", "pmids", "interactionType", "sourcedb", "interactionIdentifiers"
        )

    def init(self):
        self.files = {}

    def close(self):
        for f in self.files.values():
            f.close()
        self.files = {}

    def get_filename(self, key):
        return join(self.directory, "mitab_%s%stxt" % (key, extsep))

    def start(self, filename):
        Writer.start(self, filename)

        # Use the filename for specific MPIDB sources.

        if file_format == "custom_mpidb":
            self.source = split(filename)[-1]
        else:
            self.source = self.input_source

        if self.files:
            return

        if not os.path.exists(self.directory):
            os.mkdir(self.directory)

        for key in self.filenames:
            self.files[key] = open(self.get_filename(key), "w")

    def append(self, data):


        "Write iRefIndex-compatible output from the 'data'."

        # Interactor-specific records.

        positionA = data["int_line_no"]
        positionB = data["int_line_no"] + 1

        # Only write the principal interactor of multi-line interactions once.

        if positionA == 0:
            self.write_line(self.files["uid"], (self.source, self.filename, get_interaction_id(data), positionA) + get_one_uid(data["uidA"]) + (split_taxid(data["taxA"])[0],))
        self.write_line(self.files["uid"], (self.source, self.filename, get_interaction_id(data), positionB) + get_one_uid(data["uidB"]) + (split_taxid(data["taxB"])[0],))

        for filename, fields in (
            ("alternatives", ("altA", "altB")),
            ("alias", ("aliasA", "aliasB"))
            ):
            for position, key in enumerate(fields):
                for entry, s in enumerate(data[key]):
                    if not s:
                        continue
                    prefix, value = colon_2_tuple(s)

                    # Only write the details of the principal interactor once.

                    if position != 0 or positionA == 0:
                        self.write_line(self.files[filename], (self.source, self.filename, get_interaction_id(data), positionA + position, prefix, value, entry))

        # Experiment-specific records.

        if positionA == 0:
            self.append_lists(self.get_experiment_data(data))

    def append_lists(self, list_data):
        for data in list_data:
            for key in ("authors",):
                for entry, s in enumerate(data[key]):
                    self.write_line(self.files[key], (self.source, self.filename, data["line"], get_interaction_id(data), s, entry))

            for key in ("method", "interactionType", "sourcedb"):
                for entry, s in enumerate(data[key]):
                    term, description = split_vocabulary_term(s)
                    self.write_line(self.files[key], (self.source, self.filename, data["line"], get_interaction_id(data), term, description, entry))

            for key in ("pmids",):
                for entry, s in enumerate(data[key]):
                    prefix, value = colon_2_tuple(s)
                    if prefix == "pubmed":
                        if value.isdigit():
                            
                            self.write_line(self.files[key], (self.source, self.filename, data["line"], get_interaction_id(data), value, entry))
                        if "doi" in value:
                            ##fixing issues with VIRUSHOST example pubmed:https(//doi.org/...)
                            value = value.rstrip(")")
                            value = value.lstrip("https(//doi.org/")
                            self.write_line(self.files[key], (self.source, self.filename, data["line"], get_interaction_id(data), value, entry))
                       
                        else:
                            print >> sys.stderr, " %s is not a pubmedID or DOI: citation suppressed" % (value)
                            #value = '-'
                            #self.write_line(self.files[key], (self.source, self.filename, data["line"], get_interaction_id(data), value, entry))
                    if (prefix.lower()  == "doi"):
                        self.write_line(self.files[key], (self.source, self.filename, data["line"], get_interaction_id(data), value, entry))

            for key in ("interactionIdentifiers",):
                for entry, s in enumerate(data[key]):
                    prefix, value = colon_2_tuple(s)
                    self.write_line(self.files[key], (self.source, self.filename, data["line"], get_interaction_id(data), prefix, value, entry))




##################
# Global functions
##################

def get_interaction_id(data):
    int_id = data.get("interaction") or colon_2_tuple(data["interactionIdentifiers"][0])[-1]
    if int_id == "-" or int_id == "" or int_id == " ":
        int_id = data["file_line_no"]
    return int_id

# Value processing.

def pipe_2_list(s, return_na):
    '''
    converts a pipe-delimited list of strings to a list
    the string "-" is used to denote NA (not available) data in mitab files
    for a string "a|b|-|d" return the list [a,b,"-",d] or only [a,b,d] if return_na is false
    apply fix for uniprot entries with pipe symbol in name
    ''' 
    match =uniprotkb_protein_regexp.search(s)
    if match:
        s = fix_pipe_in_string(s)
    l = s.split("|")
    if return_na:
        return [i for i in l]
    else:
        return [i for i in l if i != "-"]

def list_2_pipe(l):
    '''for the list [a,b,c] return the string "a|b|c" or "-" for empty list'''
    if not l:
        return "-"
    else:
        return "|".join(l)

def fix_pipe_in_string(s):
    '''
    this is a fix for a QUICKGO entry where a pipe symbol is used in the name of the protein
    it originates from UNIPROT Q8VZ68 
    obviously, it could also happen anywhere else
    '''
    for regexp in [uniprotkb_protein_regexp]:
        match = regexp.search(s)
        if match:
            temp_proteins = s.split("|")
            proteins = []

            for item in range(len(temp_proteins)):
                if ':' in temp_proteins[item]:
                    proteins.append(temp_proteins[item])
                else:
                    proteins[item-1] = temp_proteins[item-1] + temp_proteins[item]
         
            return "|".join(proteins) 

def fix_vocabulary_term(s):
    '''
    this is a historical fix for MPIDB and is not required for MPIDB retrieved using psicquic
    '''
    for regexp in term_regexps:
        match = regexp.match(s)

        if match:
            return "%s(%s)" % (match.group("term"), match.group("description"))
    raise ValueError, "Term %r is not well-formed.  Fail in fix_vocabulary_term." % s

def fix_alias(s):
    '''
    this is a historical fix for MPI sources that incorrectly label aliases as from
    uniprotkb when they are in fact entrezgene/locuslink gene names
    '''
    if file_format == "custom_mpidb" and s.startswith("uniprotkb:"):
        prefix, symbol = s.split(":")[:2]
        return "entrezgene/locuslink:" + symbol
    else:
        return s

def colon_2_tuple(s):
    '''
    for the string "a:b" return the tuple (a,b)
    interpret null or - value of s as -:-
    '''
    if s == " " or s == "" or s == "-":
        s = "-:-"
    parts = s.split(":", 1)
    return tuple(parts)

def split_vocabulary_term(s):
    for regexp in term_regexps:
        match = regexp.match(s)
        if match:
            return (match.group("term"), match.group("description"))
    raise ValueError, "Term %r is not well-formed. Fail in split_vocabulary_term." % s

def get_one_uid(s):
    '''
    given a uid for an interactor in the form of 
    a string "a:b" or "a:b|c:d" return just the last tuple (c,d) for InnateDB
    or just the first tuple for all other databases 
    if - is encountered in the original file (i.e. an NA) then return the tuple
    ("-", "-")
    '''
    uids = pipe_2_list(s, return_na = False)
    if uids:
        if source.startswith("INNATE"): # NOTE: Hack for InnateDB MITAB.
            return colon_2_tuple(uids[-1])
        else:
            return colon_2_tuple(uids[0]) 
    else:
        return ('-','-')

def split_taxid(s):
    '''
    given the string 'taxid:1234("some cool organism")' 
    return the tuple taxid=1234 and description="some cool organism" 
    '''
    match = taxid_regexp.match(s)
    if match:
        taxid = match.group("taxid")
        description = match.group("description") or "NA" #allow missing description
        return (taxid, description)
    else:
        raise ValueError, "Taxonomy %r is not well-formed." % s

def detect_file_format(filenames):
    '''
    open and inspect the first of the files to be parsed
    based on the number of columns, a filetype is returned
    one of: mitab25, mitab26, mitab27, irefindex, unknown
    '''
    filename = filenames[0]
    leafname = split(filename)[-1]
    basename, ext = splitext(leafname)

    if ext.endswith("gz"):
        opener = gzip.open
    else:
        opener = open
    f = opener(filename)

    try:
        line = f.readline()
        #skip over comment lines
        if line.startswith("#"):
            line = f.readline()
        data = line.strip().split("\t")
        if len(data) == len(mitab25_fields):
            return("mitab25")
        if len(data) == len(mitab26_fields):
            return("mitab26")
        if len(data) == len(mitab27_fields):
            return("mitab27")
        if len(data) == len(irefindex_fields):
            return("irefindex")
        if len(data) == len(custom_mpidb_fields):
            return ("custom_mpidb")
        else:
            return("unknown")
        
    finally:
        f.close()    

def get_all_fields(file_format):
    '''
    return the list of fields for the parser to retrieve from each line
    of the mitab file based on the file format
    '''
    if file_format == "mitab25" or file_format == "custom_mpidb":
        return mitab25_fields
    elif file_format == "mitab26":
        return mitab26_fields    
    elif file_format == "mitab27":
        return mitab27_fields
    elif file_format == "irefindex":
        return irefindex_fields


###########
# Main
###########    

if __name__ == "__main__":
    from irdata.cmd import get_progname
    import sys

    progname = get_progname()

    try:
        source = sys.argv[1]
        directory = sys.argv[2]
        filenames = sys.argv[3:]
    except IndexError:
        print >>sys.stderr, "Usage: %s <source> <output data directory> <filename>..." % progname
        sys.exit(1)

    try:
        # set up based on type of file to be parsed.
        file_format = detect_file_format(filenames)
        if file_format == "unknown":
            sys.exit(1)
        all_fields = get_all_fields(file_format)
        #change global defaults in the case that MPIDB custom file is being parsed
        if file_format == "custom_mpidb":
            corresponding_fields = (
            "method", "authors", "pmids", "interactionType", "sourcedb",
            "interactionIdentifiers", "confidence")
            all_fields = custom_mpidb_fields
            term_regexps = [mpidb_term_regexp, standard_term_regexp]
        
        #prepare parser and try parsing
        writer = iRefIndexWriter(source, directory)
        parser = Parser(writer)
        try:
            for filename in filenames:
                parser.parse(filename)
        finally:
            parser.close() # closes the writer

    except Exception, exc:
        print >>sys.stderr, "%s: Parsing failed with exception: %s" % (progname, exc)
        sys.exit(1)

# vim: tabstop=4 expandtab shiftwidth=4
