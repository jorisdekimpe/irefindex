#!/usr/bin/python

"""
A tool which reads PSI MI XML files and produces tabular data.

PSI MI XML files can provide separate experiment, interaction and interactor
lists:

  experimentList
    experimentDescription
  interactionList
    interaction
      experimentList
        experimentRef -> experimentDescription/@id
      participantList
        participant
          interactorRef -> interactor/@id
  interactorList
    interactor

Or such files can provide interaction lists containing experiment and interactor
details:

  interactionList
    interaction
      experimentList
        experimentDescription
      participantList
        participant
          interactor

When processing both kinds of files, properties of each data type can be
captured as they are read. The current interaction identifier must be retained
in order to document the relationships between interactions and the other data
types.

For the first kind of file, interaction relationships to experiments and
interactors are explicitly given in "*Ref" elements. For the second kind of
file, such relationships are implicit when an experiment or interactor is
included within an interaction.

Participant properties are defined in terms of an interactor as part of an
interaction. Participants are always implicitly referenced.

--------

Copyright (C) 2011, 2012, 2013 Ian Donaldson <ian.donaldson@biotek.uio.no>
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

from irdata.data import *
from irdata.signatures import make_signature, normalise_sequence
from irdata.xmldata import *
import os

class PSIParser(EmptyElementParser):

    """
    A class which records the properties and relationships in PSI MI XML files.
    """

    # Attributes of supported elements.
    # Values of these attributes are read for supported elements and written to
    # output files. Some attributes are actually constructed internally but are
    # listed here as part of the output specification (property and element).

    attribute_names = {
        # references    : property, reftype, id, dblabel, dbcode, reftypelabel, reftypecode
        "primaryRef"    : ("property", "element", "id", "db", "dbAc", "refType", "refTypeAc"), # also secondary and version
        "secondaryRef"  : ("property", "element", "id", "db", "dbAc", "refType", "refTypeAc"),
        # names         : property, nametype, label, code, value
        "shortLabel"    : ("property", "element", None, None, "content"),
        "fullName"      : ("property", "element", None, None, "content"),
        "alias"         : ("property", "element", "type", "typeAc", "content"),
        # organisms     : taxid
        "hostOrganism"  : ("ncbiTaxId",)
        }

    # Elements defining scopes/entities.

    scopes = {
        "entry"                 : "entry",
        "interaction"           : "interaction",
        "interactor"            : "interactor",
        "participant"           : "participant",
        "experimentDescription" : "experimentDescription",

        # PSI MI XML version 1.0 element mappings.

        "proteinInteractor"     : "interactor",
        "proteinParticipant"    : "participant",
        }

    def __init__(self, writer):
        EmptyElementParser.__init__(self)
        self.writer = writer

    def reset(self):

        # For transient identifiers.

        self.identifiers = {
            "entry"                 : 0,
            "interaction"           : 0,
            "interactor"            : 0,
            "participant"           : 0,
            "experimentDescription" : 0
            }

    def get_scopes(self, n=None):

        "Return the scopes applying to the current path."

        scopes = []

        # Go through the path from the deepest element name to the root, looking
        # for scope names.
        print >>sys.stderr, "Scoping "

        for part in self.current_path[-1::-1]:
            if part in self.scopes.values():
                scopes.append(part)

                # Stop collecting scopes if n have been found.

                if n is not None and len(scopes) == n:
                    break

        # Pad the list with None values if n is greater than the number of
        # scopes found.

        if n is not None and len(scopes) < n:
            scopes += [None] * (n - len(scopes))

        return scopes

    def is_implicit(self, name, context):

        """
        Return whether the element with the given 'name' defines an implicit
        (not externally referenced) element, given the 'context' element name.
        """

        return name == "participant" or \
            name == "interactor" and context == "participant" or \
            name == "experimentDescription" and context == "interaction"

    def characters(self, content):

        "Handle character 'content' by stripping white-space from the ends."

        EmptyElementParser.characters(self, content.strip())

    def startElement(self, name, attrs):

        """
        Start an element, converting the element 'name' to a recognised scope if
        necessary, and adding an identifier to the 'attrs' if one is missing.

        This effectively normalises the availability of identifiers on various
        elements and makes sure that identifiers seen by handleElement are ones
        that are usable for subsequent processing.
        """

        if self.scopes.has_key(name):
            name = self.scopes[name]

            if self.identifiers.has_key(name):
                context = self.get_scopes(1)[0]

                # Handle PSI MI XML 1.0 identifiers which are absent.
                # Also assign identifiers to entries.

                # Use transient participant identifiers since these might be
                # reused within interactions (seen in InnateDB).

                # Also use transient interactor identifiers where their
                # relationship to participants is implicit, since these might be
                # reused within interactions (seen in InnateDB).

                if not attrs.has_key("id") or self.is_implicit(name, context):
                    attrs = dict(attrs)
                    attrs["id"] = str(self.identifiers[name])
                    self.identifiers[name] += 1

        # Start the element using a scope as a logical name if appropriate.
        # This normalises the element names so that they can be treated like
        # PSI MI XML 2.5 elements.

        EmptyElementParser.startElement(self, name, attrs)

    def endElement(self, name):

        "End the element using a scope as a logical name in place of 'name'."

        EmptyElementParser.endElement(self, self.scopes.get(name, name))

    def handleElement(self, content):

        "Handle a completed element with the given 'content'."

        if "entry" not in self.current_path:
            return

        # Get the element names in order of decreasing locality, padding with
        # None. Here, the current path of element names descending into the
        # document is reversed and padded with None, so that the outermost
        # elements can be missing in certain cases.

        # Some examples:

        # entry/interactionList/interaction/participant/interactor
        # -> interactor (element), participant (parent), interaction (property),
        #    interactionList (section)

        # entry/experimentList/experimentDescription
        # -> experimentDescription (element), experimentList (parent),
        #    entry (property), None (section)
        
        element, parent, property, section = map(lambda x, y: x or y, self.current_path[-1:-5:-1], [None] * 4)

        # Get the element's attributes.

        attrs = dict(self.current_attrs[-1])

        # Remember the entry element's identifier value for subsequent use.

        entry = self.path_to_attrs["entry"]["id"]

        # Get mappings from experiments to interactions.
        # The "ref" attribute is from PSI MI XML 1.0.

        if element == "experimentRef":
            if parent == "experimentList":
                self.writer.append((element, entry, content or attrs["ref"], self.path_to_attrs["interaction"]["id"]))

        # And mappings from interactors to participants to interactions.
        # The "ref" attribute is from PSI MI XML 1.0.

        elif element == "interactorRef":
            if parent == "participant":
                self.writer.append((element, entry, content or attrs["ref"], "explicit", self.path_to_attrs["participant"]["id"], self.path_to_attrs["interaction"]["id"]))

        # Implicit interactor-to-participant mappings (applying only within participant elements).

        elif element == "interactor":
            if parent == "participant":
                self.writer.append((element, entry, attrs["id"], "implicit", self.path_to_attrs["participant"]["id"], self.path_to_attrs["interaction"]["id"]))

        # Implicit mappings applying only within an interaction scope.

        elif element == "experimentDescription":
            if self.path_to_attrs.has_key("interaction"):
                self.writer.append((element, entry, attrs["id"], self.path_to_attrs["interaction"]["id"]))

        # Interactor organisms.

        elif element == "organism":
            if parent == "interactor":
                implicit = self.is_implicit(parent, property) and "implicit" or "explicit"
                self.writer.append((element, entry, parent, self.path_to_attrs["interactor"]["id"], implicit, attrs["ncbiTaxId"]))

        # Sequence data.
        # Use the "legacy" mode to upper-case and strip white-space from interactors.
        # Note that this should never be used for interaction signatures/digests.

        elif element == "sequence":
            if parent == "interactor":
                implicit = self.is_implicit(parent, property) and "implicit" or "explicit"
                self.writer.append((element, entry, parent, self.path_to_attrs["interactor"]["id"], implicit, normalise_sequence(content), make_signature(content, legacy=1)))

        # Get other data. This is of the form...
        # section/property/parent/element
        # For example:
        # interactorList/interactor/xref/primaryRef

        else:
            # Only consider supported elements.

            names = self.attribute_names.get(element)
            if not names:
                return

            # Exclude certain element occurrences (as also done above).
            # Such occurrences do not define entities and are therefore not of
            # interest.

            if property == "interactor" and section not in ("participant", "interactorList") or \
                property == "participant" and section != "participantList":
                return

            # Insist on a scope.

            scope, context = self.get_scopes(2)
            if not scope or scope == "entry":
                return

            # Determine whether the information is provided as part of separate
            # (explicit) or embedded (implicit) definitions.

            implicit = self.is_implicit(scope, context) and "implicit" or "explicit"

            # Gather together attributes.

            if content:
                attrs["content"] = content

            # Get the property and element.

            attrs["property"] = property
            attrs["element"] = element

            # Copy the required attributes.

            values = []
            for key in names:
                values.append(attrs.get(key))

            # Only write data for supported elements providing data.

            if not values:
                return

            # The parent indicates the data type and is only used to select the output file.

            self.writer.append((parent, entry, scope, self.path_to_attrs[scope]["id"], implicit) + tuple(values))

    def parse(self, filename):
        self.reset()
        self.writer.start(filename)
        EmptyElementParser.parse(self, filename)

class PSIWriter(Writer):

    "A simple writer of tabular data."

    filenames = (
        "experiment", "interactor",     # mappings
        "names", "xref", "organisms",   # properties
        "sequences"                     # properties
        )

    data_type_files = {
        "experimentRef"         : "experiment",
        "experimentDescription" : "experiment",
        "interactorRef"         : "interactor",
        "interactor"            : "interactor",
        "hostOrganismList"      : "organisms",
        "organism"              : "organisms",
        "sequence"              : "sequences",
        "names"                 : "names",
        "xref"                  : "xref",
        }

    def __init__(self, directory, source):
        Writer.__init__(self, directory)
        self.source = source

    def append(self, data):

        """
        Write out the given 'data', using the first element of 'data' to
        determine the data type.
        """

        element = data[0]
        file = self.data_type_files[element]

        # Each record is prefixed with the source and filename.

        data = (self.source, self.filename) + data[1:]
        data = map(tab_to_space, data)
        data = map(bulkstr, data)
        print >>self.files[file], "\t".join(data)

if __name__ == "__main__":
    from irdata.cmd import get_progname
    import sys

    progname = get_progname()

    try:
        i = 1
        data_directory = sys.argv[i]
        source = sys.argv[i+1]
        filenames = sys.argv[i+2:]
    except IndexError:
        print >>sys.stderr, "Usage: %s <data directory> <data source name> <data file>..." % progname
        sys.exit(1)

    writer = PSIWriter(data_directory, source)
    writer.reset()

    parser = PSIParser(writer)
    try:
        try:
            for filename in filenames:
                parser.parse(filename)
        except Exception, exc:
            print >>sys.stderr, "%s: Parsing failed with exception: %s" % (progname, exc)
            sys.exit(1)
    finally:
        writer.close()

# vim: tabstop=4 expandtab shiftwidth=4
