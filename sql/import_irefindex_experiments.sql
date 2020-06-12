-- Collect experiment-related information.

-- Copyright (C) 2012, 2013 Ian Donaldson <ian.donaldson@biotek.uio.no>
-- Original author: Paul Boddie <paul.boddie@biotek.uio.no>
--
-- This program is free software; you can redistribute it and/or modify it under
-- the terms of the GNU General Public License as published by the Free Software
-- Foundation; either version 3 of the License, or (at your option) any later
-- version.
--
-- This program is distributed in the hope that it will be useful, but WITHOUT ANY
-- WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
-- PARTICULAR PURPOSE.  See the GNU General Public License for more details.
--
-- You should have received a copy of the GNU General Public License along
-- with this program.  If not, see <http://www.gnu.org/licenses/>.

begin;

-- Get all experiment-related records of interest.

insert into xml_xref_all_experiments
    select source, filename, entry, parentid as experimentid,
        case when property = 'interactionDetection' then 'interactionDetectionMethod'
             else property
        end as property,
        reftype,
        case when dblabel in ('Pub-Med', 'PUBMED', 'Pubmed', 'PubMed') then 'pubmed'
             when dblabel in ('MI', 'psimi', 'PSI-MI') then 'psi-mi'
	     when dblabel in ('doi', 'DOI') then 'doi'
             else dblabel
        end as dblabel,

        -- Fix certain psi-mi references.

        case when dblabel = 'MI' and not refvalue like 'MI:%' then 'MI:' || refvalue
             when dblabel = 'MI' and refvalue like 'MI:%' and not refvalue ~ 'MI:[0-9]{4}$' then substring(refvalue from 'MI:[0-9]{4}')
             else refvalue
        end as refvalue

    from xml_xref
    where scope = 'experimentDescription'
        and reftype in ('primaryRef', 'secondaryRef')
        and property in ('bibref', 'interactionDetection', 'interactionDetectionMethod', 'participantIdentificationMethod');

analyze xml_xref_all_experiments;

insert into xml_xref_experiment_organisms
    select distinct source, filename, entry, parentid as experimentid, taxid
    from xml_organisms
    where scope = 'experimentDescription';

analyze xml_xref_experiment_organisms;

-- NOTE: A few records maintain primary and secondary references to the same
-- NOTE: article and thus the primary reference type is chosen.

insert into xml_xref_experiment_pubmed
    select source, filename, entry, experimentid, min(reftype) as reftype, refvalue
    from xml_xref_all_experiments
    where property = 'bibref' and dblabel = 'pubmed'
        and refvalue ~ E'^[[:digit:]]\+$'
        and refvalue <> '0'
    group by source, filename, entry, experimentid, refvalue
    union all
	select source, filename, entry, experimentid, min(reftype) as reftype, refvalue
    from xml_xref_all_experiments
    where property = 'bibref' and dblabel = 'doi'
        and refvalue ~ '\d+.\d+/[a-zA-Z0-9\.\:]'
        and refvalue <> '0'
    group by source, filename, entry, experimentid, refvalue;

analyze xml_xref_experiment_pubmed;

insert into xml_xref_experiment_methods
    select distinct source, filename, entry, experimentid, property, refvalue
    from xml_xref_all_experiments
    where property in ('interactionDetectionMethod', 'participantIdentificationMethod') and dblabel = 'psi-mi';

analyze xml_xref_experiment_methods;

-- Author information originates in the short labels, but is not consistently recorded.
-- NOTE: A list of usable sources is included here.

insert into xml_names_experiment_authors
    select source, filename, entry, parentid as experimentid, name
    from xml_names
    where scope = 'experimentDescription'
        and property = 'experimentDescription'
        and nametype = 'shortLabel'
        and source in ('BIOGRID', 'INTACT', 'MINT');

analyze xml_names_experiment_authors;

-- Some method information is also available in short labels and full names.

insert into xml_names_experiment_methods
    select distinct source, filename, entry, parentid as experimentid, property, nametype, name
    from xml_names
    where scope = 'experimentDescription'
        and property in ('interactionDetectionMethod', 'participantIdentificationMethod');

analyze xml_names_experiment_methods;

commit;
