-- Produce MITAB data from the iRefIndex build database.

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

-- Create a table of interactor pairs for each source database interaction.
-- Although one could group records according to the rigid and uidA, uidB pairs,
-- producing a consolidated table of interaction components, the data remains
-- specific to the source interaction information so that related information
-- (such as experimental details) can be associated specifically with a
-- particular source record.
--
-- For interactions with one or two participants, a single line is written:
--
-- One participant, A:         uidA = A, uidB = A
-- Two participants, A and B:  uidA = A, uidB = B
--
-- For interactions with more than two participants (complexes), as many lines
-- are written as participants:
--
-- Many participants, A...N:   uidA = rigid, uidB = A; ...; uidA = rigid, uidB = N

create temporary table tmp_interactions as

    -- One participant.

    select I.source, I.filename, I.entry, I.interactionid, rigid,
        interactorid as interactoridA, interactorid as interactoridB,
        participantid as participantidA, participantid as participantidB,
        rogid as uidA, rogid as uidB,
        cast('Y' as varchar) as edgetype, numParticipants
    from irefindex_interactions as I
    inner join (
        select source, filename, entry, interactionid, count(participantid) as numParticipants
        from irefindex_interactions
        group by source, filename, entry, interactionid
        having count(participantid) = 1
        ) as Y
        on (I.source, I.filename, I.entry, I.interactionid) = (Y.source, Y.filename, Y.entry, Y.interactionid)
    union all

    -- Two participants.

    select I.source, I.filename, I.entry, I.interactionid, rigid,
        detailsA[2] as interactoridA, detailsB[2] as interactoridB,
        detailsA[3] as participantidA, detailsB[3] as participantidB,
        detailsA[1] as uidA, detailsB[1] as uidB,
        cast('X' as varchar) as edgetype, numParticipants
    from (
        select source, filename, entry, interactionid, rigid, count(participantid) as numParticipants,
            min(array[rogid, interactorid, participantid]) as detailsA,
            max(array[rogid, interactorid, participantid]) as detailsB
        from irefindex_interactions
        group by source, filename, entry, interactionid, rigid
        having count(participantid) = 2
        ) as I
    union all

    -- Many participants.

    select I.source, I.filename, I.entry, I.interactionid, rigid,
        cast(null as varchar) as interactoridA, interactorid as interactoridB,
        cast(null as varchar) as participantidA, participantid as participantidB,
        rigid as uidA, rogid as uidB,
        cast('C' as varchar) as edgetype, numParticipants
    from irefindex_interactions as I
    inner join (
        select source, filename, entry, interactionid, count(participantid) as numParticipants
        from irefindex_interactions
        group by source, filename, entry, interactionid
        having count(participantid) > 2
        ) as C
        on (I.source, I.filename, I.entry, I.interactionid) = (C.source, C.filename, C.entry, C.interactionid);

analyze tmp_interactions;

-- Define all source database identifiers for each ROG identifier.
-- Any | characters will be replaced since the identifiers will be used in
-- pipe-separated lists.

create temporary table tmp_identifiers as
    select rogid, array_accum(distinct dblabel || ':' || replace(refvalue, '|', '_')) as names
    from irefindex_rogid_identifiers
    group by rogid;

analyze tmp_identifiers;

-- Define aliases for each ROG identifier.
-- Any | characters will be replaced since the identifiers will be used in
-- pipe-separated lists.
-- Gene symbols are prefixed with "hgnc" in order to avoid confusion with gene
-- identifiers which are prefixed with "entrezgene/locuslink".

create temporary table tmp_aliases as
    select rogid, array_accum(distinct dblabel || ':' || replace(refvalue, '|', '_')) as aliases
    from (
        select rogid, dblabel, uniprotid as refvalue
        from irefindex_rogid_identifiers
        inner join uniprot_accessions
            on dblabel = 'uniprotkb'
            and refvalue = accession
        union all
        select rogid, 'hgnc' as dblabel, symbol as refvalue
        from irefindex_rogid_identifiers
        inner join gene_info
            on dblabel = 'entrezgene/locuslink'
            and cast(refvalue as integer) = geneid
        ) as X
    group by rogid;

analyze tmp_aliases;

-- Accumulate role collections.
-- Each role is encoded as "MI:NNNN(...)".

create temporary table tmp_participants as
    select source, filename, entry, participantid, property, array_accum(distinct refvalue || '(' || coalesce(name, '-') || ')') as refvalues
    from xml_xref_participants
    left outer join psicv_terms
        on refvalue = code
        and nametype = 'preferred'
    group by source, filename, entry, participantid, property;

alter table tmp_participants add primary key(source, filename, entry, participantid, property);
analyze tmp_participants;

-- Accumulate methods.
-- Each role is encoded as "MI:NNNN(...)".

create temporary table tmp_methods as
    select source, filename, entry, experimentid, property, array_accum(distinct refvalue || '(' || coalesce(name, '-') || ')') as refvalues
    from xml_xref_experiment_methods
    left outer join psicv_terms
        on refvalue = code
        and nametype = 'preferred'
    group by source, filename, entry, experimentid, property;

alter table tmp_methods add primary key(source, filename, entry, experimentid, property);
analyze tmp_methods;

-- Accumulate PubMed identifiers.

create temporary table tmp_pubmed as
    select source, filename, entry, experimentid, array_accum(distinct 'pubmed:' || refvalue) as refvalues
    from xml_xref_experiment_pubmed
    group by source, filename, entry, experimentid;

alter table tmp_pubmed add primary key(source, filename, entry, experimentid);
analyze tmp_pubmed;

-- Consolidate assignment information to get full details of preferred assignments.

create temporary table tmp_assignments as
    select A.*, score
    from irefindex_assignments_preferred as P
    inner join irefindex_assignments as A
        on (P.source, P.filename, P.entry, P.interactorid, P.sequencelink, P.dblabel, P.refvalue, P.finaldblabel, P.finalrefvalue) =
           (A.source, A.filename, A.entry, A.interactorid, A.sequencelink, A.dblabel, A.refvalue, A.finaldblabel, A.finalrefvalue)
    inner join irefindex_assignment_scores as S
        on (P.source, P.filename, P.entry, P.interactorid) = (S.source, S.filename, S.entry, S.interactorid);

alter table tmp_assignments add primary key(source, filename, entry, interactorid);
analyze tmp_assignments;

-- Accumulate confidence score information.

create temporary table tmp_confidence as
    select I.rigid, array_accum(distinct confI.scoretype || ':' || cast(confI.score as varchar)) as confidence
    from tmp_interactions as I
    left outer join irefindex_confidence as confI
        on I.rigid = confI.rigid
    group by I.rigid;

alter table tmp_confidence add primary key(rigid);
analyze tmp_confidence;

-- Collect all interaction-related information.

create temporary table tmp_named_interactions as
    select I.*,

        -- interactionIdentifier (includes rigid, irigid, and edgetype as "X", "Y" or "C",
        -- with "source:-" if no database-specific identifier is provided)

        case when nameI.dblabel is null then coalesce(sourceI.name, lower(I.source)) || ':-'
             else nameI.dblabel || ':' || nameI.refvalue
        end || '|rigid:' || I.rigid || '|edgetype:' || I.edgetype as interactionIdentifier,

        -- sourcedb (as "MI:code(name)" using "MI:0000(name)" for non-CV sources)

        coalesce(sourceI.code, 'MI:0000') || '(' || coalesce(sourceI.name, lower(I.source)) || ')' as sourcedb,

        -- confidence

        array_to_string(confI.confidence, '|') as confidence

    from tmp_interactions as I

    -- Interaction identifier information.

    left outer join xml_xref_interactions as nameI
        on (I.source, I.filename, I.entry, I.interactionid) = (nameI.source, nameI.filename, nameI.entry, nameI.interactionid)

    -- Source database information.

    left outer join psicv_terms as sourceI
        on lower(I.source) = sourceI.name
        and sourceI.nametype = 'preferred'

    -- Confidence score information.

    inner join tmp_confidence as confI
        on I.rigid = confI.rigid;

analyze tmp_named_interactions;

-- Combine interaction and experiment information.

create temporary table tmp_interaction_experiments as
    select I.*,

        -- hostOrganismTaxid (as "taxid:...")

        case when taxidE is null then '-'
             else 'taxid:' || taxidE.taxid || '(' || coalesce(taxnamesE.name, '-') || ')'
        end as hostOrganismTaxid,

        -- method (interaction detection method as "MI:code(name)")

        case when methodE.refvalues is null or array_length(methodE.refvalues, 1) = 0 then '-'
             else array_to_string(methodE.refvalues, '|')
        end as method,

        -- authors (as "name-[year[-number]]")
        -- NOTE: Not converting the form of the authors.

        coalesce(authorE.name, '-') as authors,

        -- pmids (as "pubmed:...")

        case when pubmedE.refvalues is null or array_length(pubmedE.refvalues, 1) = 0 then '-'
             else array_to_string(pubmedE.refvalues, '|')
        end as pmids,

        -- interactionType (interaction type as "MI:code(name)")

        case when typeI.refvalue is null then '-'
             when typenameI.name is null then 'MI:0000(' || typeI.refvalue || ')'
             else typeI.refvalue || '(' || typenameI.name || ')'
        end as interactionType,

        -- imexid

        case when imexI.refvalue is null then '-'
             else 'imex:' || imexI.refvalue
        end as imexid

    from tmp_named_interactions as I
    inner join xml_experiments as E
        on (I.source, I.filename, I.entry, I.interactionid) = (E.source, E.filename, E.entry, E.interactionid)

    -- Host organism.

    left outer join xml_xref_experiment_organisms as taxidE
        on (I.source, I.filename, I.entry, E.experimentid) = (taxidE.source, taxidE.filename, taxidE.entry, taxidE.experimentid)
    left outer join taxonomy_names as taxnamesE
        on taxidE.taxid = taxnamesE.taxid
        and taxnamesE.nameclass = 'scientific name'

    -- Interaction detection method.

    left outer join tmp_methods as methodE
        on (I.source, I.filename, I.entry, E.experimentid) = (methodE.source, methodE.filename, methodE.entry, methodE.experimentid)
        and methodE.property = 'interactionDetectionMethod'

    -- Interaction type.

    left outer join xml_xref_interaction_types as typeI
        on (I.source, I.filename, I.entry, I.interactionid) = (typeI.source, typeI.filename, typeI.entry, typeI.interactionid)
    left outer join psicv_terms as typenameI
        on typeI.refvalue = typenameI.code
        and typenameI.nametype = 'preferred'

    -- PubMed identifiers.

    left outer join tmp_pubmed as pubmedE
        on (I.source, I.filename, I.entry, E.experimentid) = (pubmedE.source, pubmedE.filename, pubmedE.entry, pubmedE.experimentid)

    -- Authors.

    left outer join xml_names_experiment_authors as authorE
        on (I.source, I.filename, I.entry, E.experimentid) = (authorE.source, authorE.filename, authorE.entry, authorE.experimentid)

    -- IMEX identifier.

    left outer join xml_xref_interactions as imexI
        on (I.source, I.filename, I.entry, I.interactionid) = (imexI.source, imexI.filename, imexI.entry, imexI.interactionid)
        and imexI.dblabel = 'imex';

analyze tmp_interaction_experiments;

-- Combine interactor information.

create temporary table tmp_interactor_experiments as
    select I.*,

        -- finalReferenceA (the original reference for A, or a corrected/complete/updated/unambiguous reference)
        -- NOTE: This actually appears as "-" in the iRefIndex 9 MITAB output for complexes.

        case when edgetype = 'C' then 'complex'
             else nameA.finaldblabel
        end as finaldblabelA,

        case when edgetype = 'C' then I.rigid
             else nameA.finalrefvalue
        end as finalrefvalueA,

        -- finalReferenceB (the original reference for B, or a corrected/complete/updated/unambiguous reference)

        nameB.finaldblabel as finaldblabelB,
        nameB.finalrefvalue as finalrefvalueB,

        -- originalReferenceA (original primary or secondary reference for A, the rigid of any complex as 'complex:...')
        -- NOTE: This actually appears as "-" in the iRefIndex 9 MITAB output for complexes.

        case when edgetype = 'C' then 'complex:' || I.rigid
             else nameA.originaldblabel || ':' || nameA.originalrefvalue
        end as originalReferenceA,

        -- originalReferenceB (original primary or secondary reference for B)

        nameB.originaldblabel || ':' || nameB.originalrefvalue as originalReferenceB,

        -- taxA (as "taxid:...")

        case when edgetype = 'C' or nameA.taxid is null then '-'
             else 'taxid:' || nameA.taxid || '(' || coalesce(taxnamesA.name, '-') || ')'
        end as taxA,

        -- taxB (as "taxid:...")

        case when nameB.taxid is null then '-'
             else 'taxid:' || nameB.taxid || '(' || coalesce(taxnamesB.name, '-') || ')'
        end as taxB,

        -- mappingScoreA (operation characters describing the original-to-final transformation, "-" for complexes)

        case when edgetype = 'C' then '-' else nameA.score end as mappingScoreA,

        -- mappingScoreB (operation characters describing the original-to-final transformation)

        nameB.score as mappingScoreB

    from tmp_interaction_experiments as I

    -- Information for interactor A.

    left outer join tmp_assignments as nameA
        on (I.source, I.filename, I.entry, I.interactoridA) = 
           (nameA.source, nameA.filename, nameA.entry, nameA.interactorid)
        and I.edgetype <> 'C'
    left outer join taxonomy_names as taxnamesA
        on nameA.taxid = taxnamesA.taxid
        and taxnamesA.nameclass = 'scientific name'

    -- Information for interactor B.

    inner join tmp_assignments as nameB
        on (I.source, I.filename, I.entry, I.interactoridB) = 
           (nameB.source, nameB.filename, nameB.entry, nameB.interactorid)
    left outer join taxonomy_names as taxnamesB
        on nameB.taxid = taxnamesB.taxid
        and taxnamesB.nameclass = 'scientific name';

analyze tmp_interactor_experiments;

-- Collect all participant-, interactor- and interaction-related information.

create temporary table tmp_mitab_interactions as
    select I.*,

        -- biologicalRoleA

        case when edgetype = 'C' or bioroleA.refvalues is null or array_length(bioroleA.refvalues, 1) = 0 then 'MI:0000(unspecified)'
             else array_to_string(bioroleA.refvalues, '|')
        end as biologicalRoleA,

        -- biologicalRoleB

        case when bioroleB.refvalues is null or array_length(bioroleB.refvalues, 1) = 0 then 'MI:0000(unspecified)'
             else array_to_string(bioroleB.refvalues, '|')
        end as biologicalRoleB,

        -- experimentalRoleA

        case when edgetype = 'C' or exproleA.refvalues is null or array_length(exproleA.refvalues, 1) = 0 then 'MI:0000(unspecified)'
             else array_to_string(exproleA.refvalues, '|')
        end as experimentalRoleA,

        -- experimentalRoleB

        case when exproleB.refvalues is null or array_length(exproleB.refvalues, 1) = 0 then 'MI:0000(unspecified)'
             else array_to_string(exproleB.refvalues, '|')
        end as experimentalRoleB

    from tmp_interactor_experiments as I

    -- Information for participant A.

    left outer join tmp_participants as bioroleA
        on (I.source, I.filename, I.entry, I.participantidA) = (bioroleA.source, bioroleA.filename, bioroleA.entry, bioroleA.participantid)
        and bioroleA.property = 'biologicalRole'
        and I.edgetype <> 'C'
    left outer join tmp_participants as exproleA
        on (I.source, I.filename, I.entry, I.participantidA) = (exproleA.source, exproleA.filename, exproleA.entry, exproleA.participantid)
        and exproleA.property = 'experimentalRole'
        and I.edgetype <> 'C'

    -- Information for participant B.

    left outer join tmp_participants as bioroleB
        on (I.source, I.filename, I.entry, I.participantidB) = (bioroleB.source, bioroleB.filename, bioroleB.entry, bioroleB.participantid)
        and bioroleB.property = 'biologicalRole'
    left outer join tmp_participants as exproleB
        on (I.source, I.filename, I.entry, I.participantidB) = (exproleB.source, exproleB.filename, exproleB.entry, exproleB.participantid)
        and exproleB.property = 'experimentalRole';

analyze tmp_mitab_interactions;

-- Combine with ROG-related information to produce MITAB-appropriate records.

create temporary table tmp_mitab_all as
    select

        -- uidA (identifier, preferably uniprotkb accession, refseq, complex as 'complex:...')

        case when edgetype = 'C' then 'complex:' || I.rigid
             else prefA.dblabel || ':' || prefA.refvalue
        end as uidA,

        -- uidB (identifier, preferably uniprotkb accession, refseq)

        prefB.dblabel || ':' || prefB.refvalue as uidB,

        -- altA (alternatives for A, preferably uniprotkb accession, refseq, entrezgene/locuslink identifier, including rogid, irogid)

        case when edgetype = 'C' then 'rogid:' || I.rigid
             else array_to_string(
                array_cat(
                    rognameA.names[1:3],
                    array[
                        'rogid:' || I.uidA,
                        'irogid:' || cast(
                            case when edgetype = 'C' then irig.rig
                                 else irogA.rog
                            end as varchar)
                        ]
                    ), '|')
        end as altA,

        -- altB (alternatives for B, preferably uniprotkb accession, refseq, entrezgene/locuslink identifier, including rogid, irogid)

        array_to_string(
            array_cat(
                rognameB.names[1:3],
                array[
                    'rogid:' || I.uidB,
                    'irogid:' || cast(irogB.rog as varchar)
                    ]
                ), '|') as altB,

        -- aliasA (aliases for A, preferably uniprotkb identifier/entry, entrezgene/locuslink symbol, including crogid, icrogid)

        case when edgetype = 'C' then 'crogid:' || crigid.crigid || '|icrogid:' || cast(icrig.rig as varchar)
             else array_to_string(
                array_cat(
                    case when aliasA.aliases is null or array_length(aliasA.aliases, 1) = 0 then cast(array[] as varchar[])
                         else aliasA.aliases[1:6]
                    end,
                    array[
                        'crogid:' || crogidA.crogid,
                        'icrogid:' || cast(icrogA.rog as varchar)
                        ]
                    ), '|')
        end as aliasA,

        -- aliasB (aliases for B, preferably uniprotkb identifier/entry, entrezgene/locuslink symbol, including crogid, icrogid)

        array_to_string(
            array_cat(
                case when aliasB.aliases is null or array_length(aliasB.aliases, 1) = 0 then cast(array[] as varchar[])
                     else aliasB.aliases[1:6]
                end,
                array[
                    'crogid:' || crogidB.crogid,
                    'icrogid:' || icrogB.rog
                    ]
                ), '|') as aliasB,

        -- method (interaction detection method as "MI:code(name)")

        method,

        -- authors (as "name-[year[-number]]")

        authors,

        -- pmids (as "pubmed:...")

        pmids,

        -- taxA (as "taxid:...")

        taxA,

        -- taxB (as "taxid:...")

        taxB,

        -- interactionType (interaction type as "MI:code(name)")

        interactionType,

        -- sourcedb (as "MI:code(name)" using "MI:0000(name)" for non-CV sources)

        sourcedb,

        -- interactionIdentifier (includes rigid, irigid, and edgetype as "X", "Y" or "C",
        -- with "source:-" if no database-specific identifier is provided)

        interactionIdentifier,

        -- confidence

        case when confidence = '' then '-' else confidence end as confidence,

        -- expansion

        case when edgetype = 'C' then 'bipartite' else 'none' end as expansion,

        -- biologicalRoleA

        biologicalRoleA,

        -- biologicalRoleB

        biologicalRoleB,

        -- experimentalRoleA

        experimentalRoleA,

        -- experimentalRoleB

        experimentalRoleB,

        -- interactorTypeA

        case when edgetype = 'C' then 'MI:0315(protein complex)' else 'MI:0326(protein)' end as interactorTypeA,

        -- interactorTypeB

        cast('MI:0326(protein)' as varchar) as interactorTypeB,

        -- xrefsA (always "-")

        cast('-' as varchar) as xrefsA,

        -- xrefsB (always "-")

        cast('-' as varchar) as xrefsB,

        -- xrefsInteraction (always "-")

        cast('-' as varchar) as xrefsInteraction,

        -- annotationsA (always "-")

        cast('-' as varchar) as annotationsA,

        -- annotationsB (always "-")

        cast('-' as varchar) as annotationsB,

        -- annotationsInteraction (always "-")

        cast('-' as varchar) as annotationsInteraction,

        -- hostOrganismTaxid (as "taxid:...")

        hostOrganismTaxid,

        -- parametersInteraction (always "-")

        cast('-' as varchar) as parametersInteraction,

        -- creationDate (the iRefIndex release date as "YYYY/MM/DD")

        cast(current_date as varchar) as creationDate,

        -- updateDate (the iRefIndex release date as "YYYY/MM/DD")

        cast(current_date as varchar) as updateDate,

        -- checksumA (the rogid for interactor A as "rogid:...")

        'rogid:' || I.uidA as checksumA,

        -- checksumB (the rogid for interactor B as "rogid:...")

        'rogid:' || I.uidB as checksumB,

        -- checksumInteraction (the rigid for the interaction as "rigid:...")

        'rigid:' || I.rigid as checksumInteraction,

        -- negative (always "false")

        cast('false' as varchar) as negative,

        -- originalReferenceA (original primary or secondary reference for A, the rigid of any complex as 'complex:...')
        -- NOTE: This actually appears as "-" in the iRefIndex 9 MITAB output for complexes.

        originalReferenceA,

        -- originalReferenceB (original primary or secondary reference for B)

        originalReferenceB,

        -- finalReferenceA (the original reference for A, or a corrected/complete/updated/unambiguous reference)
        -- NOTE: This actually appears as "-" in the iRefIndex 9 MITAB output for complexes.

        finaldblabelA || ':' || finalrefvalueA as finalReferenceA,

        -- finalReferenceB (the original reference for B, or a corrected/complete/updated/unambiguous reference)

        finaldblabelB || ':' || finalrefvalueB as finalReferenceB,

        -- mappingScoreA (operation characters describing the original-to-final transformation, "-" for complexes)

        mappingScoreA,

        -- mappingScoreB (operation characters describing the original-to-final transformation)

        mappingScoreB,

        -- irogidA (the integer identifier for the rogid for A)

        cast(case when edgetype = 'C' then irig.rig
                  else irogA.rog
             end as varchar) as irogidA,

        -- irogidB (the integer identifier for the rogid for B)

        cast(irogB.rog as varchar) as irogidB,

        -- irigid (the integer identifier for the rigid for the interaction)

        cast(irig.rig as varchar) as irigid,

        -- crogidA (the canonical rogid for A, not prefixed)

        case when edgetype = 'C' then crigid.crigid
             else crogidA.crogid
        end as crogidA,

        -- crogidB (the canonical rogid for B, not prefixed)

        crogidB.crogid as crogidB,

        -- crigid (the canonical rigid for the interaction, not prefixed)

        crigid.crigid as crigid,

        -- icrogidA (the integer identifier for the canonical rogid for A)

        cast(case when edgetype = 'C' then icrig.rig
                  else icrogA.rog
             end as varchar) as icrogidA,

        -- icrogidB (the integer identifier for the canonical rogid for B)

        cast(icrogB.rog as varchar) as icrogidB,

        -- icrigid (the integer identifier for the canonical rigid for the interaction)

        cast(icrig.rig as varchar) as icrigid,

        -- imexid (as "imex:..." or "-" if not available)

        imexid,

        -- edgetype (as "X", "Y" or "C")

        I.edgetype,

        -- numParticipants (the number of participants)

        cast(I.numParticipants as varchar) as numParticipants

    from tmp_mitab_interactions as I
    left outer join tmp_identifiers as rognameA
        on I.uidA = rognameA.rogid
        and I.edgetype <> 'C'
    left outer join irefindex_rogid_identifiers_preferred as prefA
        on I.uidA = prefA.rogid
        and I.edgetype <> 'C'
    inner join tmp_identifiers as rognameB
        on I.uidB = rognameB.rogid
    inner join irefindex_rogid_identifiers_preferred as prefB
        on I.uidB = prefB.rogid
    left outer join tmp_aliases as aliasA
        on I.uidA = aliasA.rogid
        and I.edgetype <> 'C'
    left outer join tmp_aliases as aliasB
        on I.uidB = aliasB.rogid

    -- Incorporate canonical information.

    left outer join irefindex_rogids_canonical as crogidA
        on I.uidA = crogidA.rogid
        and I.edgetype <> 'C'
    inner join irefindex_rogids_canonical as crogidB
        on I.uidB = crogidB.rogid
    inner join irefindex_rigids_canonical as crigid
        on I.rigid = crigid.rigid

    -- Incorporate integer identifiers.

    left outer join irefindex_rog2rogid as irogA
        on I.uidA = irogA.rogid
        and I.edgetype <> 'C'
    inner join irefindex_rog2rogid as irogB
        on I.uidB = irogB.rogid
    inner join irefindex_rig2rigid as irig
        on I.rigid = irig.rigid

    -- Incorporate integer identifiers for canonical identifiers.
    -- These identifiers are a subset of the general integer numbering schemes
    -- for interactions and interactors.

    left outer join irefindex_rog2rogid as icrogA
        on crogidA.crogid = icrogA.rogid
        and I.edgetype <> 'C'
    inner join irefindex_rog2rogid as icrogB
        on crogidB.crogid = icrogB.rogid
    inner join irefindex_rig2rigid as icrig
        on crigid.crigid = icrig.rigid;

analyze tmp_mitab_all;

-- Final output production.

-- The complete data set for all organisms. Note that the headers will need
-- adding.

\copy tmp_mitab_all to '<directory>/mitab_all'

-- Make special organism-specific files.

create temporary table tmp_organisms (
    taxid integer not null
);

\copy tmp_organisms from '<directory>/organisms.txt'
analyze tmp_organisms;

create temporary table tmp_interaction_taxids as
    select distinct 'rigid:' || rigid as checksumInteraction, cast(substring(rogid from 28) as integer) as taxid
    from irefindex_interactions;

analyze tmp_interaction_taxids;

-- The concatenated data set for each of the selected organisms. Note that the
-- first taxonomy column will need stripping and the headers adding.

create temporary table tmp_mitab_all_organisms as
    select coalesce(cast(O.taxid as varchar), 'other') as taxid, M.*
    from tmp_mitab_all as M
    inner join tmp_interaction_taxids as T
        on T.checksumInteraction = M.checksumInteraction
    left outer join tmp_organisms as O
        on O.taxid = T.taxid
    order by T.taxid;

\copy tmp_mitab_all_organisms to '<directory>/mitab_all_organisms'

rollback;
