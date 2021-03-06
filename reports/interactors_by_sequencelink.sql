-- Show interactor origins.

-- Copyright (C) 2011 Ian Donaldson <ian.donaldson@biotek.uio.no>
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

-- Show the number of distinct interactors by connection to sequence database.

create temporary table tmp_interactors_by_sequencelink as
    select dblabel, sequencelink, count(distinct refvalue)
    from xml_xref_sequences
    group by dblabel, sequencelink
    order by dblabel, sequencelink;

\copy tmp_interactors_by_sequencelink to '<directory>/interactors_by_sequencelink'

rollback;
