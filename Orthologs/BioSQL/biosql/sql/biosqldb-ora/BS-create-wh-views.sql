--
-- SQL script to create the views for SYMGENE/BioSQL that utilize the
-- warehouse materialized views.
--
-- $GNF: projects/gi/symgene/src/DB/BS-create-wh-views.sql,v 1.5 2003/06/10 20:06:30 hlapp Exp $
--

--
-- Copyright 2002-2003 Genomics Institute of the Novartis Research Foundation
-- Copyright 2002-2008 Hilmar Lapp
-- 
--  This file is part of BioSQL.
--
--  BioSQL is free software: you can redistribute it and/or modify it
--  under the terms of the GNU Lesser General Public License as
--  published by the Free Software Foundation, either version 3 of the
--  License, or (at your option) any later version.
--
--  BioSQL is distributed in the hope that it will be useful,
--  but WITHOUT ANY WARRANTY; without even the implied warranty of
--  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
--  GNU Lesser General Public License for more details.
--
--  You should have received a copy of the GNU Lesser General Public License
--  along with BioSQL. If not, see <http://www.gnu.org/licenses/>.
--

--
-- Genome mappings
--
PROMPT
PROMPT Creating view SG_Ent_Chr_Maps

CREATE OR REPLACE VIEW SG_Ent_Chr_Maps
AS
SELECT
	ECM.EntSeg_Start_Pos	EntSeg_Start_Pos
	, ECM.EntSeg_End_Pos	EntSeg_End_Pos
	, ECM.EntSeg_Num	EntSeg_Num
	, FType.Name		EntSeg_Type_Name
	, FSrc.Name		EntSeg_Source_Name
	, ECM.ChrSeg_Start_Pos	ChrSeg_Start_Pos
	, ECM.ChrSeg_End_Pos	ChrSeg_End_Pos
	, ECM.ChrSeg_Strand	ChrSeg_Strand
	, ECM.ChrSeg_Pct_Identity ChrSeg_Pct_Identity
	, Ent.Accession		Ent_Accession
	, Ent.Identifier	Ent_Identifier
	, Ent.Description	Ent_Description
	, Ent.Version		Ent_Version
	, DB.Name		DB_Name
	, DB.Acronym		DB_Acronym
	, EntTax.Name		Ent_Tax_Name
	, EntTNod.NCBI_Taxon_ID	Ent_Tax_NCBI_Taxon_ID
	, Chr.Name		Chr_Name
	, Chr.Accession		Chr_Accession
	, Asm.Name		Asm_Name
	, Asm.Acronym		Asm_Acronym
	, ChrTax.Name		Chr_Tax_Name
	, ChrTNod.NCBI_Taxon_ID Chr_Tax_NCBI_Taxon_ID
	, ECM.EntSeg_Loc_Oid	EntSeg_Loc_Oid
	, ECM.EntSeg_Oid	EntSeg_Oid
	, ECM.Ent_Oid		Ent_Oid
	, Ent.Tax_Oid		Ent_Tax_Oid
	, Ent.DB_Oid		DB_Oid
	, ECM.ChrSeg_Loc_Oid	ChrSeg_Loc_Oid
	, ECM.ChrSeg_Oid	ChrSeg_Oid
	, ECM.Chr_Oid		Chr_Oid
	, Chr.Tax_Oid		Chr_Tax_Oid
	, Chr.DB_Oid		Asm_Oid
	, ECM.EntSeg_Type_Oid	EntSeg_Type_Oid
	, ECM.EntSeg_Source_Oid	EntSeg_Source_Oid
FROM SG_Ent_Chr_Map ECM, 
     SG_Bioentry Ent, SG_Bioentry Chr,
     SG_Taxon EntTNod, SG_Taxon_Name EntTax, 
     SG_Taxon ChrTNod, SG_Taxon_Name ChrTax, 
     SG_Biodatabase DB, SG_Biodatabase Asm,
     SG_Term FType, SG_Term FSrc
WHERE
     ECM.Ent_Oid	 = Ent.Oid
AND  ENT.DB_Oid		 = DB.Oid
AND  ENT.Tax_Oid	 = EntTNod.Oid (+)
AND  EntTax.Tax_Oid (+)	 = EntTNod.Oid
AND  ECM.Chr_Oid	 = Chr.Oid
AND  Chr.DB_Oid	 	 = Asm.Oid
AND  Chr.Tax_Oid	 = ChrTNod.Oid
AND  ChrTax.Tax_Oid	 = ChrTNod.Oid
AND  ECM.EntSeg_Type_Oid = FType.Oid
AND  ECM.EntSeg_Source_Oid = FSrc.Oid
;

--
-- Single entry point name searching
--
PROMPT
PROMPT Creating view SG_Bioentry_Names

CREATE OR REPLACE VIEW SG_Bioentry_Names
AS
SELECT
	BEN.Ent_Name		Ent_Name_Key
	, Ent.Accession		Ent_Accession
	, Ent.Identifier	Ent_Identifier
	, Ent.Version		Ent_Version
	, Ent.Name		Ent_Name
	, Ent.Description	Ent_Description
	, DB.Name		DB_Name
	, Tax.Tax_Name		Tax_Name
	, Tax.Tax_Common_Name	Tax_Common_Name
	, Tax.Tax_NCBI_Taxon_ID	Tax_NCBI_Taxon_ID
	, BEN.Ent_Oid		Ent_Oid
	, Ent.DB_Oid		DB_Oid
	, Ent.Tax_Oid		Tax_Oid
-- legacy mappings
	, Ent.Name		Ent_Display_ID
FROM SG_Bioentry_Name BEN, SG_Bioentry Ent,
     -- we take advantage of the taxa view here
     SG_Taxa Tax, SG_Biodatabase DB
WHERE
     BEN.Ent_Oid     = Ent.Oid
AND  Ent.DB_Oid      = DB.Oid
AND  Ent.Tax_Oid     = Tax.Tax_Oid
;
