--
-- This SQL script is to add tables to the pre-Singapore Oracle version of the
-- BioSQL schema that were added in the so-called Singapore version.
--
-- Disclaimer: This script and scripts it launches will modify the schema. It
-- will modify and drop tables, indexes, column names, triggers, sequences,
-- and possibly more. YOU SHOULD BACKUP YOUR DATABASE BEFORE RUNNING THIS
-- SCRIPT, or any of the scripts it launches. You should also verify that
-- you can actually restore from your backup. If you fail to properly take
-- a backup that restores the database to its prior status, serious
-- consequences up to and including complete loss of your data may result
-- from the operation of this script. THIS PACKAGE IS PROVIDED WITHOUT ANY
-- WARRANTIES WHATSOEVER. Please read the license under which you may use
-- this script and those that come with it.
--
-- $GNF: projects/gi/symgene/src/sql/migrate/singapore/migrate-add-tables.sql,v 1.2 2003/05/21 06:47:24 hlapp Exp $
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

PROMPT Adding new tables
--
-- Ontology term to DBXref associations
--
PROMPT     - Term DBXref associations

CREATE TABLE SG_Term_DBXRef_Assoc ( 
	Trm_Oid		INTEGER NOT NULL, 
	DBX_Oid		INTEGER NOT NULL , 
	Rank 		NUMBER(3), 
	CONSTRAINT XPKTerm_DBXRef_Assoc
		PRIMARY KEY (Trm_Oid , DBX_Oid)
	USING INDEX
	TABLESPACE &biosql_index
	--
)
;

CREATE INDEX XIF1Term_DBXRef_Assoc ON SG_Term_DBXRef_Assoc
(
	DBX_Oid
)
    	 TABLESPACE &biosql_index
;

ALTER TABLE SG_Term_DBXRef_Assoc
       ADD  ( CONSTRAINT FKTrm_TrmDBXA
              FOREIGN KEY (Trm_Oid)
                             REFERENCES SG_Term (Oid)
			     ON DELETE CASCADE ) ;


ALTER TABLE SG_Term_DBXRef_Assoc
       ADD  ( CONSTRAINT FKDBX_TrmDBXA
              FOREIGN KEY (DBX_Oid)
                             REFERENCES SG_DBXRef (Oid)
			     ON DELETE CASCADE ) ;


--
-- Ontology Term synonyms
--
PROMPT     - Ontology Term Synonyms

CREATE TABLE SG_Term_Synonym (
	Name			VARCHAR2(128) NOT NULL,
	Trm_Oid			INTEGER NOT NULL,
       	CONSTRAINT XPKTerm_Synonym
		PRIMARY KEY (Name, Trm_Oid)
	USING INDEX
    	TABLESPACE &biosql_index
	-- 
)
;

ALTER TABLE SG_Term_Synonym
       ADD  ( CONSTRAINT FKTrm_TSyn
              FOREIGN KEY (Trm_Oid)
                             REFERENCES SG_Term (Oid)
			     ON DELETE CASCADE ) ;


--
-- DBXref to Ontology Term associations (qualifier/value pairs)
--
PROMPT     - DBXref qualifier associations

CREATE TABLE SG_DBXRef_Qualifier_Assoc ( 
	DBX_Oid			INTEGER NOT NULL, 
	Trm_Oid			INTEGER NOT NULL, 
	Rank			NUMBER(3) DEFAULT 0 NOT NULL, 
	Value			VARCHAR2(256), 
	CONSTRAINT XPKDBXRef_Qualifier_Assoc
		PRIMARY KEY (DBX_Oid, Trm_Oid, Rank)
	USING INDEX
    	TABLESPACE &biosql_index
	--
)
; 

CREATE INDEX XIF1DBXRef_Qualifier_Assoc ON SG_DBXRef_Qualifier_Assoc
(
	Trm_Oid
)
    	TABLESPACE &biosql_index
; 

ALTER TABLE SG_DBXRef_Qualifier_Assoc
       ADD  ( CONSTRAINT FKDBX_DBXTrmA
              FOREIGN KEY (DBX_Oid)
                             REFERENCES SG_DBXRef (Oid) 
                             ON DELETE CASCADE ) ;

ALTER TABLE SG_DBXRef_Qualifier_Assoc
       ADD  ( CONSTRAINT FKTrm_DBXTrmA
              FOREIGN KEY (Trm_Oid)
                             REFERENCES SG_Term (Oid) 
                             ON DELETE CASCADE ) ;


--
-- Seqfeature to DBXref associations
--
PROMPT     - Seqfeature to DBXref associations

CREATE TABLE SG_Seqfeature_DBXref_Assoc ( 
	Fea_Oid			INTEGER NOT NULL, 
	DBX_Oid			INTEGER NOT NULL, 
	Rank			NUMBER(3), 
	CONSTRAINT XPKSeqfeature_DBXref_Assoc
		PRIMARY KEY (Fea_Oid, DBX_Oid)
	USING INDEX
    	TABLESPACE &biosql_index
	--
)
; 

CREATE INDEX XIF1Seqfeature_DBXref_Assoc On SG_Seqfeature_DBXref_Assoc
(
	DBX_Oid
)
    	TABLESPACE &biosql_index
; 

ALTER TABLE SG_Seqfeature_DBXref_Assoc
       ADD  ( CONSTRAINT FKFea_DbxFeaA
              FOREIGN KEY (Fea_Oid)
                             REFERENCES SG_Seqfeature (Oid) 
                             ON DELETE CASCADE ) ;


ALTER TABLE SG_Seqfeature_DBXRef_Assoc
       ADD  ( CONSTRAINT FKDBX_DbxFeaA
              FOREIGN KEY (DBX_Oid)
                             REFERENCES SG_DBXRef (Oid)
                             ON DELETE CASCADE ) ;


PROMPT Done with new tables.
