CREATE SCHEMA phylogeny;

SET search_path = phylogeny;

-- reference tree

CREATE TABLE phylogeny.species_tree (
  branch_id SMALLINT PRIMARY KEY,
  parent_branch_id SMALLINT,
  branch_name VARCHAR(50) DEFAULT NULL,   -- only for tips and clade/populationancestors
  is_tip BOOL 
);

CREATE TYPE phylogeny.eventtype AS ENUM ('D', 'T', 'S', 'SL', 'L', 'O');

CREATE TABLE phylogeny.species_tree_events (
  event_id INT PRIMARY KEY,
  event_type eventtype NOT NULL,
  don_branch_id SMALLINT,          -- refers to species_tree (branch_id)
  rec_branch_id SMALLINT NOT NULL  -- refers to species_tree (branch_id)
);

CREATE INDEX ON species_tree_events (event_type);
CREATE INDEX ON species_tree_events (don_branch_id);
CREATE INDEX ON species_tree_events (rec_branch_id);

CREATE TABLE phylogeny.gene_lineage_events ( --to be a large table
  event_id INT,
  replacement_label_or_cds_code VARCHAR(50) NOT NULL, -- refers to genome.coding_sequences (cds_code) and phylogeny.replaced_gene_tree_clades (replacement_label)
  freq INT NOT NULL,
  reconciliation_id SMALLINT DEFAULT NULL    -- to distinguish reconciliation sets; can be NULL if not to be redundant
);

CREATE TABLE phylogeny.reconciliation_collections (
  reconciliation_id SMALLINT NOT NULL,
  reconciliation_name VARCHAR NOT NULL,
  software VARCHAR NOT NULL,
  version VARCHAR NOT NULL,
  algorithm VARCHAR,
  reconciliation_date TIMESTAMP,
  notes TEXT
);

-- gene trees

CREATE TABLE phylogeny.criteria_collapse_gene_tree_clades (
  criterion_id INT PRIMARY KEY,
  criterion_name VARCHAR(50) NOT NULL,
  criterion_definition TEXT,
  collapsed_clade_collection_creation DATE
);

CREATE TABLE phylogeny.collapsed_gene_tree_clades (
  gene_family_id VARCHAR(20) NOT NULL,
  col_clade VARCHAR(10) NOT NULL,
  cds_code VARCHAR(20) NOT NULL,
  collapse_criterion_id INT DEFAULT NULL
);

CREATE TABLE phylogeny.criteria_replace_gene_tree_clades (
  criterion_id INT PRIMARY KEY,
  criterion_name VARCHAR(50) NOT NULL,
  criterion_definition TEXT,
  replaced_clade_collection_creation DATE
);

CREATE TABLE phylogeny.replaced_gene_tree_clades (
  gene_family_id VARCHAR(20) NOT NULL,
  col_clade_or_cds_code VARCHAR(20) NOT NULL,
  replacement_label VARCHAR(60) DEFAULT NULL,
  replace_criterion_id INT DEFAULT NULL
);

-- after filling the tables

-- ~ CREATE INDEX ON gene_lineage_events (reconciliation_id);
-- ~ CREATE INDEX ON gene_lineage_events (cds_code);
-- ~ CREATE INDEX ON gene_lineage_events (event_id);
-- ~ CREATE INDEX ON gene_lineage_events (freq);
-- ~ ALTER TABLE gene_lineage_events ADD PRIMARY KEY (event_id, cds_code, reconciliation_id);

-- ~ CREATE INDEX ON collapsed_gene_tree_clades (gene_family_id);
-- ~ CREATE INDEX ON collapsed_gene_tree_clades (cds_code);
-- ~ CREATE INDEX ON collapsed_gene_tree_clades (gene_family_id, col_clade);

-- ~ CREATE INDEX ON replaced_gene_tree_clades (gene_family_id, col_clade_or_cds_code);
-- ~ CREATE INDEX ON replaced_gene_tree_clades (replacement_label);

-- ~ CREATE TABLE replacement_label_or_cds_code2gene_families AS (
-- ~ SELECT cds_code as replacement_label_or_cds_code, gene_family_id FROM coding_sequences
-- ~ UNION
-- ~ SELECT replacement_label as replacement_label_or_cds_code, gene_family_id FROM replaced_gene_tree_clades
-- ~ );
-- ~ CREATE INDEX ON replacement_label_or_cds_code2gene_families (replacement_label_or_cds_code);
-- ~ CREATE INDEX ON replacement_label_or_cds_code2gene_families (gene_family_id);

-- ~ CREATE TABLE gene_tree_label2cds_code (replacement_label_or_cds_code, cds_code) AS
  -- ~ SELECT replacement_label_or_cds_code, replacement_label_or_cds_code as cds_code
   -- ~ FROM gene_lineage_events
   -- ~ WHERE char_length(replacement_label_or_cds_code) < 20
 -- ~ UNION 
  -- ~ SELECT rgtc1.replacement_label as replacement_label_or_cds_code, rgtc1.col_clade_or_cds_code as cds_code
   -- ~ FROM replaced_gene_tree_clades AS rgtc1
   -- ~ WHERE rgtc1.col_clade_or_cds_code NOT LIKE 'clade%'
 -- ~ UNION
  -- ~ SELECT rgtc2.replacement_label as replacement_label_or_cds_code, cgtc.cds_code
   -- ~ FROM replaced_gene_tree_clades AS rgtc2
   -- ~ INNER JOIN collapsed_gene_tree_clades AS cgtc ON rgtc2.col_clade_or_cds_code=cgtc.col_clade AND rgtc2.gene_family_id=cgtc.gene_family_id
-- ~ ;
-- ~ CREATE INDEX ON gene_tree_label2cds_code (replacement_label_or_cds_code);
-- ~ CREATE INDEX ON gene_tree_label2cds_code (cds_code);
