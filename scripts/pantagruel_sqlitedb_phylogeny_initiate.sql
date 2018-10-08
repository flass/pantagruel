-- reference tree

CREATE TABLE   species_tree (
  branch_id INT PRIMARY KEY,
  parent_branch_id INT,
  branch_name VARCHAR(50) DEFAULT NULL,   -- only for tips and clade/populationancestors
  is_tip BOOL 
);

CREATE TABLE   species_tree_events (
  event_id INT PRIMARY KEY,
  event_type char(1) NOT NULL,
  don_branch_id INT,          -- refers to species_tree (branch_id)
  rec_branch_id INT NOT NULL  -- refers to species_tree (branch_id)
);

CREATE INDEX species_tree_events_ev_type_idx ON species_tree_events (event_type);
CREATE INDEX species_tree_events_don_br_idx ON species_tree_events (don_branch_id);
CREATE INDEX species_tree_events_rec_br_idx ON species_tree_events (rec_branch_id);

CREATE TABLE gene_lineage_events ( --to be a large table
  event_id SERIAL,
  replacement_label_or_cds_code VARCHAR(60) NOT NULL,    -- refers to genome.coding_sequences (cds_code) and phylogeny.replaced_gene_tree_clades (replacement_label)
  freq INT NOT NULL,
  reconciliation_id INT NOT NULL DEFAULT 0    -- to distinguish reconciliation sets; can be left to default if not to be redundant
);

CREATE TABLE reconciliation_collections (
  reconciliation_id INT NOT NULL,
  reconciliation_name VARCHAR NOT NULL,
  software VARCHAR NOT NULL,
  version VARCHAR NOT NULL,
  algorithm VARCHAR,
  reconciliation_date TIMESTAMP,
  notes TEXT
);

-- gene trees

CREATE TABLE criteria_collapse_gene_tree_clades (
  criterion_id INT PRIMARY KEY,
  criterion_name VARCHAR(50) NOT NULL,
  criterion_definition TEXT,
  collapsed_clade_collection_creation DATE
);

CREATE TABLE collapsed_gene_tree_clades (
  gene_family_id VARCHAR(20) NOT NULL,
  col_clade VARCHAR(10) NOT NULL,
  cds_code VARCHAR(20) NOT NULL,
  collapse_criterion_id INT DEFAULT NULL
);

CREATE TABLE criteria_replace_gene_tree_clades (
  criterion_id INT PRIMARY KEY,
  criterion_name VARCHAR(50) NOT NULL,
  criterion_definition TEXT,
  replaced_clade_collection_creation DATE
);

CREATE TABLE replaced_gene_tree_clades (
  gene_family_id VARCHAR(20) NOT NULL,
  col_clade_or_cds_code VARCHAR(20) NOT NULL,
  replacement_label VARCHAR(60) DEFAULT NULL,
  replace_criterion_id INT DEFAULT NULL
);

CREATE VIEW replacement_label_or_cds_code2gene_families AS 
SELECT cds_code as replacement_label_or_cds_code, gene_family_id FROM coding_sequences
UNION
SELECT replacement_label as replacement_label_or_cds_code, gene_family_id FROM replaced_gene_tree_clades;


CREATE TABLE ortholog_collections (
  ortholog_col_id INT PRIMARY KEY,
  ortholog_col_name VARCHAR(50) NOT NULL,
  reconciliation_id INT NOT NULL,
  software VARCHAR NOT NULL,
  version VARCHAR NOT NULL,
  algorithm VARCHAR,
  ortholog_col_date TIMESTAMP,
  notes TEXT
);

CREATE TABLE orthologous_groups (
  cds_code VARCHAR(20) NOT NULL,
  gene_family_id VARCHAR(20) NOT NULL,
  og_id INT NOT NULL,
  ortholog_col_id INT
);
