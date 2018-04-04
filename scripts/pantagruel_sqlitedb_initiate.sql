CREATE TABLE assemblies (
        assembly_id CHAR(16) PRIMARY KEY,
        assembly_name VARCHAR(50) NOT NULL,
        organism VARCHAR(300),
		species VARCHAR(200),
		subspecies VARCHAR(200),
		serovar VARCHAR(200),
		strain VARCHAR(50),
		taxid INTEGER NOT NULL,
		primary_pubmed_id INTEGER,
		country VARCHAR(200),
		isolation_source VARCHAR(200),
		host VARCHAR(200),
		clinical_source VARCHAR(200),
		collection_year INTEGER,
		collection_month INTEGER,
		collection_day INTEGER,
		sequencing_technology VARCHAR(200),
		sequencing_coverage VARCHAR(50),
		note TEXT
);

CREATE TABLE replicons (
        genomic_accession CHAR(14) PRIMARY KEY,
        replicon_name VARCHAR(500),
        replicon_size INTEGER NOT NULL,
        replicon_type CHAR(16) NOT NULL,
        assembly_id CHAR(16) NOT NULL
);

CREATE TABLE coding_sequences (
        cds_serial_id SERIAL PRIMARY KEY,
        genbank_cds_id VARCHAR(50) UNIQUE NOT NULL,
        genomic_accession CHAR(14) NOT NULL,
        locus_tag VARCHAR(200),
        cds_begin INTEGER NOT NULL,
        cds_end INTEGER NOT NULL,
        cds_strand CHAR(1) NOT NULL,
        location_long VARCHAR(200),
        nr_protein_id CHAR(15),
        gene_family_id CHAR(13)
);

-- CREATE INDEX cds_genbank_cds_id_key ON coding_sequences (genbank_cds_id);
-- CREATE INDEX cds_gene_family_id_key_key ON coding_sequences (gene_family_id);
-- CREATE INDEX cds_genomic_accession_key ON coding_sequences (genomic_accession);
-- CREATE INDEX cds_genomic_accession_cds_begin_key ON coding_sequences (genomic_accession, cds_begin);
-- CREATE INDEX cds_nr_protein_id_key ON coding_sequences (nr_protein_id);

CREATE TABLE proteins (
        protein_serial_id SERIAL PRIMARY KEY,
        nr_protein_id VARCHAR(20),
        product TEXT,
        protein_family_id CHAR(13)
);

CREATE TABLE nr_protein_families (
        protein_family_id CHAR(13) PRIMARY KEY,
        is_singleton BOOL NOT NULL DEFAULT 0
);

CREATE TABLE gene_families (
        gene_family_id CHAR(13) PRIMARY KEY,
        is_orfan BOOL NOT NULL DEFAULT 0,
        protein_family_id CHAR(13) NOT NULL
);

CREATE TABLE criteria_collapse_gene_tree_clades (
		criterion_id INT PRIMARY KEY,
		criterion_name	VARCHAR(50) NOT NULL,
		criterion_definition TEXT,
		collapsed_clade_collection_creation DATE
);


CREATE TABLE collapsed_gene_tree_clades (
		gene_family_id CHAR(13) NOT NULL,
		col_clade VARCHAR(8) NOT NULL,
		cds_code VARCHAR(20) NOT NULL,
		criterion_id INT NOT NULL,
		replacement_label VARCHAR(60) DEFAULT NULL
);

CREATE INDEX col_clades_gene_family_id_key ON collapsed_gene_tree_clades (gene_family_id);
CREATE INDEX col_clades_cds_code_key ON collapsed_gene_tree_clades (cds_code);

CREATE TABLE uniptrotcode2taxid (
		code varchar(5) UNIQUE NOT NULL,
		taxid integer UNIQUE NOT NULL
);

CREATE INDEX uniprot_code_key ON uniptrotcode2taxid (code);
CREATE INDEX uniprot_taxid ON uniptrotcode2taxid (taxid);



