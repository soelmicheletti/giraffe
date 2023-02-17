# Chromatin accessibility on mus musculus

Analysis of simultaneous profiling of transcriptome and chromatin accessibility within single cells based on data provided in the paper

*An ultra high-throughput method for single-cell joint analysis of open chromatin and transcriptome* by C. Zhu et al. 

The organism under study is [mus musculus](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=10090), with the reference genome GRCm38. Here a review of the data in `raw`:

- Data in `raw/Adult_CTX_DNA` and `raw/Adult_CTX_RNA` can be found on NCBI Gene Expression Omnibus (GEO) (http://www.ncbi.nlm.nih.gov/geo/) under accession number [GSE130399](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130399). Since data are sparse, they are conveniently stored in  [mtx format](https://math.nist.gov/MatrixMarket/formats.html) to save memory. 
   - `raw/Adult_CTX_DNA/matrix.mtx` contains the cell-specific number of counts in every region of the genome.  The cell identifiers are to be found in `raw/Adult_CTX_DNA/barcodes.tsv`), and the genome regions (indexed by chromosome and an offset) are in `raw/Adult_CTX_DNA/raw/Adult_CTX_DNA/genes.tsv`.
   - `raw/Adult_CTX_RNA/matrix.mtx` contains the cell-specific gene expression for each protein coding gene in mice. The genes names are saved in `raw/Adult_CTX_DNA/barcodes.tsv` (both ensembl ID and alias), while the cell identifiers (that correspond to the ones in the bullet point above) are in `raw/Adult_CTX_RNA/barcodes.tsv`. 
- Data in `raw/mart_export.txt` contains the location of each gene (including chromosome and strand). This dataset was collected from [ENSEMBL](http://useast.ensembl.org/biomart/martview/) by selecting multiple attributes (Gene stable ID, Gene stable ID version, Transcript stable ID , Transcript stable ID version, Gene start (bp), Gene end (bp), Ensembl Canonical, Gene type, Gene name, Strand, Chromosome/scaffold name) for mouse genes (GRCm38). 
- The PPI is taken from [MIPPIE](http://cbdm-01.zdv.uni-mainz.de/~galanisl/mippie/download.php), a project for mouse integrated protein-protein interaction reference. 
   - `raw/mippie_ppi_v1_0.tsv` contains 42,610 interactions between 10,886 proteins. 
   - `mippie_proteins_v1_0.tsv` contains the translation between entrez and protein number. 
- The motif is taken from [here](https://sites.google.com/a/channing.harvard.edu/kimberlyglass/tools/resources?authuser=0), where data have properly been preprocessed for version GRCm38. In particular, the mappings were produced by using FIMO to scan a genome for sets of  CIS-BP Single Species DNA-motifs curated in the MEME suite. 
    - `raw/mm10_reseq_-750_250_1e-4.txt` contains the edges in our TF-edge prior. 
    - `raw/Mus_musculus_motifinfo.txt` contains the mapping to the TF names. 

In `preprocessing.py` we prepare the data that we will use in our analysis. In particular:

- We transform the motif to an adjacency matrix, which is the format required by GIRAFFE. 
- We extract the PPI in adjacency matrix form, using a consistent naming with the motif file.
- For each gene and each cell, we save the number of counts in the promoter region. 