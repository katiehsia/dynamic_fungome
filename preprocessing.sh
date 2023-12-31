
#AMPLICON SEQUENCE VARIANT PRE-PROCESSING AND ALIGNMENT (VIA QIIME2) 

#variables
module load qiime2/2021.11
source activate qiime2-2021.11
MANIFEST=<path/to/manifest/file>
METADATA=<path/to/metadata/file>
OUTPUT=<path/to/output/folder>
TRUNC_FORWARD=<forward read cutoff>
TRUNC_REVERSE=<reverse read cutoff>
SAMPLE_DEPTH=<sample depth>
CORE=<linik/to/core/folder>

qiime tools import \--type 'SampleData[PairedEndSequencesWithQuality]' \--input-path $MANIFEST \--output-path $OUTPUT/paired-end-demux2.qza \--input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \--i-data $OUTPUT/paired-end-demux2.qza \--o-visualization $OUTPUT/paired-end-demux2.qzv

qiime dada2 denoise-paired \--i-demultiplexed-seqs $OUTPUT/paired-end-demux2.qza \--p-trunc-len-f 280 \--p-trunc-len-r 261 \--o-representative-sequences $OUTPUT/rep-seqs-dada2.qza \--o-table $OUTPUT/table-dada2.qza \--o-denoising-stats $OUTPUT/stats-dada2.qza

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-dada2.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
  
qiime feature-table filter-samples \
  --i-table table-dada2.qza \
  --m-metadata-file $METADATA \
  --o-filtered-table filtered-table.qza \
  
  qiime feature-table summarize \
  --i-table filtered-table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file $METADATA \
  
  
qiime diversity core-metrics-phylogenetic \--i-phylogeny $OUTPUT/rooted-tree.qza \--i-table $OUTPUT/filtered-table.qza \--p-sampling-depth $SAMPLE_DEPTH \--m-metadata-file $METADATA \--output-dir $BASE_PWD/$CORE \--verbose

#FUNGAL CLASSIFIER TRAINING (VIA UNITE)

qiime tools import \
 --type FeatureData[Sequence] \
 --input-path $INPUT/sh_qiime_release_10.05.2021/sh_refs_qiime_ver8_dynamic_10.05.2021.fasta \
 --output-path $OUTPUT/qiime_ver8_dynamic_10.05.2021.qza 
 
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path $INPUT/sh_qiime_release_10.05.2021/sh_taxonomy_qiime_ver8_dynamic_10.05.2021.txt \
--output-path $OUTPUT/ref-taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads $INPUT/qiime_ver8_dynamic_10.05.2021.qza  \
--i-reference-taxonomy $INPUT/ref-taxonomy.qza \
--o-classifier $OUTPUT/classifier.qza

# Creation of ASV Table

qiime feature-classifier classify-sklearn \
  --i-classifier $INPUT/classifier.qza \
  --i-reads $OUTPUT/rep-seqs-dada2.qza \
  --o-classification $OUTPUT/taxonomy.qza

qiime metadata tabulate \
  --m-input-file $OUTPUT/taxonomy.qza \
  --o-visualization $OUTPUT/taxonomy.qzv

qiime taxa barplot \
  --i-table $OUTPUT/filtered-table.qza \
  --i-taxonomy $OUTPUT/taxonomy.qza \
  --m-metadata-file $METADATA \
  --o-visualization $OUTPUT/taxa-bar-plots.qzv


qiime tools export \
--input-path $INPUT/filtered-table.qza \
--output-path $OUTPUT/filtered-table

biom convert \
  -i $OUTPUT/filtered-table/feature-table.biom \
  -o $OUTPUT/otu_table.txt \
  --to-tsv
  
qiime tools export \
--input-path $INPUT/taxonomy.qza \
--output-path $OUTPUT/taxonomy
  
qiime tools export \
--input-path $INPUT/unrooted-tree.qza \
--output-path $OUTPUT/unrooted-tree


#PROCESSING AND SAMPLE CLASSIFICATION USING RANDOM FOREST 

qiime taxa collapse \
  --i-table $INPUT/filtered-table.qza \
  --i-taxonomy $INPUT/taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table $OUTPUT/genus_table.qza
  
qiime feature-table summarize \
  --i-table genus_table.qza \
  --o-visualization genus_table.qzv \
  --m-sample-metadata-file $METADATA \

qiime feature-table filter-features \
--i-table $INPUT/genus_table.qza \
--p-min-frequency 5000 \
--p-min-samples 25 \
--o-filtered-table $OUTPUT/table-prevalent.qza \

qiime feature-table summarize \
  --i-table $INPUT/table-prevalent.qza \
  --o-visualization $OUTPUT/prevalent_genus_table.qzv \
  --m-sample-metadata-file $METADATA \

# converted to table that could be imported into R 

qiime tools export \
--input-path $INPUT/table-prevalent.qza \
--output-path $OUTPUT/table-prevalent

biom convert \
-i $OUTPUT/table-prevalent/feature-table.biom \
-o $OUTPUT/otu_table_prevalent.txt --to-tsv

# performed further filtering in R then converted back to qza file 

biom convert \
-i $INPUT/otu_table_prevalentbest.tsv \
-o $OUTPUT/prevalent2.biom --table-type "Table" --to-hdf5

qiime tools import \
--input-path $INPUT/prevalent2.biom \
--type FeatureTable[Frequency] \
--output-path $OUTPUT/prevalent2.qza

# sample classification

qiime sample-classifier classify-samples \
  --i-table $INPUT/prevalent2.qza \
  --m-metadata-file $INPUT/metadatafinalclass4.tsv \
  --m-metadata-column Remission_cat2 \
  --p-optimize-feature-selection \
  --p-step 0.025\
  --p-test-size 0.25 \
  --p-cv 10 \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 50 \
  --output-dir $OUTPUT/ML2
  
#export for ASV analysis
Qiime tools export \
--input-path rep-seqs.qza \
--output-path phyloseq
