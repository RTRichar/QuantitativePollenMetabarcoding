## Database Curation and Classifier Training
##### Download all rbcL, trnL, trnH and whole chloroplast sequences from NCBI Genbank
##### Download all Viridiplantae and Rhodophyta ITS2 sequences from the ITS2 Database
##### Remove next line entries from sequences.
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < GenBank_rbcL_030416.fasta > GenBank_rbcL_030416_rmNextLine.fasta 
```
##### Remove everything but the Genbank Gi number from the header of each entry
```
cat file.fasta | perl -pe 's/^ >gi\|(\d+)\|.*/$1/' > file_clean.fasta
```
##### For rbcL, trnL and trnH, run the Python_script_1.py script to remove sequences containing three or more consecutive uncalled base pairs 
```
python Python_script_1.py file_clean.fasta
```
##### For trnL and trnH, remove entries longer than 20,000 and 1,500 bp, respectively (in the below command ‘LENGTH’ should be substituted for the threshold length of the specific marker being curated)
```
awk '!/^>/ { next } { getline seq } length(seq) <= LENGTH { print $0 "\n" seq }' file_clean.fasta > file_clean_trim.fasta
```
##### Combine trnL and rbcL data with whole chloroplast sequence data so that the barcode regions of these entries can be extracted during a preliminary Metaxa2 training step
##### Get list of Gi numbers for each marker
```
grep '>' file_clean_trim.fasta | perl -pe 's/>(\d+).*/$1/g' > file.gis
```
##### Use Perl script from Sickel et al. (2016) along with the NCBI Taxonomy module to get the Linnaean lineage of each entry 
```
perl ~/PATH/TO/RDP_Akenbrand/meta-barcoding-dual-indexing/code/gi2taxonomy.pl --gis file.gis --out file.tax --species file.species.taxids --genus file.genus.taxids
```
##### For trnL and rbcL, split the entries up based on their length using the Python_script_2.py
###### For rbcL, partition sequences into 100-2500 bp and 2500+ bp files. For trnL, partition files according to the following list: 0 bp – 500 bp, 501 bp – 800 bp, 801 bp – 1600 bp, 1601 bp – 148000 bp, 148001 bp – 200000 bp. This step was necessary for us since the Metaxa2 training process relies on multiple sequence alignment. Performing the training without dividing up the sequences resulted in numerous misalignments and hindered our ability to extract the exact barcode region of interest from each entry. 
```
python Python_script_2.py file_clean_trim.fasta minimum_length max_length file_minbp_maxbp.fasta
```
##### For trnL and rbcL conduct a preliminary run of the Metaxa2 database builder tool on the sequences from each partition file
###### Doing this preliminary analysis across multiple smaller files facilitates the multiple sequence alignment and extraction of each barcode region of interest. When executing this step, representative sequences are indicated using the ‘-r’ option (representatives shown below). 
```
perl ~/PATH/TO/Metaxa2_2.2/metaxa2_dbb -o MARKER_NAME -g MARKER_NAME -t TAXONOMIES.tax -i SEQUENCES.fasta -r REPRESENTATIVE_GI_NUMBER

>83700918 (trnL representative)
TTGGATTGAGCCTTGGTATGGAAACTTACTAAGTGAAAACTTTCAAATTCAGAGAAACCCTGGAATTAAAAAGGGGCAATCCTGAGCCAAATCCTTCTTTCCGAAAACAAATAAAAGTTCAGAAAGTTAAAATCAAAAAAGGATAGGTGCAGAGACTCAAT
> 9967384 (rbcL representative)
ATGTCACCACAAACGGAGACTAAAGCAGGTGTTGGATTTAAAGCTGGTGTTAAAGATTACAGATTAACCTATTACACCCCAGATTATCAGACTAAAGACACTGATATTTTGGCAGCATTTCGGATGACGCCTCAACCAGGAGTACCCGCTGAAGAGGCAGGAGCTGCAGTAGCTGCGGAATCCTCCAGCGGTACATGGACCACTGTTTGGACCGATGGACTTACTAGTCTTGATCGTTACAAAGGACGATGCTATGATCTTGAAGCAGTTCCTGGAGAAGAAAATCAATATATTGCTTATGTTGCTTATCCATTAGATTTATTTGAAGAAGGTTCTGTTACCAATTTATTCACCTCTATTGTTGGTAATGTTTTTGGATTCAAAGCTCTACGAGCTTTACGTCTAGAAGATTTACGTATTCCTCCAGCTTATTCCAAAACTTTCCAAGGCCCACCTCATGGTATTCAAGTTGAAAGAGATAAATTAAATAAATATGGTCGTCCATTATTGGGATGTACTATTAAGCCAAAATTGGGTTTATCCGCTAAATACTATGGTAGAGCTGTATATGAATGTCTCCGTGGTGGACTTGATTTCACAAAAGATGATGAAAACGTAAATTCTCAACCTTTTATGCGTTGGAGAGATCGTTTC
```
##### Concatenate extracted sequences and taxonomies into one file for each marker and convert sequences back to upper case
```
awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' file_Prelim1.fasta > file_Prelim1_upper.fasta
```
##### Curate the taxonomies to remove undefined ranks at the leaves of the lineages as well as any extraneous artifacts such as open nomenclature or informal demarcations
###### These commands represent what was required for our databases but, in general, we recommend that you devise ways to rigorously seek out and fix these types of issues in your data because the artifacts found are highly variable from one dataset to the next.
```
cat taxonomies.tax | perl -pe 's/^(\d+).E\tRoot;/$1\t/' | perl -pe 's/s__undef__\d+;$//' | perl -pe 's/g__undef__\d+;$//' | perl -pe 's/f__undef__\d+;$//' | perl -pe 's/o__undef__\d+;$//' | perl -pe 's/([a-z]+)_\d+;/$1;/g' | sed '/k__undef/d' > taxonomies_0.tax
cat taxonomies_0.tax | perl -pe 's/s__[A-Za-z]+\ssp\..*;$//' | perl -pe 's/s__[A-Za-z]+.\saff\.\s.*;$//' | perl -pe 's/s__[A-Za-z]+\sgen\.\s.*;$//' | perl -pe 's/s__[A-Za-z]+\scf\.\s.*;$//' > taxonomies_1.tax
cat taxonomies_1.tax | perl -pe 's/s__[A-Za-z]+\sx\s[A-Za-z].*//' | perl -pe 's/s__[A-Za-z]+\s[a-z]+\sx\s[A-Za-z].*//' | perl -pe 's/s__[A-Za-z]+.*\d.*//' > taxonomies_2.tax
cat taxonomies_2.tax | perl -pe 's/s__x\s.*//' | sed '/environmental\ssample/d' | sed '/uncultured\sChlorophyta/d' | perl -pe 's/s__[A-Za-z]+\shybrid\scultivar;//' | sed '/uncultured\sChlorella/d' | perl -pe 's/s__\(.*//' > taxonomies_3.tax
```
##### Use the Python_script_3.py script to revise the taxonomic lineages of entries identified at high resolution ranks (family, genus or species) but undefined at intermediate ranks (class, order or family)
###### Performing this step provides each lineage with distinct tags so that undefined Magnoliales can be distinguished from undefined Ranuculales, for example. 
```
python Python_script_3.py taxonomies.txt
```
##### Remove duplicate sequences from each reference database
###### For this step, sequence headers must contain the taxonomic lineage of the entry, which can be done using a Perl script from Sickel et al. (2016)

```
perl ~/PATH/TO/RDP_Akenbrand/meta-barcoding-dual-indexing/code/tax2rdp_utax.pl MarkerName.tax MarkerName.fasta MarkerName

java -Xmx4g -jar /PATH/TO/RDP_Akenbrand/rdp_classifier_2.11/dist/classifier.jar rm-dupseq --infile MarkerName.rdp.fa --outfile MarkerName.rmDS.fasta --duplicates --min_seq_length 50
```
##### Perform final Metaxa2 database construction using the curated sequence and taxonomy information 
###### For rbcL and trnL, the same representative sequences shown previously were used to indicate the exact barcode region of interest. For the highly divergent trnH and ITS2 markers, the ‘-r’ option was not used and the ‘--divergent’ option was set to ‘T’.
```
perl ~/PATH/TO/Metaxa2_2.2/metaxa2_dbb -o DATABASE_NAME -g GENE_NAME -t taxonomies.tax –i input.fasta -r 83700918
```
## Classifier Evaluations
##### For evaluation of classifier performance, the methods described in Richardson et al. 2017 were used (DOI: 10.1111/1755-0998.12628; https://github.com/RTRichar/evaluating-DNA-metabarcoding)
###### For this section of the methods, Python_script_4.py is needed to crop reference sequences to 150 bp. 
```
python Python_script_4.py file.fasta 150 output.fasta
```
## Classification and Analysis of Multi-locus Pollen Metabarcoding Data
##### Use Metaxa2 to classify sequences.
```
perl ~/PATH/TO/Metaxa2_2.2/metaxa2 -1 forward_reads.fastq -2 reverse_reads.fastq -o output_name -g DATABASE_NAME -R 50
```
##### use Python_script_5 to consensus filter the classifications, returning only the families detected in at least two of the four libraries per sample, and calculate the median proportion of each consensus family across the four libraries
###### In this analysis, the median proportion is calculated from only the libraries in which the family was detected (i.e. zeros are excluded). The files must be named according to the following convention: *'sample_name.marker_name.taxonomy.txt'*. There can be no periods in the sample or marker names as the period is used to delimit *sample* from *marker*. A *‘base_names.txt’* file containing all of the sample names must be designated. Lastly, a *‘markers.txt’* containing the name of each marker used, must be designated. 
```
python Analyze_Multi-Locus.py -b base_names.txt -m marker_names.txt -pf pre-filter_threshold -mt masking_threshold -o output_name.csv -r 4
```
