# Sialylation project of microbiota

![Bacteria%20sia](https://github.com/edtankian97/microbiota_sialylation/blob/teste/Bacteria%20sia.gif)

This project is part of my Ph.D. thesis, which investigates the incorporation of sialic acid by the intestinal microbiota into their cell walls. The main objective of the bioinformatic analysis is to identify, within bacterial proteomes available in the NCBI database, proteins that may be involved in the process of sialic acid incorporation. Here, proteome refers to the translated amino acid sequences derived from bacterial genomes. For this purpose, the bioinformatic analysis was divided into the following steps:

**1. Genome processing:** Retrieval of genomic data from NCBI, removal of incomplete genomes, and extraction of genome quality metrics using CheckM.

**2. HMMER models:** Construction of protein datasets and development of hidden Markov models (HMMs).

**3. Protein analysis:** Identification of sialylation pathways in filtered genomes/proteomes and in control proteomes (with known absence of sialylation pathways).

**4. Downstream analysis:** Comparative genomic analysis of bacterial genomes exhibiting sialylation pathways, with generation of figures and visualizations.

**5. Metagenomic analysis:** Detection of sialylation pathways in metagenomic datasets from diseases associated with sialylation-related processes.

## Structure of folders before everything (still in construction, ignore for while)

```
.
└── genomes_download
    ├── Protein_database
    │   ├── CMP_synthase_mixed_database.fasta
    │   ├── CMP_synthase_review_database.fasta
    │   ├── CMP_synthase_unreview_database.fasta
    │   ├── KpsM_mixed.fasta
    │   ├── KpsT_mixed.fasta
    │   ├── oacetil_plus_poli_mixed_database.fasta
    │   ├── polisialil_database.fasta
    │   ├── sialiltransferase_mixed_database_old_gold.fasta
    │   ├── sialiltransferase_review_database.fasta
    │   └── sialiltransferase_unreview_database.fasta
    ├── control_proteins
    │   ├── Bac_fragilis_ATCC_faa.faa
    │   ├── E_coli_dh10B.faa
    │   ├── F_nucleatum_faa.faa
    │   ├── M4.faa
    │   ├── MAPA1.faa
    │   ├── P_putida_kt2440.faa
    │   ├── campylobacter_jejuni_ATCC_faa.faa
    │   ├── control_proteomes.txt
    │   └── files.txt
    ├── plots_data
    │   └── itol
    ├── proteins
    └── scripts
        ├── CD_HIT_script.sh
        ├── CMP_hmm.sh
        ├── KpsM_hmm.sh
        ├── KpsT_hmm.sh
        ├── RfaH_hmm.sh
        ├── Sialiltrans_hmm.sh
        ├── hmm_models.sh
        ├── jupyter_scripts
        │   ├── Checkm_refseq_Reanalise_V2.ipynb
        │   ├── hmm_process.ipynb
        │   ├── host_distribution.ipynb
        │   ├── itol_notation.ipynb
        │   ├── pie_data.ipynb
        │   └── retrieve_assembly_info.ipynb
        ├── mafft_align.sh
        ├── o_acetiltrans_poli_hmm.sh
        ├── phylo.sh
        ├── polisialiltrans_hmm.sh
        ├── rename_control_files.sh
        ├── rename_fasta.sh
        ├── rename_fasta_control.sh
        ├── rename_file.sh
        └── teste_hmm_control.sh
```
## Recommendations

To download this repository, run:

```
git clone https://github.com/edtankian97/microbiota_sialylation.git
```
Observation: **assembly_summary.txt** may be corrupted when downloaded with the entire github repository. An easy fix is to remove the file and then download it directly from Github website.

## 1. Genome processing

### 1.1 Retrieve genome information

First of all, genomes were downloaded with the wget command from [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt). The original dataset was named **assembly_summary.txt** and it can be found in the **./genomes_download** directory.

### 1.2 Filtering the retrieved dataset from NCBI
```
mv assembly_summary.txt ./genomes_download && mv CheckM_report_prokaryotes.txt ./genomes_download
cd genomes_download
grep -c “Complete” assembly_summary.txt
grep "Complete" assembly_summary.txt > assembly_complete
cut -f1,8,9,20 assembly_complete > assembly_complete_summary.tsv #retrieve info that I want to get
```

### 1.3 Retrieve CheckM information
Please refer to script: **“01.Checkm_refseq_Reanalise_V2_R.ipynb”**, allocated at /genomes_download/scripts/jupyter_scripts/. 
Run Part 01, then go back to Readme and run the following code chunks. 

If you haven’t installed miniconda or anaconda yet, please follow the instructions in this [link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).

**-Creation of conda environment and installation of ncbi_datasets. More information about ncbi_datasets, click this [link](https://github.com/ncbi/datasets)**
```
conda create -n ncbi_datasets python=3.8 #creation of the anaconda environment: Digit y or yes to continue the installation. If error occurs, you might update the python version.
conda activate ncbi_datasets #Activation of the environment. Do this after creation of the environment
conda install -c conda-forge ncbi-datasets-cli
```
**-retrieve missing data for completeness from ncbi_datasets**
```
sed -i '1d' GCF_complete_without_checkM.txt
xargs -a GCF_complete_without_checkM.txt -I {} datasets summary genome accession {} --as-json-lines | dataformat tsv genome --fields organism-name,accession,checkm-completeness,checkm-contamination > remain_CheckM_data.tsv
mv remain_CheckM_data.tsv remain_CheckM_data_complete.tsv
```
Please refer to script: **“01.Checkm_refseq_Reanalise_V2_R.ipynb”**, allocated at /genomes_download/scripts/jupyter_scripts/. 
Run Part 02, then go back to Readme and run the following code chunks. 

Run CheckM to get completeness and contamination of missing data.  Information on instalation can be found [here](https://github.com/Ecogenomics/CheckM/wiki/Installation). Below, you can find the commands for installation
```
conda create -n checkm python=3.9
conda activate checkm
conda install -c bioconda numpy matplotlib pysam
conda install -c bioconda hmmer prodigal pplacer
pip3 install checkm-genome
```
Download CheckM database and substitute /path/to/my_checkm_data to where the database is located. Then, execute CheckM script
```
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
export CHECKM_DATA_PATH=/path/to/my_checkm_data
```
**Download missing genomes for contamination analysis and run CheckM**
```
pwd #you should be in the directory genomes_download
sed '1d' remain_CheckM_data_complete_with_NA.tsv > remain_CheckM_data_complete_with_NA_ID.tsv
datasets download genome accession --inputfile remain_CheckM_data_complete_with_NA_ID.tsv  #Download genomes files ".fna"
unzip ncbi_dataset.zip && mv ncbi_dataset/ remain_CheckM/
rm ncbi_dataset.zip
find ./remain_CheckM/data/GCF_000*/ -type f -iname "*.fna" -exec mv -v "{}" ./remain_CheckM/ \;
 bash ./scripts/checkM_ncbi.sh
```
Output folder is **checkm_result_ncbi** and output file **bin_stats_ext.tsv** will be located at **checkm_result_ncbi/storage** folder. 
But **bin_stats_ext.tsv** is also located at **genomes_download/plots_data** folder.
```
cd checkm_result_ncbi
awk -F',' '/^GCF_/ { print $1, $11, $12 }' bin_stats_ext.tsv > checkm_GCF_delim.txt #delimiter rows with comma
awk -F' ' 'BEGIN{OFS="\t"} /^GCF_/ { print $1, $6, $8 }' checkm_GCF_delim.txt > quality_report.tsv #choose right colummns that I want
```
Please refer to script: **“01.Checkm_refseq_Reanalise_V2_R.ipynb”**, allocated at /genomes_download/scripts/jupyter_scripts/. 
Run Part 03. 

## 2 HMMER models

The complete workflow for sequence retrieval and duplicate removal is described in detail at [Thais_github](https://github.com/ThaisAFM/sialic_acid_catalog)
Briefly, protein sequences were downloaded and subsequently filtered to remove redundant entries.

Duplicate sequences were removed using CD-HIT, which can be obtained following the installation instructions available at [guide](https://github.com/weizhongli/cdhit). Sequence clustering was performed with a sequence identity threshold of 100% to retain only unique protein sequences.
(Example usage of cd-hit: cd-hit -i [PROTEIN_FASTA_NAME] -o [CD_HIT_ENZYME_NAME_MODE_TYPE_OUTPUT_FILE]  -c 1.00 -n 5 ).

```
#go to where tar file of sialylation sequences fasta is located 
cd genomes_download/Protein_database

#untar the file, then delete tar file after
tar -xf fastas_sialylation_final.tar.gz --strip-components=1 
rm fastas_sialylation_final.tar.gz
tar -xf fastas_others_final.tar.gz --strip-components=1 
rm fastas_others_final.tar.gz
```
After this, it's time to do run the protein alignment. For this purpose, follow [mafft](https://mafft.cbrc.jp/alignment/software/) for installation instructions.  
(Example of mafft usage: mafft --auto [CD_HIT_ENZYME_NAME_MODE_TYPE_OUTPUT_FILE] > [ENZYME_NAME_MODE_TYPE_OUTPUT_FILE]_mafft.fasta).
Results will be located at genomes_download/Protein_database/mafft_align/ director.
```
bash ../scripts/mafft_align.sh 
cd genomes_download/Protein_database/mafft_align/
ls
```
In the end, let's construct protein models with [HMMER](https://github.com/EddyRivasLab/hmmer).In this link aside, I've installed the 3.4 version for creation of models and sialylation protein's annotation.
```
bash ../../scripts/hmm_models.sh #results will be located in this PATH: genomes_download/Protein_database/mafft_align/hmm_models 
```
Let's move hmm models for external rings annotation into another directory
```
cd ../../../ #you must be located in Protein_database
mkdir external_rings_models
mv ./mafft_align/hmm_models/neu*hmm ./external_rings_models
mv ./mafft_align/hmm_models/kps*hmm ./external_rings_models
ls ./external_rings_models #check if files are there
cd ../ #you must be located in genomes_download

```
## 3. Protein analysis ### StOPPED HERE

### 3.1 Protein analysis: Control proteomes

First, we will download some proteomes based on NCBI ID which are presented in **control_proteomes.txt**. Other ones are already available in **control_proteins** folder. Those are from ATCC web and can only be download with authorizated access. Please, just use
these fasta files for academic purposes.
```
conda activate ncbi_datasets
cd ./control_proteins/ 
datasets download genome accession --inputfile control_proteomes.txt --include protein --filename control.zip
unzip control.zip -d control
ls control
```
After download, process the files
```
bash ../scripts/rename_control_files.sh #rename the files based on their directories
find ./control/ncbi_dataset/data/GCF*/ -type f -iname "*.faa" -exec mv -v "{}" ./ \; #move files
ls #see moved files 
while read line; do eval mv $line; done < files.txt #rename with species names
bash ../scripts/rename_fasta_control.sh #rename fasta header with filename
less GCF_004015025.1_Akker_munciph_NEG.faa #see content of a file
```
Create directories to organize the results.
```
mkdir HMMER_CONTROL_RESULTS && cd HMMER_CONTROL_RESULTS
bash ../../scripts/teste_hmm_control.sh
```
Join all output files for each enzyme
```
cat neuA*_output.tsv > all_CMP_neuA_control_output.tsv
cat lic3X*_output.tsv > all_lic3X_sialil_control_output.tsv
cat lst*_output.tsv > all_lst_sialil_control_output.tsv
cat pm0188*_output.tsv > all_pm0188_sialil_control_output.tsv
cat PF06002*_output.tsv > all_PF06002_sialil_control_output.tsv
cat PF11477*_output.tsv > all_PF11477_sialil_control_output.tsv
cat IPR010866*_output.tsv > all_IPR010866_polisialil_control_output.tsv
cat neuS*_output.tsv > all_neuS_polisialil_control_output.tsv
```
Do the same thing for coverage files
```
cat neuA*_output.tsv_coverage > all_CMP_neuA_control_output_coverage.tsv
cat lic3X*_output.tsv_coverage > all_lic3X_sialil_control_output_coverage.tsv
cat lst*_output.tsv_coverage > all_lst_sialil_control_output_coverage.tsv
cat pm0188*_output.tsv_coverage > all_pm0188_sialil_control_output_coverage.tsv
cat PF06002*_output.tsv_coverage > all_PF06002_sialil_control_output_coverage.tsv
cat PF11477*_output.tsv_coverage > all_PF11477_sialil_control_output_coverage.tsv
cat IPR010866*_output.tsv_coverage > all_IPR010866_polisialil_control_output_coverage.tsv
cat neuS*_output.tsv_coverage > all_neuS_polisialil_control_output_coverage.tsv
```
Process output file
```
sed -i '/#/d' *_output.tsv
sed -i 's/ \{1,\}/\t/g' *_output.tsv 
```

### 3.2 Protein analysis: NCBI analysis
After filtration with CheckM, the proteomes were downloaded from NCBI
```
cd ../../ #go to genomes_download folder and then move checkm_filter_v2_complete.tsv to this folder
cut -f2 checkm_filter_v2_complete.tsv > checkm_filter_v2_complete_ID.tsv
sed -i '1d' checkm_filter_v2_complete_ID.tsv
datasets download genome accession --inputfile checkm_filter_v2_complete_ID.tsv --include protein --dehydrated  --filename proteins.zip
rm -rfv proteins/
unzip proteins.zip -d proteins
datasets rehydrate --directory proteins
```
After the download, proteins files were renamed with their own directory name 
```
bash ../scripts/rename_files.sh
```
With the proteins with their respective names, you can move them to the **proteins** directory
```
find proteins/ncbi_dataset/data/GCF*/ -type f -iname "*.faa" -exec mv -v "{}" ./proteins/ \;
mkdir proteins/proteins_sialylation/final_complete_sialylation/ #create more directories
```
Proteomes are now in **proteins** directory and then you can edit their fasta header, so we can identify them later on HMMER analysis.
For this, do the following command:
```
cd proteins
for f in *.faa; do sed -i "s/^>/>${f}_/" "$f"; done
```
Now the proteins files are ready to be analised. Now, let's organize subdiretories to store the results for each enzyme of sialylation process
```
cd .. #must be at genomes_download folder
mkdir HMMER_analysis # you must be located at **genomes_download** folder
cd HMMER_analysis
mkdir neuA_out pm0188_out PF11477_out neuS_out PF06002_out lst_out lic3X_out IPR010866_out
cd ../../ #must be at genomes_download folder
```
Scripts for each enzyme model will be available with their respective names in **scripts** directory

```
cd scripts/
bash ncbi_teste_own_hmmer.sh
```

With all done, now it's time to start the process of coverage, e-value and bit-score filter. First, let's cat coverage results for each protein
```
cd ../HMMER_analysis/
find ./ -type f -name 'neuA*coverage' -exec cat {} + > CMP_coverage.tsv
find ./ -type f -name 'neuS*coverage' -exec cat {} + > neuS_coverage.tsv
find ./ -type f -name 'lic3X*coverage' -exec cat {} + > lic3X_coverage.tsv
find ./ -type f -name 'lst*coverage' -exec cat {} + > lst_coverage.tsv
find ./ -type f -name 'pm0188*coverage' -exec cat {} + > pm0188_coverage.tsv
find ./ -type f -name 'PF11477*coverage' -exec cat {} + > PF11477_coverage.tsv
find ./ -type f -name 'PF06002*coverage' -exec cat {} + > PF06002_coverage.tsv
find ./ -type f -name 'IPR010866*coverage' -exec cat {} + > IPR010866_coverage.tsv
```
Now, go to the script **cover_hmm_filter_NCBI.ipynb** which is loccated in the path: microbial_sialylation/genomes_download/scripts/jupyter_scripts/ and follow the script. We will filter first by coverage value. Results of coverage are alrealdy available at **plots_data** folder.

After part 01 with coverage assessment, follow this for each core enzyme. First we are going to remove the header of each file, then get filenames to move into their respective directories.  The, we will cat output files and processed them.

**NeuA**
```
cd ./HMMER_analysis
sed -i '1d' CMP_NCBI_ID_cover_filter_40.tsv
cut -f2 CMP_NCBI_ID_cover_filter_40.tsv > CMP_NCBI_ID_cover_filter_40_modified.tsv
for file in $(cat ./CMP_NCBI_ID_cover_filter_40_modified.tsv); do mv "$file" ./neuA_out/; done
cd neuA_out
# join all output tsv files into one
find ./ -type f -name '*protein_output*' -exec cat {} + > new_file.tsv
#process output file. Remove each line with "#" and delim the file with tab
sed -i '/#/d' new_file.tsv
sed  's/ \{1,\}/\t/g' new_file.tsv > file_output_CMP.tsv
cd ../ #must be located at HMMER_analysis folder
```
**lic3X**
```
sed -i '1d' Lic3_NCBI_ID_cover_filter_40.tsv
cut -f2 Lic3_NCBI_ID_cover_filter_40.tsv > Lic3_NCBI_ID_cover_filter_40_modified.tsv
for file in $(cat ./Lic3_NCBI_ID_cover_filter_40_modified.tsv); do mv "$file" ./lic3X_out/; done
cd lic3X_out
# join all output tsv files into one
find ./ -type f -name '*protein_output*' -exec cat {} + > new_file.tsv
#process output file 
sed -i '/#/d' new_file.tsv
sed  's/ \{1,\}/\t/g' new_file.tsv > lic3X_cover_filter_all.tsv
```
**lst**
```
cd ../ #must be located at HMMER_analysis folder
sed -i '1d' LST_NCBI_ID_cover_filter_40.tsv
cut -f2 LST_NCBI_ID_cover_filter_40.tsv > LST_NCBI_ID_cover_filter_40_modified.tsv
for file in $(cat ./LST_NCBI_ID_cover_filter_40_modified.tsv); do mv "$file" ./lst_out/; done
cd lst_out
# join all output tsv files into one
find . -type f -name '*protein_output*' -exec cat {} + > new_file.tsv
#process output file. Remove each line with "#" and delim the file with tab
sed -i '/#/d' new_file.tsv
sed  's/ \{1,\}/\t/g' new_file.tsv > lst_cover_filter_all.tsv
```
**neuS**
```
cd ../ #must be located at HMMER_analysis folder
sed -i '1d' neuS_thais_NCBI_ID_cover_filter_40.tsv
cut -f2 neuS_thais_NCBI_ID_cover_filter_40.tsv > neuS_thais_NCBI_ID_cover_filter_40_modified.tsv
for file in $(cat ./neuS_thais_NCBI_ID_cover_filter_40_modified.tsv); do mv "$file" ./neuS_out/; done
cd neuS_out
# join all output tsv files into one
find . -type f -name '*protein_output' -exec cat {} + > new_file.tsv
#process output file. Remove each line with "#" and delim the file with tab
sed -i '/#/d' new_file.tsv
sed  's/ \{1,\}/\t/g' new_file.tsv > neuS_cover_filter_all.tsv
```
**pm0188**
```
cd ../ #must be located at HMMER_analysis folder
sed -i '1d' pm1088_NCBI_ID_cover_filter_40.tsv
cut -f2 pm1088_NCBI_ID_cover_filter_40.tsv > pm1088_NCBI_ID_cover_filter_40_modified.tsv
for file in $(cat ./pm1088_NCBI_ID_cover_filter_40_modified.tsv); do mv "$file" ./pm0188_out/; done
cd pm0188_out
# join all output tsv files into one
find . -type f -name '*protein_output*' -exec cat {} + > new_file.tsv
#process output file. Remove each line with "#" and delim the file with tab
sed -i '/#/d' new_file.tsv
sed  's/ \{1,\}/\t/g' new_file.tsv > pm0188_cover_filter_all.tsv
```
**PF06002**
```
cd ../ #must be located at HMMER_analysis folder
sed -i '1d' PF06002_NCBI_ID_cover_filter_40.tsv
cut -f2 PF06002_NCBI_ID_cover_filter_40.tsv > PF06002_NCBI_ID_cover_filter_40_modified.tsv
for file in $(cat ./PF06002_NCBI_ID_cover_filter_40_modified.tsv); do mv "$file" ./PF06002_out/; done
cd PF06002_out
# join all output tsv files into one
find . -type f -name '*protein_output*' -exec cat {} + > new_file.tsv
#process output file. Remove each line with "#" and delim the file with tab
sed -i '/#/d' new_file.tsv
sed  's/ \{1,\}/\t/g' new_file.tsv > PF06002_cover_filter_all.tsv
```

**PF11477**
```
cd ../ #must be located at HMMER_analysis folder
sed -i '1d' PF11477_NCBI_ID_cover_filter_40.tsv
cut -f2 PF11477_NCBI_ID_cover_filter_40.tsv > PF11477_NCBI_ID_cover_filter_40_modified.tsv
for file in $(cat ./PF11477_NCBI_ID_cover_filter_40_modified.tsv); do mv "$file" ./PF11477_out/; done
cd PF11477_out
# join all output tsv files into one
find . -type f -name '*protein_output*' -exec cat {} + > new_file.tsv
#process output file. Remove each line with "#" and delim the file with tab
sed -i '/#/d' new_file.tsv
sed  's/ \{1,\}/\t/g' new_file.tsv > PF11477_cover_filter_all.tsv
```
**IPR010866**
```
cd ../ #must be located at HMMER_analysis folder
sed -i '1d' IPR01086_NCBI_ID_cover_filter_40.tsv 
cut -f2 IPR01086_NCBI_ID_cover_filter_40.tsv > IPR01086_NCBI_ID_cover_filter_40_modified.tsv
for file in $(cat ./IPR01086_NCBI_ID_cover_filter_40_modified.tsv); do mv "$file" ./IPR010866_out/; done
cd IPR010866_out
# join all output tsv files into one
find . -type f -name '*protein_output*' -exec cat {} + > new_file.tsv
#process output file. Remove each line with "#" and delim the file with tab
sed -i '/#/d' new_file.tsv
sed  's/ \{1,\}/\t/g' new_file.tsv > IPR010866_cover_filter_all.tsv
cd ../ #must be located at HMMER_analysis folder
```

Now, go to the script **hmm_process.ipynb** which is loccated in the path: microbial_sialylation/genomes_download/scripts/jupyter_scripts/ and follow the script to process output file. We will filter by bit-score and e-value.

With the result, let's move faa files that passed the final filter.
```
cd ../proteins/
sed -i '1d' completesialylation_GCF_ID.tsv #remove header
```
Move now the files into **proteins_sialylation** folder
```
 while read id; do
  mv "${id}"* proteins_sialylation/ 
done < completesialylation_GCF_ID.tsv
```
Check files and how many passed
```
ls ./proteins_sialylation
ls proteins_sialylation | wc -l 
```

### 3.3 Protein analysis: Interproscan analysis
First, download Interproscan tar file. For this purpose, I followed the manual by this [link](https://interproscan-docs.readthedocs.io/en/v5/UserDocs.html#obtaining-a-copy-of-interproscan).

```
cd ../../ #must be at genome_download folder
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.76-107.0/interproscan-5.76-107.0-64-bit.tar.gz
tar -pxvzf interproscan-5.76-107.0-*-bit.tar.gz
conda install bioconda::seqkit #download seqkit, which will retrieve fasta sequences for interproscan analysis
```
Tsv files with sequences ID for fasta sequences to retrieve will be at **/proteins/proteins_sialylation/** folder after hmm_process R's scripts had been finished. Before take the sequences, you must rename fasta header of proteins files to retrieve the sequences successfully. 

```
cd ./proteins/proteins_sialylation/
sed -i 's/ .*//' *.faa
cd ../../../scripts # must be at genomes_download folder
```
So now, execute the script below to retrieve sequences
```
bash retrieve_sequences_for_interpro.sh
```
Sequences retrieved will have **_retrieved_now** tag and will be at **scripts** folder. Let's move to another place
```
mkdir ../Interpro_analysis/Interpro_results
mv *_retrieved_now ../Interpro_analysis
conda deactivate
bash interpro_analysis.sh
conda activate ncbi_datasets
cd .. #must be at genomes_download folder
```
Results will be at this path **Interpro_analysis/Interpro_results**. Final results will be available at **plots_data/Interpro_results/** folder, which will be useful to execute Interproscan R'script
```
ls ./Interpro_analysis/Interpro_results
```
Now it's time to execute **Interpro_results** R'script to filter sequences based on signatures. This script is available at **scripts/jupyter_scripts/** folder


Final result with all proteomes that passed interproscan are available in the file **complete_sialylation_interpro_filtration_final.tsv** which can be encounter at **genomes_download/plots_data/** folder. This is going to be used to extract info for the next topic **Downstream analysis**

# 4. Downstream analysis

## 4.1 Datasets for plots
This topic and subtopics forwards are about how to get data that will be important to create the plots.

### 4.1.1 Information of genomes with sialylation pathway

**After Interproscan analysis of core enzymes, do the following to get information of genomes that have sialylation pathway**
```
#get file with protein ID with whole sialylation pathway
cd plots_data/
less complete_sialylation_interpro_filtration_final.tsv #see the data
#process and download zip file to extract information
sed -i '1d' complete_sialylation_interpro_filtration_final.tsv
datasets download genome accession --inputfile complete_sialylation_interpro_filtration_final.tsv 

#select desired fields
dataformat tsv genome --package ncbi_dataset.zip --fields accession,assminfo-biosample-geo-loc-name,assminfo-biosample-host,assminfo-biosample-host-disease,assminfo-biosample-source-type,assmstats-gc-percent,assmstats-total-sequence-len,organelle-assembly-name,organism-name,organism-tax-id > accession_complete_fields.tsv
```
Final file **accession_complete_fields.tsv** is already at **genomes_download/plots_data/** folder

### 4.1.2 Taxonomy information

Final file **accession_complete_fields.tsv** is already at **genomes_download/plots_data/** folder

### 4.1.3 Phylogenetic tree

```
cp complete_sialylation_interpro_filtration_final.tsv ../proteins/proteins_sialylation/
cd ../proteins/proteins_sialylation/
for file in $(cat ./complete_sialylation_interpro_filtration_final.tsv); do cp "$file" ./final_complete_sialylation/; done
cd ./final_complete_sialylation/
ls ./final_complete_sialylation/*faa > representative_species.txt
sed 's/_protein.faa//g' representative_species.txt > representative_species_modified.txt
cp representative_species_modified.txt ../../../genomes_download/plots_data/
```

To generate a tree from phylophlan, first you must download it. Check this [link](https://github.com/biobakery/phylophlan) with the procedures.
```
conda create -n "phylophlan" -c bioconda phylophlan=3.1.1
conda activate phylophlan
```
**phylophlan database setup and installation**
I followed instructions upon this [link](https://github.com/biobakery/phylophlan/wiki#databases). I followed the option 2 and installed the **phylophlan** database. Download **phylophlan_databases.txt** and follow the instructions below.
```
cd .. #must be in protein folder
mkdir protein_tree && cd protein_tree
mkdir phylophlan_database && cd phylophlan_database
cat phylophlan_databases.txt # copy and paste one of the links to **phylophlan** database (not amphora)
tar -xf phylophlan.tar
bunzip2 -k phylophlan/phylophlan.bz2
cd ..
```
**phylophlan configuration file**
```
phylophlan_write_config_file.py --db_aa diamond --map_aa diamond --msa mafft \
--trim trimal --tree1 fasttree --tree2 raxml -o genome_ed_config.cfg \
--db_type a
```
**Run script of phylophlan analysis**
```
cd ../../../scripts/
bash phylo.sh #make sure phylophlan's conda environment is activated
```
phylophlan generates a lot of files, but the most important is refine tree called **RAxML_result.proteins_unique_comm_sia_refined.tre** which was used for tree annotation with iTOL.
This output is present in **genomes_download** folder

## 4.2 Genome information

Follow the script **retrieve_genome_info.ipynb** which is loccated in the path: microbial_sialylation/genomes_download/scripts/jupyter_scripts/

## 4.3 Host distribution

Follow the script **host_distribution.ipynb** which is loccated in the path: microbial_sialylation/genomes_download/scripts/jupyter_scripts/

## 4.4 Species distribution

Follow the script **species_distribution.ipynb** which is loccated in the path: microbial_sialylation/genomes_download/scripts/jupyter_scripts/

## 4.5 iTOL annotation

### 4.5.1 Hmmer analysis of external rings

```
cd ../ #If you are one level above at jupyter_scripts folder, do this command to go to scripts folder
bash external_rings_hmmer.sh #do the annotation of proteins for external rings
```

Results for each protein: KpsM, KpsT, neuO and neuD will be located at **../Protein_database/external_rings_models/external_rings_output/**

```
cd ../Protein_database/external_rings_models/external_rings_output/
mkdir kpsM_out kpsT_out kpsD_out neuD_out neuO_out
find ./ -type f -name 'kpsM*coverage' -exec cat {} + > all_kpsM_coverage.tsv
find ./ -type f -name 'kpsT*coverage' -exec cat {} + > all_kpsT_coverage.tsv
find ./ -type f -name 'kpsD*coverage' -exec cat {} + > all_KpsD_coverage.tsv
find ./ -type f -name 'neuO*coverage' -exec cat {} + > all_neuO_coverage.tsv
find ./ -type f -name 'neuD*coverage' -exec cat {} + > all_neuD_coverage.tsv
```
Final resulted file are located at **microbiota_sialylation/genomes_download/plots_data/hmmer_out/**
First process resulted of coverage file in **hmm_process_external_rings.ipynb** script **Part 01** and then continue

**KpsM**
```
sed -i '1d' cover_kpsM_ncbi_filter.tsv
cut -f2 cover_kpsM_ncbi_filter.tsv > kpsM_cover_filter_40_ID.tsv 
#transfer files to another folder
for file in $(cat ./kpsM_cover_filter_40_ID.tsv); do mv "$file" ./kpsM_out/; done
cd kpsM_out
find . -type f -name '*protein_output*' -exec cat {} + > file_output_kpsM.tsv
#process output file
sed -i '/#/d' file_output_kpsM.tsv
sed -i 's/ \{1,\}/\t/g' file_output_kpsM.tsv
cd ..
```

**KpsT**
```
#remove first line and first columm to get only ID
sed -i '1d' cover_kpsT_ncbi_filter.tsv
cut -f2 cover_kpsT_ncbi_filter.tsv > kpsT_cover_filter_40_ID.tsv 
#transfer files to another folder
for file in $(cat ./kpsT_cover_filter_40_ID.tsv); do mv "$file" ./kpsT_out/; done
cd kpsT_out
find . -type f -name '*protein_output*' -exec cat {} + > file_output_kpsT.tsv
#process output file
sed -i '/#/d' file_output_kpsT.tsv
sed -i 's/ \{1,\}/\t/g' file_output_kpsT.tsv
cd ..
```
**KpsD**
```
#remove first line and first columm to get only ID
sed -i '1d' cover_kpsD_ncbi_filter.tsv
cut -f2 cover_kpsD_ncbi_filter.tsv > kpsD_cover_filter_40_ID.tsv 
#transfer files to another folder
for file in $(cat ./kpsD_cover_filter_40_ID.tsv); do mv "$file" ./kpsD_out/; done
cd kpsT_out
find . -type f -name '*protein_output*' -exec cat {} + > file_output_kpsD.tsv
#process output file
sed -i '/#/d' file_output_kpsD.tsv
sed -i 's/ \{1,\}/\t/g' file_output_kpsD.tsv
cd ..
```
**NeuD**
```
#remove first line and first columm to get only ID
sed -i '1d' cover_neuD_ncbi_filter.tsv
cut -f2 cover_neuD_ncbi_filter.tsv > cover_neuD_ncbi_filter_modified.tsv
#transfer files to another folder
for file in $(cat ./cover_neuD_ncbi_filter_modified.tsv); do mv "$file" ./neuD_out/; done
cd neuD_out
find . -type f -name '*protein_output*' -exec cat {} + > file_output_neuD.tsv
#process output file
sed -i '/#/d' file_output_neuD.tsv
sed -i 's/ \{1,\}/\t/g' file_output_neuD.tsv
cd ..
```

**NeuO**
```
#remove first line and first columm to get only ID
sed -i '1d' cover_neuO_ncbi_filter.tsv
cut -f2 cover_neuO_ncbi_filter.tsv > neuO_cover_filter_40_ID.tsv 
#transfer files to another folder
for file in $(cat ./neuO_cover_filter_40_ID.tsv); do mv "$file" ./neuO_out/; done
cd neuO_out
find . -type f -name '*protein_output*' -exec cat {} + > file_output_neuO.tsv
#process output file
sed -i '/#/d' file_output_neuO.tsv
sed -i 's/ \{1,\}/\t/g' file_output_neuO.tsv
cd ..
```
Return to the **hmm_process_external_rings.ipynb** script again for **Part 02**. Final resulted files are already located at **microbiota_sialylation/genomes_download/plots_data/hmmer_out/**

### 4.5.2 Interpro_analysis

Retrieve sequences for analysis.
```
conda deactivate
cd ../../../../scripts #must be at scripts folder
bash retrieve_sequences_for_interpro_external_rings.sh 
ls ../Interpro_analysis # files will be at Interpro_analysis folder
```
Run Interpro analysis
```
bash interpro_analysis_external_rings.sh
cd .. # must be at genomes_download
```
Check output files at **/Interpro_analysis/Interpro_results/** PATH. Our results are already available at **plots_data/Interpro_results/** PATH
```
ls ./Interpro_analysis/Interpro_results/
```
Now, it's time to process the results with interproscan'R script for external rings: **Interpro_analysis_external_rings_final_script**. Final output files wil be used for iTOL annotation of external rings.


### 4.5.3 Virulence factors

Virulence factors were detected by abricate software with VFDB database.
```
#install abricate
conda activate ncbi_datasets
conda install -c conda-forge -c bioconda -c defaults abricate
```

After install abricate, download genomes.
```
cd ../../../proteins/proteins_sialylation/final_complete_sialylation/
datasets download genome accession --inputfile ./complete_sialylation_interpro_filtration_final.tsv --include genome --dehydrated --filename genomes_unique.zip
unzip genomes_unique.zip -d genomes_unique
datasets rehydrate --directory genomes_unique
find ./genomes_unique/ncbi_dataset/data/GCF*/ -type f -iname "*.fna" -exec mv -v "{}" ./genomes_unique \; #move files
ls ./genomes_unique # list files
```
Start the process of analysis
```
abricate ./genomes_unique/*fna --db vfdb --csv --minid 70 --mincov 60 > out_70id_60cov.csv
mv out_70id_60cov.csv ../../../../plots_data/itol/
```
### 4.5.4 Processing of resulted files for iTOL annotation

Follow the script **itol_notation.ipynb** which is loccated in the path: microbial_sialylation/genomes_download/scripts/jupyter_scripts/
After this, each dataset created in the script was concatenated with iTol dataset
```
cd ../../../../plots_data/itol/
#remove header
sed -i '1d' kpsT_represent_itol_EANS_paper.tsv
sed -i '1d' kpsM_itol_represent_EANS_paper.tsv
sed -i '1d' kpsD_itol_represent_EANS_paper.tsv
sed -i '1d' neuD_itol_represent_EANS_paper.tsv
sed -i '1d' neuO_itol_represent_EANS_paper.tsv
sed -i '1d' vfdb_itol_represent_EANS.tsv
sed -i '1d' representative_labels_itol.tsv
```
Now, join files to get a dataset for itol
```
 #Join files
cat kpst_itol_representative.txt kpsT_represent_itol_EANS_paper.tsv > final_kpsT_itol.txt
cat kpsM_itol_representative.txt kpsM_itol_represent_EANS_paper.tsv > final_kpsM_itol.txt
cat kpsD_itol_representative.txt kpsD_itol_represent_EANS_paper.tsv > final_kpsD_itol.txt
cat neuD_itol_representative.txt neuD_itol_represent_EANS_paper.tsv > final_neuD_itol.txt
cat neuO_itol_representative.txt neuO_itol_represent_EANS_paper.tsv > final_neuO_itol.txt
cat vfdb_itol_representative_EANS.txt vfdb_itol_represent_EANS.tsv > final_vfdb_itol.txt
cat ed_tree_label.txt representative_labels_itol.tsv > ed_tree_label_represent.txt
cat representative_phyli_sialylation.txt representative_phylum_itol.tsv > final_representative_phylum.txt
```
Now you can upload the files in iToL site. First submit and open final tree file which is inside **plots_data/itol/** folder
with the name **RAxML_result.final_complete_sialylation_refined.tre**

# 5. Metagenomic analysis

## 5.1 Study 01 colon cancer - France

### 5.1.1 Download of fastq files

First download  fastq files
```
conda install -c bioconda sra-tools #Install sra-tools
cd ../../../metagen_files/Study01_france_cancer/
cat *_get_ID > all_cancer_download.sh
bash cancer_01_download.sh
```
### 5.1.2 Quality filtering

After this, it's time to do a filtering with fastp
```
conda install -c bioconda fastp #download fastp
bash script_fastp_filtering.sh
ls outputs_fastp #folder with result of filtering
```
Results will be at **outputs_fastp** folder

### 5.1.3 Download of human genome
```
conda install -c bioconda bowtie2 #download bowtie2 via conda
bash script_get_human_genome_GRch38.sh
```
Output will be at **human_genome** folder.

### 5.1.4 Human genome alignment

We will align all fastq files against human genome 
```
bash script_align_reads_human_genome.sh
ls outputs_bowtie2_FR # list output files
```
Files will be saved at **outputs_bowtie2_FR** folder.

### 5.1.5 Remove human genome alignment

After alignment, we will remove all human genome alignment to our fastq files, as we want only microbiota genomes'source

```
conda install -c bioconda samtools #download samtools via conda
bash script_remove_human_genome.sh
```
Results will be saved at **clean_reads_FR** folder

### 5.1.6 Merge sequences with FLASH
```
conda install -c bioconda flash #install flash via conda
bash script_flash.sh
```
### 5.1.7 Generate contigs with MEGAHIT
```
conda install -c bioconda megahit
bash script_megahit.sh
```
### 5.1.8 Generate bins with metabat
```
conda install bioconda::metabat2
bash create_index.sh #create index
bash prepare_to_binning.sh #coverage count
bash generate_bins.sh #generate bins
bash rename_bins_files.sh #rename files based on their directories
```

### 5.1.9 Generate proteins faa files with prokka
```
cd bins_paired
mkdir all_prokka
bash prokka_script.sh
bash mv_prokka_files.sh 
bash rename_fasta_header.sh 
```
Files with faa extension are now inside **all_prokka** folder. Let's see if everything is ok. Check if files are renamed corrected with their respective directory with **less file_name.faa** command
```
cd all_prokka 
ls
```
### 5.1.10 Hmmer analysis

Now, it's time to do a hmmer annotation

```
cd ../../ #must be at Study01_france_cancer folder
bash metagen_hmmer.sh
```
Process hmmer results. Output files are in this path **metagen_files/Study01_france_cancer/output_data**
```

```
### 5.1.11 Interpro analysis

Output files with sequences fasta to run Interpro will be at this path **metagen_files/Study01_france_cancer/output_data/Interproscan_analysis**. 
```
bash ../scripts/interpro_analysis_FR.sh
```
Final results of Interproscan are available at **metagen_files/Study01_france_cancer/output_data/Interpro_results/** folder, which will be useful to execute Interproscan R'script.
Now it's time to execute **Interpro_results_FR** R'script to filter sequences based on signatures. This script is available at **metagen_files/Study01_france_cancer/** folder. Final result of  Interproscan analysis are available at  **metagen_files/Study01_france_cancer/output_data/Interpro_results/**

### 5.1.12 Phylogeny Identification of bins

Now, it's time to identify bins that passed the Interproscan analysis filtering. File with bins_ID are presented at **metagen_files/Study01_france_cancer/output_data/Interpro_results/** with the name **bins_for_identification.tsv**. Let's handle this file and extract the genomes file for identification.

```
cd ./output_data/Interpro_results/ #considering that you are at Study01_france_cancer folder
less bins_for_identification.tsv
sed -i '1d' bins_for_identification.tsv
cut -d'_' -f2-3 bins_for_identification.tsv > bins_for_identification_modified.tsv #extract right ID of bins filename
mv bins_for_identification_modified.tsv ../../bins_paired #move file to desired folder
cd ../../bins_paired #go to folder of extraction
mkdir bins_paired_sialylation #create folder to move the files

#execute command to move the right files
while IFS= read -r file; do
    find . -type f -name "$file" -exec mv -t bins_paired_sialylation {} +
done < bins_for_identification_modified.tsv

ls ./bins_paired_sialylation/
```
Download GTDB-tk for analysis. I followed GTDB-tk [manual](https://ecogenomics.github.io/GTDBTk/installing/index.html#installing)
```
cd ../../ #must be at metagen_files folder
mkdir GTDB_DATA && cd GTDB_DATA
#download GTDB splitted files
wget -r -np -nH --cut-dirs=8 \
https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/auxillary_files/gtdbtk_package/split_package/
cd split_package #go to files' folder
 cat gtdbtk_r226_data.tar.gz.part_* > gtdbtk_r226_data.tar.gz #cat all files into one
 tar xvzf gtdbtk_r226_data.tar.gz #untar file
```
Also, download GTDB-tk via conda
```
conda create -n gtdbtk-2.6.1 -c conda-forge -c bioconda gtdbtk=2.6.1
conda activate gtdbtk-2.6.0
```
Execute GTDB script.
```
bash metagen_FR_GTDB.sh
```
Summary results are at **metagen_files/Study01_france_cancer/output_data** folder with file named **gtdbtk.bac120.summary_FR.tsv**

### 5.1.13 Check quality of bins with CheckM
```
conda activate checkm
export CHECKM_DATA_PATH=YOUR/PATH/TO/CHECKM/DATA
bash checkm_script_metagen_FR.sh

```
Results will be at **checkm_result_bins_with_sia** folder
