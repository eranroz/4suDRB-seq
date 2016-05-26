This is a simple script for calling the Pol II wave frontier in 4su DRB seq.

* Input files: read coverage in matlab cell array, organized by chromosomes
* Output: csv file with transcripts and their Pol II wave frontier

# Install
```
python setup.py install
```
Or manually install [hmm_kit](https://github.com/eranroz/hmm) and clone the repository.


# Usage
1. Get 4su DRB-seq data. Refer to https://doi.org/10.1186/gb-2014-15-5-r69 for more info
2. Extract the reads to matlab array and save it in format of REPLICATE_CONDITION_TIME.mat
   Do it for all replicates and conditions and store them in same directory
3. Change config.py:
 3.1. TRANSCRIPTION_DATA_DIR - Should point to the above directory
 3.2. jump - Should be compatible to the size of the bins in the matlab array. (the size of the matlab array of chromsome should be compatible to UCSC chrom.sizes/jump
4. Run the script:
```
python hmm_polwave.py CONDITION EXPERIMENT_TIME --rep_sample REPLICATE --genome hg19
```

# References
* Gilad Fuchs, Eran Rosenthal, Debora-Rosa Bublik, Tommy Kaplan*, and Moshe Oren* [Gene body H2B monoubiquitylation regulates gene-selective RNA Polymerase II pause release and is not rate limiting for transcription elongation](http://biorxiv.org/content/early/2015/12/25/035386.full.pdf) [bioRxiv, 2015](http://biorxiv.org/content/early/2015/12/25/035386)
* Gilad Fuchs†, Yoav Voichek†, Sima Benjamin, Shlomit Gilad, Ido Amit and Moshe Oren [4sUDRB-seq: measuring genomewide transcriptional elongation rates and initiation frequencies within cells](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-5-r69), Genome biology 15.5 (2014): 1-11.‏

