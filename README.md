This repository describes steps to prepare and annotate VCF files for loading into HPDS as in https://github.com/hms-dbmi/pic-sure-hpds-copdgene.

# Preparing VCF

The input VCF file needs to be normalized for variant representaiton: to split multi-allelic variants into separate lines in VCF.
There are various ways eeto do this, but one of the way is to use bcftools (http://samtools.github.io/bcftools/).

`bcftools norm -m -any -f [fasta file for referenge genome] [input VCF file] | bgzip -c > [normalized VCF file]`


Also make the normalized VCF file tabix-indexed to enable fast-accessing during annotation steps.

`tabix -p vcf [normalized VCF file]`

(tabix is part of htslib. http://www.htslib.org/download/)


# Annotating VCF using VEP (variant effect predictor: https://www.ensembl.org/vep)

Here, we use VEP as Docker container, without direct installation to local system.
Before using VEP, we need to:
- docker image pulled (obviously)
- download approprite VEP cache files from ftp://ftp.ensembl.org/pub/current_variation/indexed_vep_cache/ for the latest version for GRCh38 reference genome. Previous versions or files for GRCh37 referenge genome can also be found in the FTP site.
  (but it needs to be the "merged" cache file which contains both RefSeq & Ensembl transcripts)
- download reference genome fasta files from ftp://ftp.ensembl.org/pub/current_fasta.

```
docker run --rm -it -v [local VEP cache directory]:/cache -v [directory for other files as needed]:/data -v [directory for input VCF file]:/work ensemblorg/ensembl-vep /opt/vep/vep \
	--cache --offline --merged \
	--species homo_sapiens \
	--compress_output bgzip \
	--input_file /work/[input VCF path/filename] \
	--output_file /work/[annotated VCF path/filename] \
	--no_stats \
	--force_overwrite \
	--assembly GRCh38 (or GRCh37) \
	--dir_cache /cache/ \
	--fasta /cache/[to reference genome fasta file] \
	--everything \
	--total_length \
	--allele_number \
	--hgvsg \
	--shift_hgvs 1 \
	--transcript_version \
	--canonical \
	--vcf \
	--flag_pick
```

For more detail, refer to Ensembl VEP documentation: https://ensembl.org/info/docs/tools/vep/script/vep_download.html#installer.


# Post-processing annotated VEP for loading into HPDS

The python script "transform_csq.v2.py" removes complex and bulky VEP annotation from VCF file and leaves only the following informations, reformatted for loading into HPDS.

The VEP annotation field in VCF can vary by exact options used in VEP annotation.
The script can detect VEP annotation format from the header line in the VCF file, if it follows the style "##INFO=<ID=CSQ... Format: ...>."

~~Thus, before running the python script:~~
~~1) Open the annotated VCF file to find the line begins with "##INFO=<ID=CSQ" and copy the string between "Format: " and "">".~~
~~2) Then open the python script and edit the line begins with 'csq_headerL' into `csq_headerL = '[string copied from annotated VCF file]'.split('|')`~~

~~Finally, run the python script.~~

`python transform_csq.v2.py [options] [VEP annotated VCF path/filename] [new filename]`

## New option: `--pick`

If present, use only the most severe consequences from VEP annotation (flagged as 'PICK', by VEP option `--flag_pick`)

## New option: `--cds`

If present, use only the variants in coding sequence (CDS). 
Specifically, this option will keep only variants whose rate of variant impact by VEP (https://ensembl.org/info/genome/variation/prediction/predicted_data.html) is not "MODIFIER."
