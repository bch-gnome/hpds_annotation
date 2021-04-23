![Docker](https://github.com/bch-gnome/hpds_annotation/workflows/Docker%20Image%20CI/badge.svg)

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

`python transform_csq.v2.py [options] [VEP annotated VCF path/filename] [new filename]`

or use docker image.

```
docker run --rm -it -v [directory for input VCF file]:/work ikarus97/hpds_annotation:latest \
	python /transform_csq.v2.py [options] /work/[input VCF path/filename] /work/[output VCF path/filename]
```

Docker image for the main script `transform_csq.v2.py` is available from [Docker Hub](https://hub.docker.com/r/ikarus97/hpds_annotation).

Image with the tag `latest` contains the most up-to-date version.

Previous versions are archived with the creation dates  as tags.

## Options

`--pick`

If present, use only the most severe consequences from VEP annotation (flagged as 'PICK', by VEP option `--flag_pick`)

`--cds`

If present, use only the variants in coding sequence (CDS). 
Specifically, this option will keep only variants whose rate of variant impact by VEP (https://ensembl.org/info/genome/variation/prediction/predicted_data.html) is not "MODIFIER."

`--vep-gnomad-af <string>`

Specify which field in VEP annotation will be extracted for gnomAD allele frequency. Use this if custom file (e.g., gnomAD genomes file) is used for gnomAD allele frequency.
Default value: `gnomAD_AF`

For example, if you want to use gnomAD genome allele frequency from the following VEP argument:
`--custom /path/to/custom/file.vcf.gz,CUSTOM_TAG,vcf,exact,0,My_Field`

then, add `--vep-gnomad-af CUSTOM_TAG_My_Field` to options to use the value of "My_Field" as gnomAD allele frequency.

`--allow-modifier`

If present, output "Variant_severity" for variants that are "MODIFIER". As of 2021-04-23, by default, such variants do not have "Variant_severity" in the INFO column to reduce overhead in HPDS.

# Change Log

## 2021-04-23

### Added
- Added option `--allow-modifier` to output tags for "MODIFIER" variants.

### Changed
- Changed the default behavior for "Variant_severity" tag: not to output the tag if the value is "MODIFIER".
- Modified header line for VCF to meet VCF specification.
