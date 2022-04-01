![Docker](https://github.com/bch-gnome/hpds_annotation/workflows/Docker%20Image%20CI/badge.svg)

This repository describes steps to prepare and annotate VCF files for loading into HPDS as in https://github.com/hms-dbmi/pic-sure-hpds-copdgene.

# Recommended steps

## 1. Pre-process VCF

The input VCF file needs to be normalized for variant representaiton: to split multi-allelic variants into separate lines in VCF.
There are various ways to do this, but one of the way is to use bcftools (http://samtools.github.io/bcftools/).

`bcftools norm -m -any -f [fasta file for referenge genome] [input VCF file] | bgzip -c > [normalized VCF file]`

Also make the normalized VCF file tabix-indexed to accelerate the next step: annotation by VEP.

`tabix -p vcf [normalized VCF file]`

(tabix is part of htslib. http://www.htslib.org/download/)


## 2. Annotating VCF using VEP ([variant effect predictor](https://www.ensembl.org/vep))

For annotating VCF with VEP, we need:
- VEP docker image (ensemblorg/ensembl-vep) or local installation
- VEP cache files from Ensembl [FTP](ftp://ftp.ensembl.org/pub/current_variation/indexed_vep_cache): **release 102 or later, GRCh38, merged for RefSeq and Ensembl.** ("homo_sapiens_merged_vep_NNN_GRCh38.tar.gz")
- GRCh38 reference fasta file from Ensembl [FTP](ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna_index).
- gnomAD genomes release v3 ([link](https://gnomad.broadinstitute.org/downloads#v3-variants)). Variant files are split by chromosome, needs to be merged into single file (e.g., `bcftools concat`).

Recommended VEP command:

```
docker run --rm -it -v [local VEP cache directory]:/cache -v [local gnomAD directory]:/gnomAD -v [directory for input VCF file]:/work ensemblorg/ensembl-vep /opt/vep/vep \
	--cache --offline --merged \
	--species homo_sapiens \
	--compress_output bgzip \
	--input_file /work/[input VCF path/filename] \
	--output_file /work/[annotated VCF path/filename] \
	--no_stats \
	--force_overwrite \
	--assembly GRCh38 \
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
	--custom /gnomAD/[to gnomAD genomes as single file],GNOMAD_G,vcf,exact,0,AF \\
	--flag_pick
```

The `--custom` option line makes "AF" values from gnomAD genomes file added to the matching variant as "GNOMAD_G_AF" value.
For more detail, refer to Ensembl VEP documentation: https://ensembl.org/info/docs/tools/vep/script/vep_download.html#installer.


## 3. Post-processing annotated VEP for loading into HPDS

The python script "transform_csq.v2.py" removes complex and bulky VEP annotation from VCF file and leaves only the following informations, reformatted for loading into HPDS.

| VEP field | Output tag | Values |
| --------- | ---------- | ------ |
| SYMBOL | Gene_with_variant | Official gene symbol for the variant |
| IMPACT | Variant_severity | The severity for the calculated consequence of the variant. Possible values: `HIGH`,`MODERATE`,`LOW`, or empty (for `MODIFIER`).  |
| Consequence | Variant_consequence_calculated | Stardardized description for the calculated consequence of the variant. |
| VARIANT_CLASS | Variant_class | Type of variant. Possible values: `SNV`,`deletion`,`insertion`. |
| GNOMAD_G_AF | Variant_frequency_in_gnomAD | Variant allele frequency from gnomAD genomes. |
| GNOMAD_G_AF | Variant_frequency_as_text | Discretized variant allele frequency from gnomAD genomes. Possble values: `Novel`, `Rare` (less than 1%), `Common` (1% or greater) |

The VEP annotation field in VCF can vary by exact options used in VEP annotation.
The script can detect VEP annotation format from the header line in the VCF file, if it follows the style "##INFO=<ID=CSQ... Format: ...>."

`python transform_csq.v2.py [options] [VEP annotated VCF path/filename] [new filename] --vep-gnomad-af GNOMAD_G_AF`

or use docker image.

```
docker run --rm -it -v [directory for input VCF file]:/work ikarus97/hpds_annotation:latest \
	python /transform_csq.v2.py [options] /work/[input VCF path/filename] /work/[output VCF path/filename] --vep-gnomad-af GNOMAD_G_AF
```

If gnomAD genomes were used with different name from `GNOMAD_G` (case sensitive), then you need to provide the correct field name with opton `--vep-gnomad-af`.
For example, if VEP was run with `--custom /[gnomad file],GG,vcf,exact,0,AF`, then the correct field name for gnomAD allele frequencies would be `GG_AF`, thus we'll need `--vep-gnomad-af GG_AF`.

Docker image for the main script `transform_csq.v2.py` is available from [Docker Hub](https://hub.docker.com/r/ikarus97/hpds_annotation).

Image with the tag `latest` contains the most up-to-date version.

Previous versions are archived with the creation dates  as tags.

### Options

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

## 2022-04-01
- Updated readme description to add more details and recommended steps.

## 2021-06-21

### Added
- Added a routine to check if the input VEP annotation includes all the required fields. If found any missing fields, the script will show error message on the missing field name and exit.

## 2021-04-23

### Added
- Added option `--allow-modifier` to output tags for "MODIFIER" variants.

### Changed
- Changed the default behavior for "Variant_severity" tag: not to output the tag if the value is "MODIFIER".
- Modified header line for VCF to meet VCF specification.
