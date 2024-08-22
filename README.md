![Docker](https://github.com/bch-gnome/hpds_annotation/workflows/Docker%20Image%20CI/badge.svg)
This repository contains scripts and outlines steps for preparing genomic data in VCF format to be loaded into [HPDS](https://github.com/hms-dbmi/pic-sure-hpds-copdgene).

# Recommended Steps for Preparing VCF Files
The VCF preparation for HPDS consists of three steps: (1) variant normalization, (2) variant annotation, and (3) post-processing of annotations.

## Note: Splitting VCF files by chromosome
Depending on available resources and the scale of the input VCF file(s), you can split the file by chromosome to speed up processing or reduce resource requirements. `bcftools` can be used to split VCF file by chromosomes: `bcftools view -O z -o [output VCF] [input VCF]`. It can be done in parallel using `xargs` as follows:

```
xargs --max-procs [# of threads] --arg-file chromosome_list -i bcftools view -O z -o [prefix].{}.vcf.gz [input VCF]
```

The `chromosome_list` file contains all chromosome names, with each name on a separate line. An example of this file is available in the repository.

Also, it can be done at any stage before or after steps 1–3, or integrated with step 1 as shown below:

```
xargs --max-procs [# of threads] --arg-file chromosome_list -i bcftools norm -m -any -f [FASTA] -O z -o [prefix].{}.vcf.gz [input VCF]
```

## 1. Variant normalization

The Input VCF file needs normalization for variant representation: to split multi-allelic variants into separate lines. While there are various methods to do this, here's an example using [bcftools](http://samtools.github.io/bcftools):

```
bcftools norm -m -any -f [reference genome FASTA file] [input VCF] | bgzip -c > [normalized VCF]
```

It is recommended to build a tabix index for the output file: `tabix -p vcf [normalized VCF]`.

## 2. Variant annotation using [VEP](https://www.ensembl.org/vep) (Variant Effect Predictor)

Requirements for this step:
- VEP Docker image (ensemblorg/ensembl-vep) or local installation.
  - Could be used with Singularity or similar.

- VEP cache files from [Ensembl FTP](https://ftp.ensembl.org/pub) (we recommend using the indexed cache, [url](https://ftp.ensembl.org/pub/current_variation/indexed_vep_cache/) for the latest version): **release 102 or later, GRCh38, merged for both RefSeq and Ensembl**.

- GRCh38 reference genome FASTA file, preferably from Ensembl FTP.
  - For the command below, we assume the file is in the root directory of the VEP cache files.

- gnomAD release v3 (genomes) or v4 (joint frequency).
  - These files are divided by chromosome and need merging into a single file (e.g., `bcftools concat`).

Recommended VEP command (using Docker):

```
docker run --rm -it -v [local VEP cache directory]:/cache \
	-v [local gnomAD directory]:/gnomAD \
	-v [input VCF directory]:/work \
	ensemblorg/ensembl-vep /opt/vep/vep \
	--cache --offline --compress_output bgzip --dir_cache /cache/ \
	--species homo_sapiens --merged --assembly GRCh38 \
	--input_file /work/[input VCF] \
	--output_file /work/[output VCF] \
	--fasta /cache/[to reference genome FASTA file] \
	--no_stats --force_overwrite --everything --vcf \
	--total_length --allele_number --hgvsg --transcript_version --canonical --shift_hgvs 1 \
	--custom /gnomAD/[gnomAD file as a single file],GNOMAD_G,vcf,exact,0,AF \
	--flag_pick
```

Note:
- The arguments `GNOMAD_G` and `AF` used with `--custom` will be relevant in the next step.

- For gnomAD release v4 files with joint frequency, use `AF_joint` as the last argument of the `--custom` option. The full option would be: `--custom /gnomAD/[gnomAD file],GNOMAD_G,vcf,exact,0,AF_joint`.

- The `--flag_pick` option will impact in the next step.

## 3. Post-processing of VEP annotations

The Python script `transform_csq.v3.py` in this repository simplifies VCF files by removing complex VEP annotations. It preserves only essential information (detailed in the table below), reformatting it for HPDS loading.

| VEP field | Output tag | Values |
| --------- | ---------- | ------ |
| SYMBOL | Gene_with_variant | The official gene symbol associated with the variant |
| IMPACT | Variant_severity | The severity of the variant's calculated consequence. Possible values: `HIGH`, `MODERATE`, `LOW`, or empty (for `MODIFIER`). |
| Consequence | Variant_consequence_calculated | A standardized description of the variant's calculated consequence. |
| VARIANT_CLASS | Variant_class | Varaint type. Possible values: `SNV`, `deletion`, or `insertion`. |
| GNOMAD_G_AF | Variant_frequency_in_gnomAD | Variant allele frequency from gnomAD. |
| GNOMAD_G_AF | Variant_frequency_as_text | Discretized variant allele frequency from gnomAD. Possible values: `Novel` (not found in gnomAD), `Rare` (less than 1% frequency), or `Common` (1% frequency or greater). |

The VEP annotation field in VCF files can vary depending on the specific options used during VEP annotation. The script automatically detects the VEP annotation format from the VCF header, but only if it follows the standard format: `##INFO=<ID=CSQ... Format: ...>`.

### Using the Python script

```
python3 transform_csq.v3.py -R [reference genome FASTA file] --vep-gnomad-af GNOMAD_G_AF [other options] [VEP annotated VCF] [new filename]
```

The argument for `--vep-gnomad-af` in the Python script should match the gnomAD file arguments used with VEP's `--custom` option. For example, if VEP was run with `--custom [gnomAD file],GG,vcf,exact,0,AF`, you should use `--vep-gnomad-af GG_AF` in the Python script to correctly retrieve gnomAD allele frequency values.

The Python script also performs several checks to ensure the input VCF file is compatible with GRCh38:
1. It checks if the provided FASTA file is compatible with GRCh38 by checking the chromosome naming style (chr1, chr2, etc.), chromosome lengths, and MD5 checksums for sequence data. This check covers autosomes (chr1 ~ chr22) and sex chromosomes (chrX and chrY) only.

2. It checks if the chromosome names in the input VCF file follow the GRCh38 style (chr1, chr2, etc.).

3. It checks the compatibility between the input VCF and FASTA files. For this, it compares the reference sequence bases in VCF and FASTA files at 100 random positions.

The script also excludes variants on chromosomes other than autosomes or sex chromosomes from the output.

### Note: Handling non-coding variants

Non-coding variants (specifically, variants with `MODIFIER` as `IMPACT` by VEP) are often of less interest. As of 2024-08-22, the script offers three modes (`--mode`) for handling these variants:
1. `all`: Processes all variants regardless of their impact. This is the default mode.

2. `cds_only`: Outputs only variants whose `IMPACT` is not `MODIFIER`.

3. `cds_rsid`: Similar to `cds_only`, but also includes `MODIFIER` variants if their `ID` column in the VCF is not empty.
   - This assumes RSIDs from [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) are annotated in the VCF's `ID` column. You can add these annotations using bcftools: `bcftools annotate -x ID -c +ID -a [dbSNP file] [input VCF]`.

### Note: bcftools versions

bcftools 1.16 had an issue where some FORMAT fields were improperly handled - left empty - when missing for certain samples. This problem didn't occur if all samples had values or if none had values. The Python script later translated these empty fields as '0' (the default behavior of the pysam library used in the script). This issue was resolved in version 1.17 and later.

To check if your VCF is affected by this issue, you can compare VCFs before and after running the Python script as follows:

```
zcat [VCF before script] | grep -v '^#' | cut -f 10- | md5sum
zcat [VCF after script] | grep -v '^#' | cut -f 10- | md5sum
```

If the two checksums differ, it indicates that FORMAT fields were affected by the issue.

### Python script options

`-R <path to FASTA file>`

*Required*. Specifies the path to the GRCh38 reference FASTA file. The corresponding index file (.fai) should be located in the same directory.

`--pick`

If specified, this option uses only the most severe consequences from VEP annotations (marked as 'PICK'). It requires the VEP option `--flag-pick` to be used during annotation.

`--cds`

If specified, this option outputs only variants in coding sequences (CDS) (variants whose impact is not `MODIFIER`). It **cannot** be used with `--mode` all or `--mode cds_rsid`. In future releases, this option will be retired because it is the same as `--mode cds_only`.

`--vep-gnomad-af <string>`

Specify the field in VEP annotations to be used as the variant allele frequency from gnomAD. This option requires the VEP `--custom` option with gnomAD release files. The default value is `gnomAD_AF`.

For example, if you used the following option in VEP:

`--custom /path/to/gnomAD/file.vcf.gz,CUSTOM_TAG,vcf,exact,0,My_field`

you should add `--vep-gnomad-af CUSTOM_TAG_My_field` (case sensitive) to use the value of "My_field" as the gnomAD variant allele frequency.

`--mode`

Specifies which variants to include in the output file. Options are: `all` (includes all variants), `cds_rsid` (coding variants and those with RSIDs), or `cds_only` (only coding variants). The default is `all`.

1. `all`: Processes all variants regardless of their impact. This is the default mode.

2. `cds_only`: Outputs only variants whose `IMPACT` is not `MODIFIER`.

3. `cds_rsid`: Similar to `cds_only`, but also includes `MODIFIER` variants if their `ID` column in the VCF file is not empty.

`--allow-modifier-tag`

Previously named as `--allow-modifier`. If specified, this option outputs the `Variant_severity` value even for `MODIFIER` variants. Note that this option **does not** determine whether such variants are included in the output file. As of 2021-04-23, by default, `MODIFIER` variants lack a `Variant_severity` value in the `INFO` column.


# Change Log

## 2024-08-22
Summary: added `--mode` option to manage `MODIFIER` variants. Updated script to eliminate redundant output lines and exclude variants outside autosomes or sex chromosomes. Renamed options for clarity.

### Added
- Added a `--mode` option to control the handling of `MODIFIER` variants in the output file. For the `cds_rsid` mode, RSIDs for dbSNP variants must be annotated in the `ID` column of the VCF.

### Changed
- When printing different consequences for the same variant, the script checks for redundancy based on (Symbol, Consequence) pair. This check is performed line by line. As a result, some redundancy may still exist if there were redundant lines in the original VCF.
- The option `--allow-modifier` has been renamed to `--allow-modifier-tag` for clarity.
- Regardless of the `--mode` setting, variants outside of autosomes and sex chromosomes will not be written to the output.

## 2022-12-06
Summary: major re-writing of code using python3 and pysam library. Also added a few basic check-up routines to ensure the compatibility with GRCh38.

### Added
- Added a *required* option (`-R`) to specify the path to GRCh38 reference FASTA file. It also requires the index file (.fai) in the same place as FASTA.
- Added routines to check the following:
	- if the chromosome names and sequences in the FASTA file is compatible with GRCh38.
	- if the chromosome names in the input VCF is compatible with GRCh38.
	- if the coordinates in the input VCF is compatible with the FASTA.

### Changed
- Code is written with python3 and pysam library. The old python2 scrypt (transform_csq.v2.py) is moved under "legacy_script/".
- The name of main script is changed to "transform_csq*v3*.py".
- Now the script requires tabix-index for input VCF (.tbi) in the same place as VCF.
- The script prints more descriptive error or information messages during the process.

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
