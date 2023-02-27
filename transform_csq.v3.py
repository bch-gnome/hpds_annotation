import os, sys, pysam, hashlib, logging, random, math, re
from optparse import OptionParser

N_TEST = 100
VEP_TAG = "CSQ"

## lengths & checksums for GRCh38 chromosomes (autosomes + X,Y). Taken from Homo_sapiens_assembly38.fasta (downloaded from Broad bucket)
b38 = {
	'chr1': {'len': 248956422, 'md5': '6aef897c3d6ff0c78aff06ac189178dd'},
	'chr2': {'len': 242193529, 'md5': 'f98db672eb0993dcfdabafe2a882905c'},
	'chr3': {'len': 198295559, 'md5': '76635a41ea913a405ded820447d067b0'},
	'chr4': {'len': 190214555, 'md5': '3210fecf1eb92d5489da4346b3fddc6e'},
	'chr5': {'len': 181538259, 'md5': 'a811b3dc9fe66af729dc0dddf7fa4f13'},
	'chr6': {'len': 170805979, 'md5': '5691468a67c7e7a7b5f2a3a683792c29'},
	'chr7': {'len': 159345973, 'md5': 'cc044cc2256a1141212660fb07b6171e'},
	'chr8': {'len': 145138636, 'md5': 'c67955b5f7815a9a1edfaa15893d3616'},
	'chr9': {'len': 138394717, 'md5': '6c198acf68b5af7b9d676dfdd531b5de'},
	'chr10': {'len': 133797422, 'md5': 'c0eeee7acfdaf31b770a509bdaa6e51a'},
	'chr11': {'len': 135086622, 'md5': '1511375dc2dd1b633af8cf439ae90cec'},
	'chr12': {'len': 133275309, 'md5': '96e414eace405d8c27a6d35ba19df56f'},
	'chr13': {'len': 114364328, 'md5': 'a5437debe2ef9c9ef8f3ea2874ae1d82'},
	'chr14': {'len': 107043718, 'md5': 'e0f0eecc3bcab6178c62b6211565c807'},
	'chr15': {'len': 101991189, 'md5': 'f036bd11158407596ca6bf3581454706'},
	'chr16': {'len': 90338345, 'md5': 'db2d37c8b7d019caaf2dd64ba3a6f33a'},
	'chr17': {'len': 83257441, 'md5': 'f9a0fb01553adb183568e3eb9d8626db'},
	'chr18': {'len': 80373285, 'md5': '11eeaa801f6b0e2e36a1138616b8ee9a'},
	'chr19': {'len': 58617616, 'md5': '85f9f4fc152c58cb7913c06d6b98573a'},
	'chr20': {'len': 64444167, 'md5': 'b18e6c531b0bd70e949a7fc20859cb01'},
	'chr21': {'len': 46709983, 'md5': '974dc7aec0b755b19f031418fdedf293'},
	'chr22': {'len': 50818468, 'md5': 'ac37ec46683600f808cdd41eac1d55cd'},
	'chrX': {'len': 156040895, 'md5': '2b3a55ff7f58eb308420c8a9b11cac50'},
	'chrY': {'len': 57227415, 'md5': 'ce3e31103314a704255f3cd90369ecce'}
}

## Fields for INFO column in the output VCF (HPDS-ready)
## TODO: will be changed as HPDS fields change
out_columnH = {
	'Gene_with_variant': {
		'Number': '.',
		'Type': 'String',
		'Description': "The official symbol for a gene affected by a variant.",
		'orgVEP': 'SYMBOL'
	},
	'Variant_severity': {
		'Number': '.',
		'Type': 'String',
		'Description': "The severity for the calculated consequence of a variant on a gene. Possible values: HIGH (frameshift, splice disrupting, or truncating variants), MODERATE (non-frameshift insertions or deletions, variants altering protein sequencing without affecting its length), LOW (other coding sequence variants including synonymous variants).",
		'orgVEP': 'IMPACT'
	},
	'Variant_consequence_calculated': {
		'Number': '.',
		'Type': 'String',
		'Description': "A standardized term from the Sequence Ontology (http://www.sequenceontology.org) to describe the calculated consequence of a variant.",
		'orgVEP': 'Consequence'
	},
	'Variant_class': {
		'Number': 1,
		'Type': 'String',
		'Description': "A standardized term from the Sequence Ontology (http://www.sequenceontology.org) to describe the type of a variant. Possible values: SNV, deletion, insertion.",
		'orgVEP': 'VARIANT_CLASS'
	},
	'Variant_frequency_in_gnomAD': {
		'Number': 1,
		'Type': 'Float',
		'Description': "The variant allele frequency in gnomAD combined population.",
		'orgVEP': 'gnomAD_AF'
	},
	'Variant_frequency_as_text': {
		'Number': 1,
		'Type': 'String',
		'Description': "The variant allele frequency in gnomAD combined population as discrete text categories. Possible values: Novel, Rare (variant frequency less than 1%), Common (variant frequency greater than or equal to 1%).",
		'orgVEP': 'gnomAD_AF'
	}
}

def check_fasta_for_genome_version(in_filename):
	logger = logging.getLogger()
	logger.info("Start checking if the reference FASTA is GRCh38.")

	input_ref = pysam.FastaFile(in_filename)

	b38_names = list(b38.keys())
	b38_names.sort()
	input_names = input_ref.references
	input_names.sort()

	## check if all b38_name is found in input_names (not necessary equal. FASTA can contain additional contigs..., but they should at least contain 1-22, X, Y)
	res = all(chrom in input_names for chrom in b38_names)
	if not res:
		logger.error("==========================================================")
		logger.error("The chromosome names in FASTA are not compatible with GRCh38, or there are chromosomes missing.")
		logger.error(" - Chromosome names in FASTA: %s." % input_names)
		logger.error(" - Chromosome names in GRCh38: %s." % b38_names)
		logger.error("Exiting...")
		logger.error("==========================================================")
		sys.exit(1)

	for chrom in b38:
		## compare chromosome length
		ref_size = input_ref.get_reference_length(chrom)
		if ref_size != b38[chrom]['len']:
			logger.error("==========================================================")
			logger.error("The chromosome length in FASTA does not match with GRCh38: %s" % chrom)
			logger.error(" - %s in FASTA vs. %s in GRCh38." % (ref_size, b38[chrom]['len']))
			logger.error("Exiting...")
			logger.error("==========================================================")
			sys.exit(1)

		## compare MD5 checksum value
		ref_md5 = hashlib.md5(input_ref.fetch(chrom).encode('utf-8'))
		if ref_md5.hexdigest() != b38[chrom]['md5']:
			logger.error("==========================================================")
			logger.error("The sequence in FASTA has different MD5 checksum than GRCh38: %s" % chrom)
			logger.error(" - %s in FASTA vs. %s in GRCh38." % (ref_md5.hexdigest(), b38[chrom]['md5']))
			logger.error("Exiting...")
			logger.error("==========================================================")
			sys.exit(1)

	logger.info("Finished checking if the reference FASTA is GRCh38.")

def check_vcf_chromosome_name(in_filename):
	logger = logging.getLogger()
	logger.info("Start checking chromosome names in input VCF.")

	input_vcf = pysam.VariantFile(in_filename)
	input_names = list(input_vcf.header.contigs)
	input_names.sort()

	b38_names = list(b38.keys())
	b38_names.sort()

	## just check if there's intersection between b38_names and input_names
	## we can't assume input VCF will always have all chromosomes (1-22,X,Y)
	res = any(chrom in input_names for chrom in b38_names)

	if res:
		#info message for continuing with 'intersecting' chromosomes
		common_names = [chrom for chrom in b38_names if chrom in input_names]
		logger.info("Continuing with the following GRCh38 chromosome names found in VCF.")
		logger.info(" - %s" % common_names)
	else:
		logger.error("==========================================================")
		logger.error("The chromosome names in input VCF are not compatible with GRCh38.")
		logger.error(" - Chromosome names in input VCF: %s" % input_names)
		logger.error(" - Chromosome names in GRCh38: %s" % b38_names)
		logger.error("Exiting...")
		logger.error("==========================================================")
		sys.exit(1)

	logger.info("Finished checking chromosome names in input VCF.")

def check_vcf_with_FASTA(infile_vcf, infile_fasta):
	logger = logging.getLogger()
	logger.info("Start checking if the given reference FASTA is compatible with the input VCF file.")

	input_vcf = pysam.VariantFile(infile_vcf)
	input_fasta = pysam.FastaFile(infile_fasta)

	## random select 1 million base intervals (until it counts N_TEST
	n_test = 0
	b38_names = list(b38.keys())
	while n_test < N_TEST:
		test_chr = random.choice(b38_names)
		test_start = random.randint(0, math.floor(b38[test_chr]['len']/1000000))*1000000
		test_end = test_start + 1000000
		logger.info(" - Checking random variants in %s:%s-%s" % (test_chr, test_start, test_end))
		tt = input_vcf.fetch(contig=test_chr, start=test_start, stop=test_end)
		for record in tt:
			if random.randrange(100) < 10:
				vcf_base = record.ref
				fasta_base = input_fasta.fetch(reference=test_chr, start = record.start, end = record.stop)
				n_test += 1

				if vcf_base.upper() != fasta_base.upper():
					logger.error("==========================================================")
					logger.error("The VCF file is not compatible with the given FASTA.")
					logger.error(" - Mismatch found in %s:%s. %s in VCF vs. %s in FASTA." % (test_chr, record.pos, vcf_base, fasta_base))
					logger.error("Exiting...")
					logger.error("==========================================================")
					sys.exit(1)
			##if
		#for record
	## while

	logger.info("Finished checking if the given reference FASTA is compatible with the input VCF file.")

def update_info(infile, outfile, options):
	logger = logging.getLogger()
	logger.info("Start updating INFO fields.")

	input_vcf = pysam.VariantFile(infile)

	## check input header if CSQ field exists.
	## And take the list of VEP annotation fields
	flag = False
	csq_headerL = []
	for rec in input_vcf.header.records:
		if rec.type == 'INFO' and rec['ID'] == VEP_TAG:
			flag = True
			csq_headerL = re.match('".*Format: (.*)"', rec['Description']).groups()[0].split('|')
			break

	## exit if CSQ is not found in the header
	if not flag:
		logger.error("==========================================================")
		logger.error("The input VCF doesn't have VEP information in the header.")
		logger.error("Exiting...")
		logger.error("==========================================================")
		sys.exit(1)

	## check if all of necessary VEP annotation fields (in out_columnH) are there.
	for key in out_columnH:
		if not out_columnH[key]['orgVEP'] in csq_headerL:
			logger.error("==========================================================")
			logger.error("THE VEP annotation in VCF is missing the required field: %s" % out_columnH[key]['orgVEP'])
			logger.error("Exiting...")
			logger.error("==========================================================")
			sys.exit(1)

	in_header = input_vcf.header
	out_header = pysam.VariantHeader()

	## copy header from input VCF, except for INFO fields.
	for rec in input_vcf.header.records:
		if rec.type != 'INFO':
			out_header.add_record(rec)
	## add INFO fields for HPDS-annotation
	for key in out_columnH:
		out_header.add_meta('INFO', items=[('ID', key), ('Number', out_columnH[key]['Number']), ('Type', out_columnH[key]['Type']), ('Description', out_columnH[key]['Description'])])
		in_header.add_meta('INFO', items=[('ID', key), ('Number', out_columnH[key]['Number']), ('Type', out_columnH[key]['Type']), ('Description', out_columnH[key]['Description'])])

	## copy list of samples to new header
	for samp in input_vcf.header.samples:
		out_header.add_sample(samp)

	output_vcf = pysam.VariantFile(outfile, "w", header = out_header)

	## main loop for all variant records in the input VCF
	for rec in input_vcf:
		## skip the line without any VEP annotation (happens with strange ALT values such as '*')
		if not VEP_TAG in list(rec.info):
			continue

		csqL = list(rec.info[VEP_TAG])

		## loop through multiple annotations (for different transcripts)
		for csq in csqL:
			new_rec = rec.copy()
			new_rec.info.clear()
			csq_valL = csq.split('|')

			if options.cds_only and csq_valL[ csq_headerL.index('IMPACT') ] == 'MODIFIER':
				continue
			if options.pick_only and csq_valL[ csq_headerL.index('PICK') ] == '':
				continue

			for key in out_columnH:
				## first take VEP annotation field value
				csq_val = csq_valL[ csq_headerL.index(out_columnH[key]['orgVEP']) ]

				## For 'Consequence', split complex terms by '&' and take only the first one
				if key == 'Variant_consequence_calculated':
					csq_val = csq_val.split('&')[0]

				## For gnomAD AFs, convert to float. If empty ('.' or ''), assign -10.
				if out_columnH[key]['orgVEP'] == options.vep_gnomad_af:
					if csq_val == '.' or csq_val == '':
						csq_val = -10
					else:
						csq_val = float(csq_val)

				## Select appropriate text description for gnomAD AF based on the previous step.
				if key == 'Variant_frequency_as_text':
					if float(csq_val) < 0:
						csq_val = 'Novel'
					elif float(csq_val) < 0.01:
						csq_val = 'Rare'
					else:
						csq_val = 'Common'
					
				## If the value is still not empty, add it into INFO field for new record.
				if csq_val != '.' and csq_val != '':
					## But if the field is 'Variant_severity' and the value is 'MODIFIER', don't add it unless modifier is allowed by --allow-modifier.
					if key == 'Variant_severity' and csq_val == 'MODIFIER' and not options.allow_modifier:
						continue
					new_rec.info[key] = csq_val

			## make the new record fit for new header
			new_rec.translate(out_header)
			## output the record to file
			output_vcf.write(new_rec)
		#for csq
	#for rec in input_vcf
	output_vcf.close()

	logger.info("Finished updating INFO fields.")

if __name__ == '__main__':
	## set up logging
	logging.basicConfig(format='[%(asctime)s] %(message)s', level='INFO')

	## set up input options
	usage = "usage: %prog [options] -R [reference FASTA file] [input VCF file (needs indexed)] [output VCF file]"
	parser = OptionParser(usage)
	parser.add_option("-R", action="store", type="string", dest="ref_fasta", help="(required) FASTA file for the reference genome version as used in VCF. Needs .fai index.")
	parser.add_option("--pick", action="store_true", dest="pick_only", help="If asserted, picks the most severe consequence only (marked by VEP). [default: %default]")
	parser.add_option("--cds", action="store_true", dest="cds_only", help="If asserted, keeps variants in coding region (CDS) only, i.e., VEP IMPACT is not 'MODIFIER'. [default: %default]")
	parser.add_option("--vep-gnomad-af", action="store", type="string", dest="vep_gnomad_af", default="gnomAD_AF", help="VEP field name to be used for 'Variant_frequency_in_gnomAD' and 'Variant_frequency_as_text'. [default: %default]")
	parser.add_option("--allow-modifier", action="store_true", dest="allow_modifier", help="If asserted, writes 'Variant_severity' INFO column for variants with 'MODIFIER' as 'VARIANT_IMPACT' by VEP (mostly intronic or intergenic variants). Not recommended. [default: %default]")
	
	(options, args) = parser.parse_args()

	## check input options
	if len(args) < 2:
		parser.error("Need input & output file names!")
		sys.exit(1)
	if options.ref_fasta == None:
		parser.error("Reference FASTA file (-R) is required!")
		sys.exit(1)

	if not os.path.exists("%s.fai" % options.ref_fasta):
		parser.error("Needs index file (.fai) for the FASTA!")
		sys.exit(1)

	if not os.path.exists("%s.tbi" % args[0]):
		parser.error("Needs index file (.tbi) for the input VCF!")
		sys.exit(1)

	if options.vep_gnomad_af:
		out_columnH['Variant_frequency_in_gnomAD']['orgVEP'] = options.vep_gnomad_af
		out_columnH['Variant_frequency_as_text']['orgVEP'] = options.vep_gnomad_af

	## 1. check if the given FASTA is compatible with GRCh38
	##    - the function will exit with 1 upon error
	check_fasta_for_genome_version(options.ref_fasta)

	## 2. check if the given VCF has GRCh38-compatible chromosome names
	check_vcf_chromosome_name(args[0])

	## 3. check if random loci in VCF matches with the given FASTA
	check_vcf_with_FASTA(args[0], options.ref_fasta)

	## 4. MAIN routine
	update_info(args[0], args[1], options)
