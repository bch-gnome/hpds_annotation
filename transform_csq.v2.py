#!/usr/bin/python

import os,sys,re
import gzip
from optparse import OptionParser

csq_headerL='Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|REFSEQ_MATCH|SOURCE|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|HGVSg|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|DownstreamProtein|ProteinLengthChange|TSSDistance|CSN|SpliceRegion|FATHMM_pred|FATHMM_score|MutationAssessor_pred|MutationAssessor_score|MutationTaster_pred|MutationTaster_score|PROVEAN_pred|PROVEAN_score|REVEL_score|ada_score|rf_score|CADD_PHRED|CADD_RAW|Condel|LoFtool|ExACpLI|HGMD|HGMD_CLASS|HGMD_MUT|HGMD_PHEN|HGMD_RANKSCORE|GNOMAD_G|GNOMAD_G_AF_AFR|GNOMAD_G_AF_AMR|GNOMAD_G_AF_ASJ|GNOMAD_G_AF_EAS|GNOMAD_G_AF_FIN|GNOMAD_G_AF_NFE|GNOMAD_G_AF_OTH|phastCons100|phyloP100|RMSK|GERP|COV_GNOMAD_E|COV_GNOMAD_G'.split('|')

out_columnH = {
    'SYMBOL': {'Name': 'Gene_with_variant', 'Number': 1, 'Type': 'String', 'Description': "The official symbol for a gene affected by a variant."},
    'IMPACT': {'Name': 'Variant_severity', 'Number': 1, 'Type': 'String', 'Description': "The severity for the calculated consequence of a variant on a gene. Possible values: HIGH (frameshift, splice disrupting, or truncating variants), MODERATE (non-frameshift insertions or deletions, variants altering protein sequencing without affecting its length), LOW (other coding sequence variants including synonymous variants), MODIFIER (all others)."},
    'Consequence': {'Name': 'Variant_consequence_calculated', 'Number': 1, 'Type': 'String', 'Description': "A standardized term from the Sequence Ontology (http://www.sequenceontology.org) to describe the calculated consequence of a variant."},
    'VARIANT_CLASS': {'Name': 'Variant_class', 'Number': 1, 'Type': 'String', 'Description': "A standardized term from the Sequence Ontology (http://www.sequenceontology.org) to describe the type of a variant. Possible values: SNV, deletion, insertion."},
    'gnomAD_AF': {'Name': 'Variant_frequency_in_gnomAD', 'Number': 1, 'Type': 'Float', 'Description': "The variant allele frequency in gnomAD combined population."}
}

csq_infoH = {
    'Allele': {'Number': 1, 'Type':'String', 'Description': "The variant allele used to calculate the consequence"},
    'Consequence': {'Number': 1, 'Type':'String', 'Description': "The consequence type of this variant"},
    'IMPACT': {'Number': 1, 'Type': 'String', 'Description': "The impact modifier for the consequence type"},
    'SYMBOL': {'Number': '.', 'Type': 'String', 'Description': "The gene symbol"},
    'Gene': {'Number': '.', 'Type': 'String', 'Description': "Ensembl stable ID of affected gene"},
    'Feature_type': {'Number': '.', 'Type': 'String', 'Description': "The type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature."},
    'Feature': {'Number': '.', 'Type': 'String', 'Description': "Ensembl stable ID of feature"},
    'BIOTYPE': {'Number': '.', 'Type': 'String', 'Description': "Biotype of transcript or regulatory feature"},
    'EXON': {'Number': '.', 'Type': 'String', 'Description': "The exon number out of total number"},
    'INTRON': {'Number': '.', 'Type': 'String', 'Description': "The intron number out of total number"},
    'HGVSc':{'Number': '.', 'Type': 'String', 'Description': "The HGVS coding sequence name"},
    'HGVSp':{'Number': '.', 'Type': 'String', 'Description': "The HGVS protein sequence name"},
    'cDNA_position':{'Number': '.', 'Type': 'String', 'Description': "Relative position of base pair in cDNA sequence"},
    'CDS_position':{'Number': '.', 'Type': 'String', 'Description': "Relative position of base pair in coding sequence"},
    'Protein_position':{'Number': '.', 'Type': 'String', 'Description': "Relative position of amino acid in protein"},
    'Amino_acids':{'Number': '.', 'Type': 'String', 'Description': "Reference and variant amino acids"},
    'Codons':{'Number': '.', 'Type': 'String', 'Description': "Reference and variant codon sequence"},
    'Existing_variation':{'Number': '.', 'Type': 'String', 'Description': "Identifier(s) of co-located known variants"},
    'ALLELE_NUM':{'Number': '.', 'Type': 'Integer', 'Description': "Allele number from input; 0 is reference, 1 is first alternate, etc."},
    'DISTANCE':{'Number': '.', 'Type': 'Integer', 'Description': "Shortest distance from variant to transcript"},
    'STRAND':{'Number': '.', 'Type': 'Integer', 'Description': "The DNA strand (1 or -1) on which the transcript/feature lies"},
    'FLAGS':{'Number': '.', 'Type': 'String', 'Description': "Transcript quality flags"},
    'PICK':{'Number': '.', 'Type': 'Integer', 'Description': "Indicates if this block of consequence data was picked by --flag_pick or --flag_pick_allele"},
    'VARIANT_CLASS':{'Number': '.', 'Type': 'String', 'Description': "Sequence Ontology variant class"},
    'SYMBOL_SOURCE':{'Number': '.', 'Type': 'String', 'Description': "The source of the gene symbol"},
    'HGNC_ID':{'Number': '.', 'Type': 'String', 'Description': "HGNC ID of the gene symbol"},
    'CANONICAL':{'Number': '.', 'Type': 'String', 'Description': "A flag indicating if the transcript is denoted as the canonical transcript for this gene"},
    'TSL':{'Number': '.', 'Type': 'String', 'Description': "Transcript support level."},
    'APPRIS':{'Number': '.', 'Type': 'String', 'Description': "Annotates alternatively spliced transcripts as primary or alternate based on a range of computational methods."},
    'CCDS':{'Number': '.', 'Type': 'String', 'Description': "The CCDS identifer for this transcript, where applicable"},
    'ENSP':{'Number': '.', 'Type': 'String', 'Description': "The Ensembl protein identifier of the affected transcript"},
    'SWISSPROT':{'Number': '.', 'Type': 'String', 'Description': "Best match UniProtKB/Swiss-Prot accession of protein product"},
    'TREMBL':{'Number': '.', 'Type': 'String', 'Description': "Best match UniProtKB/TrEMBL accession of protein product"},
    'UNIPARC':{'Number': '.', 'Type': 'String', 'Description': "Best match UniParc accession of protein product"},
    'REFSEQ_MATCH':{'Number': '.', 'Type': 'String', 'Description': "The RefSeq transcript match status; contains a number of flags indicating whether this RefSeq transcript matches the underlying reference sequence and/or an Ensembl transcript"},
    'SOURCE':{'Number': '.', 'Type': 'String', 'Description': "The source of the gene symbol"},
    'GIVEN_REF':{'Number': '.', 'Type': 'String', 'Description': "Reference allele from input"},
    'USED_REF':{'Number': '.', 'Type': 'String', 'Description': "Reference allele as used to get consequences"},
    'BAM_EDIT':{'Number': '.', 'Type': 'String', 'Description': "Indicates success or failure of edit using BAM file"},
    'GENE_PHENO':{'Number': '.', 'Type': 'String', 'Description': "Indicates if overlapped gene is associated with a phenotype, disease or trait"},
    'SIFT':{'Number': '.', 'Type': 'String', 'Description': "The SIFT prediction and score as prediction(score)"},
    'PolyPhen':{'Number': '.', 'Type': 'String', 'Description': "The PolyPhen prediction and score as prediction(score)"},
    'DOMAINS':{'Number': '.', 'Type': 'String', 'Description': "The source and identifer of any overlapping protein domains"},
    'miRNA':{'Number': '.', 'Type': 'String', 'Description': "The identifier of affected miRNA, where applicable."},
    'HGVS_OFFSET':{'Number': '.', 'Type': 'String', 'Description': "Indicates by how many bases the HGVS notations for this variant have been shifted"},
    'HGVSg':{'Number': '.', 'Type': 'String', 'Description': "The HGVS genomic sequence name"},
    'AF': {'Number': '.', 'Type': 'String', 'Description': "Frequency of existing variant in 1000 Genomes"},
    'AFR_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in 1000 Genomes combined African population"},
    'AMR_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in 1000 Genomes combined American population"},
    'EAS_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in 1000 Genomes combined East Asian population"},
    'EUR_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in 1000 Genomes combined European population"},
    'SAS_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in 1000 Genomes combined South Asian population"},
    'AA_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in NHLBI-ESP African American population"},
    'EA_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in NHLBI-ESP European American population"},
    'gnomAD_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in gnomAD exomes combined population"},
    'gnomAD_AFR_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in gnomAD exomes African/American population"},
    'gnomAD_AMR_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in gnomAD exomes American population"},
    'gnomAD_ASJ_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in gnomAD exomes Ashkenazi Jewish population"},
    'gnomAD_EAS_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in gnomAD exomes East Asian population"},
    'gnomAD_FIN_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in gnomAD exomes Finnish population"},
    'gnomAD_NFE_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in gnomAD exomes Non-Finnish European population"},
    'gnomAD_OTH_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in gnomAD exomes combined other combined populations"},
    'gnomAD_SAS_AF':{'Number': '.', 'Type': 'Float', 'Description': "Frequency of existing variant in gnomAD exomes South Asian population"},
    'MAX_AF':{'Number': '.', 'Type': 'Float', 'Description': "Maximum observed allele frequency in 1000 Genomes, ESP and gnomAD"},
    'MAX_AF_POPS':{'Number': '.', 'Type': 'String', 'Description': "Populations in which maximum allele frequency was observed"},
    'CLIN_SIG':{'Number': '.', 'Type': 'String', 'Description': "ClinVar clinical significance of the dbSNP variant"},
    'SOMATIC':{'Number': '.', 'Type': 'String', 'Description': "Somatic status of existing variant(s); multiple values correspond to multiple values in the Existing_variation field"},
    'PHENO':{'Number': '.', 'Type': 'String', 'Description': "Indicates if existing variant is associated with a phenotype, disease or trait; multiple values correspond to multiple values in the Existing_variation field"},
    'PUBMED':{'Number': '.', 'Type': 'String', 'Description': "Pubmed ID(s) of publications that cite existing variant"},
    'MOTIF_NAME':{'Number': '.', 'Type': 'String', 'Description': "The source and identifier of a transcription factor binding profile aligned at this position"},
    'MOTIF_POS':{'Number': '.', 'Type': 'String', 'Description': "The relative position of the variation in the aligned TFBP"},
    'HIGH_INF_POS':{'Number': '.', 'Type': 'String', 'Description': "A flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP)"},
    'MOTIF_SCORE_CHANGE':{'Number': '.', 'Type': 'String', 'Description': "The difference in motif score of the reference and variant sequences for the TFBP"},
    'DownstreamProtein':{'Number': '.', 'Type': 'String', 'Description': "The downstream effects of a frameshift variant on the protein sequence of a transcript"},
    'ProteinLengthChange':{'Number': '.', 'Type': 'String', 'Description': ""},
    'TSSDistance':{'Number': '.', 'Type': 'String', 'Description': "The distance from the transcription start site for upstream variants"},
    'CSN':{'Number': '.', 'Type': 'String', 'Description': "Clinical Sequencing Nomenclature (CSN) for the variant"},
    'SpliceRegion':{'Number': '.', 'Type': 'String', 'Description': "More granular predictions of splicing effects for splicing variant"},
    'FATHMM_pred':{'Number': '.', 'Type': 'String', 'Description': "The deleteriousness prediction from FATHMM"},
    'FATHMM_score':{'Number': '.', 'Type': 'String', 'Description': "The deleteriousness score from FATHMM"},
    'MutationAssessor_pred':{'Number': '.', 'Type': 'String', 'Description': "The deleteriousness prediction from MutationAssessor"},
    'MutationAssessor_score':{'Number': '.', 'Type': 'String', 'Description': "The deleteriousness score from MutationAssessor"},
    'MutationTaster_pred':{'Number': '.', 'Type': 'String', 'Description': "The deleteriousness prediction from MutationTaster"},
    'MutationTaster_score':{'Number': '.', 'Type': 'String', 'Description': "The deleteriousness score from MutationTaster"},
    'PROVEAN_pred':{'Number': '.', 'Type': 'String', 'Description': "The deleteriousness prediction from PROVEAN"},
    'PROVEAN_score':{'Number': '.', 'Type': 'String', 'Description': "The deleteriousness score from PROVEAN"},
    'REVEL_score':{'Number': '.', 'Type': 'String', 'Description': "The REVEL score for missense variants"},
    'ada_score':{'Number': '.', 'Type': 'String', 'Description': "The dbscSNV score using Adaboost algorithm for splicing variant"},
    'rf_score':{'Number': '.', 'Type': 'String', 'Description': "The dbscSNV score using Random Forest algorithm for splicing variant"},
    'CADD_PHRED':{'Number': '.', 'Type': 'String', 'Description': "The CADD score for the variant in PHRED scale"},
    'CADD_RAW':{'Number': '.', 'Type': 'String', 'Description': "The raw CADD score"},
    'Condel':{'Number': '.', 'Type': 'String', 'Description': "The Condel score for the variant"},
    'LoFtool':{'Number': '.', 'Type': 'String', 'Description': "A rank of genic intolerance and consequent susceptibility to disease for the affected gene"},
    'ExACpLI':{'Number': '.', 'Type': 'String', 'Description': "The probabililty of the affected gene being loss-of-function intolerant (pLI)"},
    'HGMD':{'Number': '.', 'Type': 'String', 'Description': "The existing variant in HGMD"},
    'HGMD_CLASS':{'Number': '.', 'Type': 'String', 'Description': "HGMD variatn classifier"},
    'HGMD_MUT':{'Number': '.', 'Type': 'String', 'Description': "HGMD variant type"},
    'HGMD_PHEN':{'Number': '.', 'Type': 'String', 'Description': "The phenotype associated with the existing HGMD variant"},
    'HGMD_RANKSCORE':{'Number': '.', 'Type': 'String', 'Description': "The score for the existing HGMD variant"},
    'GNOMAD_G':{'Number': '.', 'Type': 'String', 'Description': "The existing variant in gnomAD genomes combined population"},
    'GNOMAD_G_AF_AFR':{'Number': '.', 'Type': 'String', 'Description': "Frequency of existing variant in gnomAD genomes African/American population"},
    'GNOMAD_G_AF_AMR':{'Number': '.', 'Type': 'String', 'Description': "Frequency of existing variant in gnomAD genomes American population"},
    'GNOMAD_G_AF_ASJ':{'Number': '.', 'Type': 'String', 'Description': "Frequency of existing variant in gnomAD genomes Ashkenazi Jewish population"},
    'GNOMAD_G_AF_EAS':{'Number': '.', 'Type': 'String', 'Description': "Frequency of existing variant in gnomAD genomes East Asian population"},
    'GNOMAD_G_AF_FIN':{'Number': '.', 'Type': 'String', 'Description': "Frequency of existing variant in gnomAD genomes Finnish population"},
    'GNOMAD_G_AF_NFE':{'Number': '.', 'Type': 'String', 'Description': "Frequency of existing variant in gnomAD genomes Non-Finnish European population"},
    'GNOMAD_G_AF_OTH':{'Number': '.', 'Type': 'String', 'Description': "Frequency of existing variant in gnomAD genomes combined other combined populations"},
    'phastCons100':{'Number': '.', 'Type': 'String', 'Description': "The conservation score for variant site across 100 vertebrates by phastCons method"},
    'phyloP100':{'Number': '.', 'Type': 'String', 'Description': "The conservation score for variant site across 100 vertebrates by phyloP method"},
    'RMSK':{'Number': '.', 'Type': 'String', 'Description': "The identifier by Repeat Masker for repetitive elements overlapping variant site"},
    'GERP':{'Number': '.', 'Type': 'String', 'Description': "The GERP++ conservation score for vatiant site"},
    'COV_GNOMAD_E':{'Number': '.', 'Type': 'String', 'Description': "The percentage of gnomAD exomes covering 20x on the variant site"},
    'COV_GNOMAD_G':{'Number': '.', 'Type': 'String', 'Description': "The percentage of gnomAD genomes covering 20x on the variant site"}
}

def parse_info(col):
    infoH = {}
    tagL = []
    for item in col.split(';'):
        if '=' in item:
            tag = item.split('=')[0]
            val = '='.join(item.split('=')[1:])
            infoH[tag] = val
            tagL.append(tag)
        else:
            infoH[item] = ''
            tagL.append(item)
    #for item
    return(infoH, tagL)

usage = "usage: %prog [options] inFile outFile"
parser = OptionParser(usage)
parser.add_option("--pick", action="store_true", dest="pick_only", help="Pick the most severe consequence only (marked by VEP)")
parser.add_option("--cds", action="store_true", dest="cds_only", help="Keep variants in coding region (CDS) only, i.e., VEP IMPACT is not 'MODIFIER'")
parser.add_option("--vep-gnomad-af", action="store", type="string", dest="vep_gnomad_af", help="VEP field name to be used as 'Variant_frequency_in_gnomAD'")

(options, args) = parser.parse_args()
if len(args) < 2:
    parser.error("Need input & output file names")
    sys.exit(1)
pick_only = False
if options.pick_only:
    pick_only = True
cds_only = False
if options.cds_only:
    cds_only = True
GNOMAD_AF_FIELD='gnomAD_AF'
if options.vep_gnomad_af:
    out_columnH[ options.vep_gnomad_af ] = out_columnH[ 'gnomAD_AF' ]
    out_columnH.pop('gnomAD_AF', None)
    GNOMAD_AF_FIELD = options.vep_gnomad_af

inFile = gzip.open(args[0], 'rb')
outFile = gzip.open(args[1], 'wb')
for line in inFile:
    if line[:2] == '##':
        if 'INFO=<ID=CSQ' in line:
            csq_header = re.match('##.*Format: (.*)">',line).groups()[0]
            csq_headerL = csq_header.split('|')
        if 'INFO=<ID=' in line:
            continue
        else:
            outFile.write(line)
    elif line[:6] == '#CHROM':
        for h in out_columnH:
            id=h
            if h == 'AF':
                id='1000G_AF'
            outFile.write("##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\">\n" % (out_columnH[h]['Name'], out_columnH[h]['Number'], out_columnH[h]['Type'], out_columnH[h]['Description']))
        #for h
        ## additional field
        outFile.write("##INFO=<ID=Variant_frequency_as_text,Number=1,Type=String,Description=\"The variant allele frequency in gnomAD combined population as discrete text categories. Possible values: Novel, Rare (variant frequency less than 1%), Common (variant frequency greater than or equal to 1%).\">\n")
        outFile.write(line)
    else:
        colL=line.strip().split('\t')
        (infoH, tagL)=parse_info(colL[7])
        if 'CSQ' not in infoH:
            continue

        tagL.pop()
        csqL = infoH['CSQ'].split(',')
        annotL = []
        for csq in csqL:
            csq_valL = csq.split('|')
            if cds_only and csq_valL[ csq_headerL.index('IMPACT') ] == 'MODIFIER':
                continue
            if pick_only and csq_valL[ csq_headerL.index('PICK') ] == '':
                continue
            outText = ''
            for h in out_columnH:
                csq_val = csq_valL[ csq_headerL.index(h) ]
                if h == 'Consequence':
                    csq_val = csq_val.split('&')[0]
                if h == GNOMAD_AF_FIELD and (csq_val == '.' or csq_val == ''):
                    csq_val = -10
                if csq_val != '' and csq_val != '.':
                    outText += "%s=%s;" % (out_columnH[h]['Name'], csq_val)
                if h == GNOMAD_AF_FIELD:
                    if float(csq_val) < 0:
                        outText += "Variant_frequency_as_text=Novel;"
                    elif float(csq_val) < 0.01:
                        outText += "Variant_frequency_as_text=Rare;"
                    else:
                        outText += "Variant_frequency_as_text=Common;"
            annotL.append(outText)
        #for csq
        for annot in set(annotL):
            outFile.write('\t'.join(colL[:7]))
            outFile.write('\t%s' % annot) 
            outFile.write('\t%s\n' % '\t'.join(colL[8:]))
    #else
#for line
outFile.flush()
outFile.close()
inFile.close()
