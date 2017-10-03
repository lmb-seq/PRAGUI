#!/usr/bin/python

import gzip
import multiprocessing
import os
import random
import string
import subprocess
import sys
import uuid
import glob

import cross_fil_util as util
import numpy as np

sys.path.append('/home/paulafp/Documents/temp/crossFil/')
from readCsvFile import readCsvFile


PROG_NAME = 'RNAseq Pipeline'
DESCRIPTION = 'Fastq files to RNAseq data analysis.'

QUIET   = False
LOGGING = False

FILE_TAG = '_rnapip_' # This tag is used for formatting file names so they can be passed between the various RNAseq pipeline programs 

TEMP_ID = '%s' % uuid.uuid4()
LOG_FILE_PATH = 'rnapip-out-%s.log' % TEMP_ID
LOG_FILE_OBJ = None # Created when needed
MAX_CORES = multiprocessing.cpu_count()


ALIGNERS = ('STAR', 'hisat2', 'tophat2')
ALIGNER_STAR, ALIGNER_HISAT2, ALIGNER_TOPHAT2 = ALIGNERS
DEFAULT_ALIGNER = ALIGNER_STAR
OTHER_ALIGNERS = [ALIGNER_HISAT2, ALIGNER_TOPHAT2]


def exists_skip(filename):
  if os.path.exists(filename):
    print('%s already exists and will not be overwritten. Skipping this file...' % filename)
    return(False)
  else:
    return(True)

def append_to_file_name(file_name,extension):
  new_file_name = file_name + extension
  return(new_file_name)
  
def rm_low_mapq(in_file,out_file_name,mapq):
  cmdArgs=['samtools', 'view', '-bq',
          str(mapq), in_file]
  out_file = open(out_file_name,'wb')
  util.call(cmdArgs,stdout=out_file)
  out_file.close()

if __name__ == '__main__':

  from argparse import ArgumentParser
   
  epilog = 'For further help on running this program please email paulafp@mrc-lmb.cam.ac.uk.\n\n'
  epilog += 'Example use:\n\n'
  
  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)
  
#  arg_parse.add_argument('fastq_paths', nargs='+', metavar='FASTQ_FILES',
#                         help='File paths of FASTQ sequence read files (may contain wildcards) ') 
  
  arg_parse.add_argument('samples_csv', metavar='SAMPLES_CSV',
                         help='File path of a comma-separated file containing the samples names, the file path for read1, the file path for read2 and the experimental condition (e.g. Mutant or Wild-type). For single-ended experiments, please fill read2 slot with NA.') 
  
  arg_parse.add_argument('genome_fasta', metavar='GENOME_FASTA',
                         help='File path of genome sequence FASTA file (for use by genome aligner)') 
  
  arg_parse.add_argument('genome_gtf', metavar='GENE_ANNOTATIONS_GTF',
                         help='File path of gene annotations in gtf/gff format (for use by htseq-count)') 
  
  arg_parse.add_argument('-trim_galore', # metavar='TRIM_GALORE_OPTIONS',
                         default=None, 
                         help='options to be provided to trim_galore. They should be provided under quotes. If not provided, trim_galore will run with developer\'s default options.')
  
  arg_parse.add_argument('-fastqc_args', metavar='FASTQC', 
                         default=None,
                         help='options to be provided to fastqc. They should be provided under double quotes. If not provided, fastqc will run with developer\'s default options.')
                         
  arg_parse.add_argument('-skipfastqc', default=False, action='store_true',
                         help='Option to skip fastqc step. If this option is set, the option -fastqc_args will be ignored.')
  
  arg_parse.add_argument('-al', metavar='ALIGNER_NAME', default=DEFAULT_ALIGNER,
                         help='Name of the program to perform the genome alignment/mapping: Default: %s Other options: %s' % (DEFAULT_ALIGNER, OTHER_ALIGNERS)) 
  
  arg_parse.add_argument('-star_index', metavar='STAR_GENOME_INDEX', default=None,
                         help='Path to directory where genome indices are stored.') 
  
  arg_parse.add_argument('-mapq', default=20, type=int,
                         help='Threshold below which reads will be removed from the aligned bam file.') 
  
  arg_parse.add_argument('-barcode_csv',  metavar='BARCODE_CSV_FILE',
                         help='CSV format file containing barcode strain/sample names') 
  
  arg_parse.add_argument('-q', default=False, action='store_true',
                         help='Sets quiet mode to supress on-screen reporting.')
  
  arg_parse.add_argument('-log', default=False, action='store_true',
                         help='Log all reported output to a file.')
  
  arg_parse.add_argument('-cpu', metavar='NUM_CORES', default=util.MAX_CORES, type=int,
                         help='Number of parallel CPU cores to use. Default: All available (%d)' % util.MAX_CORES) 
 
  arg_parse.add_argument('-pe', nargs=2, metavar='PAIRED_READ_TAGS', default=['r_1','r_2'],
                         help='The subtrings/tags which are the only differences between paired FASTQ file paths. Default: r_1 r_2') 
  
  arg_parse.add_argument('-se', default=False, action='store_true',
                         help='Input reads are single-end data, otherwise defaults to paired-end.')
  
  arg_parse.add_argument('-stranded', default=False, action='store_true',
                         help='Input strand-specific protocol, otherwise defaults to non-strand-specific protocol.')
  
  args = vars(arg_parse.parse_args())

  # fastq_paths   = args['fastq_paths']
  samples_csv   = args['samples_csv']
  genome_fasta  = args['genome_fasta']
  genome_gtf    = args['genome_gtf']
  trim_galore   = args['trim_galore']
  skipfastqc    = args['skipfastqc']
  fastqc_args   = args['fastqc_args']
  aligner       = args['al']
  star_index    = args['star_index']
  mapq          = args['mapq']
  barcode_csv   = args['barcode_csv']
  pair_tags     = args['pe']
  is_single_end = args['se']
  # out_top_dir   = args['outdir']
  num_cpu       = args['cpu'] or None # May not be zero
  
  
  # Parse input comma separated file
  
  csv = readCsvFile(filename=samples_csv,separator='\t',header=False) # returns numpy array
  
  if is_single_end:
    fastq_paths = list(csv[:,1])
  else:
    fastq_paths = list(csv[:,1]) + list(csv[:,2])
  
  # Run Trim_galore followed by fastqc
  
  fastq_paths2 = []
  fastq_paths3 = []
  
  for f in fastq_paths:
    f = os.path.expanduser(f)
    f0 = f
    f=f.split(".")
    if f[-1] == 'gz':
      f = f[:-2]
    else:
      f = f[:-1]
    f = '.'.join(f)
    trimmed_filename = f+'_trimmed.fq.gz'
    if exists_skip(trimmed_filename):
      fastq_paths2.append(f0)
    fastq_paths3.append(f)
  
  
  if fastq_paths2 != []:
    cmdArgs = ['trim_galore','--gzip']
    
    if trim_galore is not None:
      trim_galore = trim_galore.split(' ')
      cmdArgs += trim_galore
    
    if '-o' not in cmdArgs or '--output_dir' not in cmdArgs:
       cmdArgs.append('-o')
       od = f.split("/")
       od = od[:-1]
       od = '/'.join(od)
       cmdArgs.append(od)
    
    
    if skipfastqc is False:
      cmdArgs += ['-fastqc']
      if fastqc_args is not None:
        cmdArgs += ['-fastqc_args',fastqc_args]
    else:
      print('Skipping fastqc step...\n')
    
    cmdArgs += fastq_paths2
    
    util.call(cmdArgs)
  
  # Run Aligner
  
  # Check whether genomes indices are present. If not, create them.
  
  # ADD STAR OPTIONS!!!
  if aligner is ALIGNER_STAR:
    if not os.path.exists(star_index):
      print('STAR indices not found. Generating STAR indices...\n')
      os.mkdir(star_index)
      cmdArgs = [ALIGNER_STAR,
                 '--runMode','genomeGenerate',
                 '--genomeDir',star_index,
                 '--genomeFastaFiles', genome_fasta,
                 '--runThreadN',str(num_cpu)]
      util.call(cmdArgs)
    
    print('\nAligning reads using STAR...\n')
    
    cmdArgs = (ALIGNER_STAR,
               '--genomeDir',star_index,
               '--runThreadN',str(num_cpu),
               '--readFilesCommand', 'zcat', '-c',
               '--outSAMtype','BAM','SortedByCoordinate',
               '--readFilesIn')
    
    trimmed_fq = [append_to_file_name(x,'_trimmed.fq.gz') for x in fastq_paths3]
    bam_files = []
    
    if is_single_end:
      
      print("Running single-end mode...\n")
      
      #trimmed_fq = glob.glob('*trimmed.fq.gz')
      #trimmed_fq = ','.join(trimmed_fq)
      
      for f in trimmed_fq:
        if mapq > 0 :
          bam = '%s.sorted_fil_%d.out.bam' % (f,mapq)
        else:
          bam = '%s.sorted.out.bam' % f
        
        bam_files.append(bam)
        if exists_skip(bam):
          cmdArgs_se = list(cmdArgs)
          cmdArgs_se.append(f)
          util.call(cmdArgs_se)
          if mapq > 0 :
            rm_low_mapq('./Aligned.sortedByCoord.out.bam',bam,mapq) # Remove reads with quality below mapq
            os.remove('./Aligned.sortedByCoord.out.bam')
          else:
            os.rename('./Aligned.sortedByCoord.out.bam',bam)
  
    else:
      
      print("Running paired-end mode...\n")
      
      trimmed_fq_r1 = '*%s*trimmed.fq.gz' % pair_tags[0]
      trimmed_fq_r1 = list(filter(lambda x:trimmed_fq_r1 in x, trimmed_fq)) # grep for python3
      #trimmed_fq_r1 = glob.glob(trimmed_fq_r1)
      
      trimmed_fq_r2 = '*%s*trimmed.fq.gz' % pair_tags[1]
      trimmed_fq_r2 = list(filter(lambda x:trimmed_fq_r2 in x, trimmed_fq)) # grep for python3
      #trimmed_fq_r2 = glob.glob(trimmed_fq_r2)
      
      if len(trimmed_fq_r1) != len(trimmed_fq_r2):
        sys.exit('ERROR: Number of fq files differs for read1 and read2... Exiting...\n')
      

      for i in range(0,len(trimmed_fq_r1)):
        if mapq > 0 :
          bam = './%s.pe.sorted_fil_%d.out.bam' % (trimmed_fq_r1[i],mapq)
        else:
          bam = './%s.pe.sorted.out.bam' % trimmed_fq_r1[i]
        
        
        bam_files = bam_files.append(bam)
        if exists_skip(bam):
          cmdArgs_pe = list(cmdArgs)
          cmdArgs_pe += [trimmed_fq_r1[i],trimmed_fq_r2[i]]
          util.call(cmdArgs_pe)
          if mapq > 0 :
            rm_low_mapq('./Aligned.sortedByCoord.out.bam',bam,mapq) # Remove reads with quality below mapq
            os.remove('./Aligned.sortedByCoord.out.bam')
          else:
            os.rename('./Aligned.sortedByCoord.out.bam',bam)
        
  
  # Generate Count matrix with HTSeq
  
  rc_file_list = []
  
  for f in bam_files:
    rc_file = '%s_count_table.txt' % f
    rc_file_list.append(rc_file)
    if exists_skip(rc_file):
      fileObj = open(rc_file,'wb')
      cmdArgs = ['htseq-count','--format=bam','--stranded=no']
      cmdArgs += [f,genome_gtf]
      util.call(cmdArgs,stdout=fileObj)
      fileObj.close()
  
  
  # Create csv file for DESeq function DESeqDataSetFromHTSeqCount
  
  csv_deseq_name = samples_csv + '_DESeq_table.txt'
  
  if exists_skip(csv_deseq_name):
    
    N = csv.shape[0]
    
    csv_deseq = np.zeros((N,3))
    csv_deseq = np.array(csv_deseq,dtype=object) # dtype=object provides an array of python object references. 
                                                 # It can have all the behaviours of python strings.
    
    csv_deseq[:,0] = csv[:,0]
    csv_deseq[:,1] = np.array(rc_file_list)
    csv_deseq[:,2] = csv[:,-1]
    
    np.savetxt(fname=csv_deseq_name,X=csv_deseq,delimiter='\t',fmt='%s')
  
  
  # Gene Expression analysis using R
  
  exploratory_analysis_plots = samples_csv + '_sclust.pdf'
  TPMs = samples_csv + '_tpm.txt'
  DESeq_summary = samples_csv + '_DESeq_summary.txt'
  DESeq_results = samples_csv + '_DESeq_results.txt'
  
  i=[]
  
  if exists_skip(exploratory_analysis_plots):  # Gene expression analysis has 3 steps.
                                               # These do not need to be repeated if they have
    i.append("ea")                             # already been run. Therefore, the script checks
                                               # whether the output files have been generated
  if exists_skip(TPMs):                        # and stores a specific flag each time that's the case.
                                               # The following R script checks which flags have been
    i.append("tpm")                            # stored and thus knows which steps to skip (if any).
    
  if exists_skip(DESeq_results):
    
    i.append("deseq")
    
  print(i)
  print(len(i))
  
  if len(i) > 0:
    i = "_".join(i)

    DESeq_out_obj = open(DESeq_summary,"wb")
    
    cmdArgs = ['Rscript','--vanilla', os.environ["RNAseq_analysis"], csv_deseq_name, i]
    
    util.call(cmdArgs,stdout=DESeq_out_obj)
    
    DESeq_out_obj.close()

    
    
  
  
  
  
