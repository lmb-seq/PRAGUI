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
import numpy as np
import shutil
import HTSeq

current_path = os.path.realpath(__file__)
current_path = os.path.dirname(current_path) + '/cell_bio_util'

sys.path.append(current_path)
import cell_bio_util as util

from readCsvFile import readCsvFile


PROG_NAME = 'RNAseq Pipeline'
DESCRIPTION = 'Process fastq files to RNAseq data analysis.'

util.init_app('rnapip') # Redefine variables from cross_fil_util.py

ALIGNERS = ('STAR', 'hisat2', 'tophat2')
ALIGNER_STAR, ALIGNER_HISAT2, ALIGNER_TOPHAT2 = ALIGNERS
DEFAULT_ALIGNER = ALIGNER_STAR
OTHER_ALIGNERS = [ALIGNER_HISAT2, ALIGNER_TOPHAT2]


def exists_skip(filename):
  if os.path.exists(filename):
    util.info('%s already exists and will not be overwritten. Skipping this folder/file...' % filename)
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
  
def new_dir(new_dir):
  if not exists_skip(new_dir):
    new_dir = util.get_temp_path(new_dir)
    util.info('Output from Cuffdiff will be saved in %s' % new_dir)
  os.makedirs(new_dir,exist_ok = True) #,mode = 0o666)
  return(new_dir)
  
def report_cuff_version(CUFF_PROG):
  cufflinks_vs_obj = open('cufflinks_version_control.txt','a')
  util.call(CUFF_PROG,stderr=cufflinks_vs_obj)
  cufflinks_vs_obj.close()
  cufflinks_vs_obj = open('cufflinks_version_control.txt','r')
  cufflinks_vs = cufflinks_vs_obj.readline()
  cufflinks_vs_obj.close()
  util.LOG_FILE_OBJ.write(cufflinks_vs)
  os.remove('cufflinks_version_control.txt')

def rm_lines(f_in,f_out,string = '> Processing Locus'):
  f_in_obj  = open(f_in,'r')
  f_out_obj = open(f_out,'a')
  line0 = None
  for line in f_in_obj:
    if not string in line:
      if line0 is not None:
        f_out_obj.write(line0)
        line0 = None
      f_out_obj.write(line)
    else:
      line0 = line
  f_in_obj.close()
  f_out_obj.flush()
  os.remove(f_in)


if __name__ == '__main__':

  from argparse import ArgumentParser
   
  epilog = 'For further help on running this program please email paulafp@mrc-lmb.cam.ac.uk.\n\n'
  epilog += 'Example use:\n\n'
  
  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)
  
  arg_parse.add_argument('samples_csv', metavar='SAMPLES_CSV',
                         help='File path of a tab-separated file containing the samples names, the file path for read1, the file path for read2, the experimental condition (e.g. Mutant or Wild-type) and any other information to be used as contrasts for differential expression calling. For single-ended experiments, please fill read2 slot with NA.') 
  
  arg_parse.add_argument('genome_fasta', metavar='GENOME_FASTA',
                         help='File path of genome sequence FASTA file (for use by genome aligner)') 

  arg_parse.add_argument('-analysis_type', metavar='ANALYSIS_TYPE',default=['DESeq','Cufflinks'][0],
                         help='Specify whether to perform analysis using DESeq2 or Cufflinks. Default is set to DESeq2.') 
  
  arg_parse.add_argument('-genome_gtf', metavar='GENOME_ANNOTATIONS_GTF', default=None,
                         help='File path of gene annotations in gtf/gff format (for use by htseq-count). This file is only required when performing an analysis using DESeq2.') 
  
  arg_parse.add_argument('-geneset_gtf', default=None,
                         help='File path of gene annotations in gtf/gff format needed to compute TPMs. If this file is not provided, GENOME_ANNOTATIONS_GTF will be used.') 
  
  arg_parse.add_argument('-trim_galore', # metavar='TRIM_GALORE_OPTIONS',
                         default=None, 
                        help='Options to be provided to trim_galore. They should be provided under quotes. If not provided, trim_galore will run with developer\'s default options.')
  
  arg_parse.add_argument('-fastqc_args', metavar='FASTQC', 
                         default=None,
                         help='options to be provided to fastqc. They should be provided under double quotes. If not provided, fastqc will run with developer\'s default options.')
                         
  arg_parse.add_argument('-skipfastqc', default=False, action='store_true',
                         help='Option to skip fastqc step. If this option is set, the option -fastqc_args will be ignored.')
  
  arg_parse.add_argument('-al', metavar='ALIGNER_NAME', default=DEFAULT_ALIGNER,
                         help='Name of the program to perform the genome alignment/mapping: Default: %s Other options: %s' % (DEFAULT_ALIGNER, OTHER_ALIGNERS)) 
  
  arg_parse.add_argument('-star_index', metavar='STAR_GENOME_INDEX', default=None,
                         help='Path to directory where genome indices are stored.') 

  arg_parse.add_argument('-star_args', default=None,
                         help='Options to be provided to STAR. They should be provided under double quotes. If not provided, STAR will be expecting the following options: --readFilesCommand zcat -c, --outSAMtype BAM, SortedByCoordinate')   
  
  arg_parse.add_argument('-mapq', default=20, type=int,
                         help='Threshold below which reads will be removed from the aligned bam file.') 
  
  arg_parse.add_argument('-cpu', metavar='NUM_CORES', default=util.MAX_CORES, type=int,
                         help='Number of parallel CPU cores to use. Default: All available (%d)' % util.MAX_CORES) 
 
  arg_parse.add_argument('-pe', nargs=2, metavar='PAIRED_READ_TAGS', default=['r_1','r_2'],
                        help='The subtrings/tags which are the only differences between paired FASTQ file paths. Default: r_1 r_2') 
  
  arg_parse.add_argument('-se', default=False, action='store_true',
                         help='Input reads are single-end data, otherwise defaults to paired-end.')
  
  arg_parse.add_argument('-stranded', default=False, action='store_true',
                         help='Input strand-specific protocol, otherwise defaults to non-strand-specific protocol.')
  
  arg_parse.add_argument('-contrast', default='condition',
                         help='Set column from SAMPLES_CSV file to be used as contrast by DESeq2 otherwise defaults to the third column')
  
  arg_parse.add_argument('-contrast_levels', nargs=2, default=None,
                         help='Set comparisons for DESeq2. By default, DESeq2 compare last level over the first level from the CONTRAST column.')

  arg_parse.add_argument('-cuff_opt', default=None,
                         help='options to be provided to cufflinks. They should be provided under quotes. If not provided, cufflinks will run with developer\'s default options.')

  arg_parse.add_argument('-cuff_gtf', default=False, action='store_true',
                         help='Set "-g" option from cufflinks and use file specified in "-geneset_gtf" option. This option should not be set if "-cuff_gtf" already incorporates a gtf file to be used.')
  
  arg_parse.add_argument('-cuffnorm', default=False, action='store_true',
                         help='Specify whether Cuffnorm should be executed besides Cuffdiff.')
  
  arg_parse.add_argument('-q', default=False, action='store_true',
                         help='Sets quiet mode to supress on-screen reporting.')
  
  arg_parse.add_argument('-log', default=False, action='store_true',
                         help='Log all reported output to a file.')

  args = vars(arg_parse.parse_args())

  samples_csv   = args['samples_csv']
  genome_fasta  = args['genome_fasta']
  analysis_type = args['analysis_type']
  genome_gtf    = args['genome_gtf']
  geneset_gtf   = args['geneset_gtf']
  trim_galore   = args['trim_galore']
  skipfastqc    = args['skipfastqc']
  fastqc_args   = args['fastqc_args']
  aligner       = args['al']
  star_index    = args['star_index']
  star_args     = args['star_args']
  mapq          = args['mapq']
  num_cpu       = args['cpu'] or None # May not be zero
  pair_tags     = args['pe']
  is_single_end = args['se']
  stranded      = args['stranded']
  contrast      = args['contrast']
  levels        = args['contrast_levels']
  cuff_opt      = args['cuff_opt']
  cuff_gtf      = args['cuff_gtf']
  cuffnorm      = args['cuffnorm']
  
  # Reporting handled by cross_fil_util.py (submodule)
  util.QUIET   = args['q']
  util.LOGGING = args['log']
  
  python_version = 'python version ' + sys.version + '\n'
  
  python_command = ' '.join(sys.argv) + '\n'
  
  util.info(python_version)
  util.info(python_command)
  
  if analysis_type == 'DESeq':
    util.info('Differential gene expression analysis using DESeq2...')
    if genome_gtf is None:
      util.critical('Expecting file with gene annotations in gtf/gff format. Please provide full file path using the "-genome_gtf" option...')
  elif analysis_type == 'Cufflinks':
    util.info('Analysis of transcript expression using Cufflinks...')
  else:
    util.critical('Expecting ANALYSIS_TYPE to be either DESeq2 or Cufflinks...')
  
  if geneset_gtf is None:
    geneset_gtf = genome_gtf
  
  
  # Parse input comma separated file
  
  csvfile = open(samples_csv,'r')                  # Get header from csv file and build new header for
  header = csvfile.readline()                      # input table needed for analysis in R.
  csvfile.close()
  header = header.split()
  header2 = header
  header = ['samplename','filename'] + header[3:]
  header = np.array(header)
  
  csv = readCsvFile(filename=samples_csv,separator='\t',header=True) # returns numpy array
  
  cmdArgs = ['trim_galore','--gzip']
  
  fastq_paths2 = []
  trimmed_fq = []
  fastq_dirs  = []
  
  if trim_galore is not None:
    trim_galore = trim_galore.split(' ')
    cmdArgs += trim_galore
  
  # Output from trim_galore will be saved in a new folder './trim_galore' 
  # if not otherwise specified in the command line
  if '-o' in cmdArgs:
    ind = cmdArgs.index('-o') + 1
    od  = cmdArgs[ind]  
  elif '--output_dir' in cmdArgs:
    ind = cmdArgs.index('--output_dir') + 1
    od  = cmdArgs[ind]    
  else:
    cmdArgs.append('-o')
    if exists_skip('./trim_galore'):
      os.makedirs('./trim_galore',exist_ok = True,mode = 0o666)
    od = './trim_galore'
    cmdArgs.append(od)
  
  if skipfastqc is False:
    cmdArgs += ['-fastqc']
    if fastqc_args is not None:
      cmdArgs += ['-fastqc_args',fastqc_args]
  else:
    util.info('Skipping fastqc step...')
  
  
  if is_single_end:
    util.info('User specified input data to be single-end... Running single-end mode...')
    fastq_paths = list(csv[:,1])
    
    for f in fastq_paths:
      f0 = os.path.expanduser(f)
      d = os.path.dirname(f0)
      f = os.path.basename(f)
      f=f.split(".")
      if f[-1] == 'gz':
        f = f[:-2]
      else:
        f = f[:-1]
      f = '.'.join(f)
      trimmed_filename = od + '/' + f +'_trimmed.fq.gz'
      if exists_skip(trimmed_filename):
        fastq_paths2.append(f0)
      trimmed_fq.append(trimmed_filename)
      fastq_dirs.append(d)
    
  else:
    util.info('User specified input data to be paired-end... Running paired-end mode with tags %s and %s...' % (pair_tags[0],pair_tags[1]) )
    cmdArgs.append('--paired')
    
    fastq_paths = []
    R = csv.shape[0]
    
    for i in range(R):
      for j in [1,2]:
        fastq_paths.append(csv[i,j])
      
    for f in fastq_paths:
      f0 = os.path.expanduser(f)
      f = os.path.basename(f)
      f=f.split(".")
      if f[-1] == 'gz':
        f = f[:-2]
      else:
        f = f[:-1]
      f = '.'.join(f)
      if pair_tags[0] in f:
        trimmed_filename = od + '/' + f + '_val_1.fq.gz'
        d = os.path.dirname(f0)             # directory where fastq file is stored
      elif pair_tags[1] in f:
        trimmed_filename = od + '/' + f + '_val_2.fq.gz'
      else:
        util.critical('Paired read tag not found... Exiting...')
      
      if exists_skip(trimmed_filename):
        fastq_paths2.append(f0)
      #fastq_paths3.append(f)
      trimmed_fq.append(trimmed_filename)
      fastq_dirs.append(d)
  
  # Run Trim_galore followed by fastqc
  
  if fastq_paths2 != []:
    
    #util.call(['trim_galore','-v'],stdout=util.LOG_FILE_OBJ) No need because Trim_galore reports version automatically
    
    if skipfastqc is False:
      util.call(['fastqc','-v'],stdout=util.LOG_FILE_OBJ)
    
    cmdArgs += fastq_paths2
    
    util.call(cmdArgs)
    
  
  # Run Aligner
  
  # Check whether genomes indices are present. If not, create them.
  
  if aligner is ALIGNER_STAR:
    if star_index is None:
      star_index = os.path.dirname(genome_fasta) + '/star-genome/'
      util.warn('Folder where STAR indices are located hasn\'t been specified. Program will default to %s...' % star_index)
    if not os.path.exists(star_index):
      util.info('STAR indices not found. Generating STAR indices...')
      os.mkdir(star_index)
      cmdArgs = [ALIGNER_STAR,
                 '--runMode','genomeGenerate',
                 '--genomeDir',star_index,
                 '--genomeFastaFiles', genome_fasta,
                 '--runThreadN',str(num_cpu)]
      util.call([ALIGNER_STAR,'--version'],stdout=util.LOG_FILE_OBJ)
      util.call(cmdArgs)
    
    util.info('Aligning reads using STAR...')
    
    cmdArgs = [ALIGNER_STAR,
               '--genomeDir',star_index,
               '--runThreadN',str(num_cpu)]
    
    if star_args is None:
      cmdArgs += ['--readFilesCommand', 'zcat', '-c',
                  '--outSAMtype','BAM','SortedByCoordinate',
                  '--readFilesIn']
    else:
      star_args = star_args.split()
      cmdArgs += star_args
      cmdArgs.append('--readFilesIn')
    
    bam_files = []
    
    k=0
    
    if is_single_end:
      
      util.info('Running single-end mode...')

      for f in trimmed_fq:
        fo = os.path.basename(f)
        fo = fastq_dirs[k]+ '/' + fo
        if mapq > 0 :
          bam = '%s.sorted_fil_%d.out.bam' % (fo,mapq)
        else:
          bam = '%s.sorted.out.bam' % fo
        bam_files.append(bam)
        
        if exists_skip(bam):
          cmdArgs_se = list(cmdArgs)
          cmdArgs_se.append(f)
          util.call([ALIGNER_STAR,'--version'],stdout=util.LOG_FILE_OBJ)
          util.call(cmdArgs_se)
          if mapq > 0 :
            util.call(['samtools','--version'],stdout=util.LOG_FILE_OBJ)
            rm_low_mapq('./Aligned.sortedByCoord.out.bam',bam,mapq) # Remove reads with quality below mapq
            os.remove('./Aligned.sortedByCoord.out.bam')
          else:
            os.rename('./Aligned.sortedByCoord.out.bam',bam)
        k+=1
  
    else:
      
      util.info('Running paired-end mode...')
      
      trimmed_fq_r1 = list(filter(lambda x:pair_tags[0] in x, trimmed_fq)) # grep for python3
      trimmed_fq_r2 = list(filter(lambda x:pair_tags[1] in x, trimmed_fq))
      
      if len(trimmed_fq_r1) != len(trimmed_fq_r2):
        util.critical('Number of fq files differs for read1 and read2... Exiting...')
      
      for i in range(0,len(trimmed_fq_r1)):
        fo = os.path.basename(trimmed_fq_r1[i])
        fo = fastq_dirs[k] + '/' + fo
        if mapq > 0 :
          bam = '%s.pe.sorted_fil_%d.out.bam' % (fo,mapq)
        else:
          bam = '%s.pe.sorted.out.bam' % fo
        
        bam_files.append(bam)
        
        if exists_skip(bam):
          cmdArgs_pe = list(cmdArgs)
          cmdArgs_pe += [trimmed_fq_r1[i],trimmed_fq_r2[i]]
          util.call([ALIGNER_STAR,'--version'],stdout=util.LOG_FILE_OBJ)
          util.call(cmdArgs_pe)
          if mapq > 0 :
            util.call(['samtools','--version'],stdout=util.LOG_FILE_OBJ)
            rm_low_mapq('./Aligned.sortedByCoord.out.bam',bam,mapq) # Remove reads with quality below mapq
            os.remove('./Aligned.sortedByCoord.out.bam')
          else:
            os.rename('./Aligned.sortedByCoord.out.bam',bam)  
        k+=1
  
  #########################
  ## Analysis with DESeq2##
  #########################
  
  if analysis_type is 'DESeq':

    # Generate Count matrix with HTSeq
    
    rc_file_list = []
    
    for f in bam_files:
      rc_file = '%s_count_table.txt' % f
      rc_file_list.append(rc_file)
      if exists_skip(rc_file):
        htseq_version = HTSeq.__version__
        util.info('HTSeq version %s' % htseq_version)
        fileObj = open(rc_file,'wb')
        cmdArgs = ['htseq-count','--format=bam','--stranded=no']
        cmdArgs += [f,genome_gtf]
        util.call(cmdArgs,stdout=fileObj)
        fileObj.close()
    
    
    # Create csv file for DESeq function DESeqDataSetFromHTSeqCount
    
    deseq_dir = rc_file_list[0].split('/')
    deseq_dir = deseq_dir[:-1]
    deseq_dir = '/'.join(deseq_dir) + '/'
    
    deseq_head = samples_csv.split('/')
    deseq_head = deseq_head[-1]
    deseq_head = deseq_dir + deseq_head
    
    csv_deseq_name = append_to_file_name(deseq_head,'_DESeq_table.txt')
    
    if exists_skip(csv_deseq_name):
      
      M = csv.shape[0]
      N = csv.shape[1] - 1
      
      csv_deseq = np.zeros((M,N))
      csv_deseq = np.array(csv_deseq,dtype=object) # dtype=object provides an array of python object references. 
                                                   # It can have all the behaviours of python strings.
      
      csv_deseq[:,0] = csv[:,0]
      csv_deseq[:,1] = np.array(rc_file_list)
      csv_deseq[:,2:] = csv[:,3:]
      
      csv_deseq_wh = np.zeros((M+1,N))
      csv_deseq_wh = np.array(csv_deseq_wh,dtype=object)
      csv_deseq_wh[0,:] = header
      csv_deseq_wh[1:,:] = csv_deseq
      
      np.savetxt(fname=csv_deseq_name,X=csv_deseq_wh,delimiter='\t',fmt='%s')
    
    
    # Gene Expression analysis using R
    
    exploratory_analysis_plots = append_to_file_name(deseq_head, '_sclust.pdf')
    TPMs = append_to_file_name(deseq_head,'_tpm.txt')
    DESeq_summary = append_to_file_name(deseq_head,'_DESeq_summary.txt')
    DESeq_results = append_to_file_name(deseq_head,'_DESeq_results.txt')
    
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
      
    
    if len(i) > 0:
      i = "_".join(i)
      
      if levels is None:
        cmdArgs = ['Rscript','--vanilla', os.environ["RNAseq_analysis"], csv_deseq_name, i, geneset_gtf, contrast]
      else:
        cmdArgs = ['Rscript','--vanilla', os.environ["RNAseq_analysis"], csv_deseq_name, i, geneset_gtf, contrast] + levels
      
      if "deseq" in i:
        DESeq_out_obj = open(DESeq_summary,"wb")
        util.call(cmdArgs,stdout=DESeq_out_obj)
        DESeq_out_obj.close()
      else:
        util.call(cmdArgs)
      
      util.logging('')
      sessionInfo_file = deseq_head + '_sessionInfo.txt'
      shutil.copyfileobj(open(sessionInfo_file, 'r'), util.LOG_FILE_OBJ)
      os.remove(sessionInfo_file)
  
  
  #############################
  ## Analysis with Cufflinks ##
  #############################
  
  if analysis_type == 'Cufflinks':
    util.info('Running Cufflinks...')
    
    out_folder = './'
    library_type = None
    no_output_folder = True
    
    # Get cufflinks options
    if cuff_opt is not None:
      cuff_opt = cuff_opt.split(' ')
      # Get output folder if specified as argument
      if '-o' in cuff_opt:
        ind = cuff_opt.index('-o') + 1
        out_folder = cuff_opt[ind] +'/'
        util.info('Output folder for Cufflinks has been specified. Saved all output in:%s' % out_folder)
        no_output_folder = False
      if '--library-type' in cuff_opt:
        ind2 = cuff_opt.index('--library-type') + 1
        library_type = ['--library-type',cuff_opt[ind2]]
      is_gtf_specified = '-g' in cuff_opt or '–GTF-guide' in cuff_opt
      if '-g' in cuff_opt:
        ind3 = cuff_opt.index('-g')+1
        cuff_gtf_file = ['-g',cuff_opt[ind3]]
      if '-GTF-guide' in cuff_opt:
        ind3 = cuff_opt.index('-GTF-guide')+1
        cuff_gtf_file = ['-g',cuff_opt[ind3]]
    
    # Create assemblies file needed for cuffmerge    
    assemblies = out_folder + 'assembly_GTF_list.txt'
    if os.path.exists(assemblies):
      os.remove(assemblies) 
    fileObj_assemblies = open(assemblies,'a')
    
    # Index bam files using samtools
    for f in bam_files:
      fi = f + '.bai'
      if exists_skip(fi):
        util.info('Indexing file %s...' % f)
        util.call(['samtools','--version'],stdout=util.LOG_FILE_OBJ)
        cmdArgs = ['samtools','index',f]
        util.call(cmdArgs)
      
    # Run Cufflinks command 
      cuff_files = ['genes.fpkm_tracking', 'isoforms.fpkm_tracking', 'skipped.gtf', 'transcripts.gtf']
      if no_output_folder:
        header_cuff = f + '_'
      else:
        f2 = f.split('/')[-1]
        header_cuff = out_folder + f2 + '_'
      f_transcripts = header_cuff + cuff_files[3]
      
      fileObj_assemblies.write(f_transcripts + '\n')
      
      if exists_skip(f_transcripts):
        
        report_cuff_version('cufflinks') # Report version of cufflinks
        
        cmdArgs = ['cufflinks','-p',str(num_cpu)]
        if cuff_opt is not None:
          cmdArgs += cuff_opt
        else:
          util.warn('No options were specified for Cufflinks. Developer\'s default options will be used...')
        if no_output_folder:
          util.info('No output folder for cufflinks has been specified. Files will be saved in the same folder as %s...' % f)
        if cuff_gtf is True:
          if not is_gtf_specified:
            cmdArgs.append('-g')
            cmdArgs.append(geneset_gtf)
          else:
             util.critical('Option "-cuff_gtf" should not be specified if "-g" option from Cufflinks has already been set in "-cuff_opt". Exiting...')
        cmdArgs.append(f)
        
        util.call(cmdArgs,stderr='cufflinks_stderr.log')
        rm_lines('cufflinks_stderr.log',util.LOG_FILE_PATH)
        
        # Rename output files 
        for i in range(4):
          ofc = out_folder + cuff_files[i]
          nn = header_cuff + cuff_files[i]
          os.rename(ofc, nn)
    
    
    fileObj_assemblies.close()
    
    # Run Cuffmerge
    
    cuff_head = samples_csv.split('/')[-1]
    ofc2 = out_folder + cuff_head + '_cuffmerge.gtf'
    
    if exists_skip(ofc2):
      util.call(['cuffmerge','--version'],stdout=util.LOG_FILE_OBJ)
      err = 0
      cmdArgs = ['cuffmerge', '-s',genome_fasta,
                 '-p',str(num_cpu),
                 '-o',out_folder]
      if is_gtf_specified:
        cmdArgs += cuff_gtf_file
        err = 1
      elif cuff_gtf is True:
        if err is 1:
          util.critical('Option "-cuff_gtf" should not be specified if "-g" option from Cufflinks has already been set in "-cuff_opt". Exiting...')
        cmdArgs.append('-g')
        cmdArgs.append(geneset_gtf)
      cmdArgs.append(assemblies)
      util.call(cmdArgs,stderr='cuffmerge_stderr.log')
      rm_lines('cuffmerge_stderr.log',util.LOG_FILE_PATH)
      os.rename(out_folder + 'merged.gtf', ofc2)
      
    # Run Cuffquant
    
    cxb_list=[]
    
    basic_options = ['-u',
                     '-b', genome_fasta,
                     '-p', str(num_cpu)]
    if library_type is not None:
      basic_options += library_type
    
    basic_options += ['-o', out_folder] # Output folder added to the end so to facilitate using this object in downstream code (cuffdiff and cuffnorm steps)
    
    for f in bam_files:
      f2 = f.split('/')[-1]
      ofc3 = out_folder + f2 + '_abundances.cxb'
      cxb_list.append(ofc3)
      
      if exists_skip(ofc3):
        report_cuff_version('cuffquant')
        cmdArgs = ['cuffquant'] + basic_options + [ofc2,f]
        util.call(cmdArgs,stderr='cuffquant_stderr.log')
        rm_lines('cuffquant_stderr.log',util.LOG_FILE_PATH)
        os.rename(out_folder + 'abundances.cxb', ofc3)
        
    
    # Set list of replicates and conditions for both Cuffnorm and Cuffdiff
    
    reps = [cxb_list[0]]
    reps_list = []
    conds = list(set(csv[:,3]))
    conds = ','.join(conds)
    
    
    for i in range(1,csv.shape[0]):
      if csv[i-1,3]==csv[i,3]:
        reps.append(cxb_list[i])
      else:
        reps_list.append(reps)
        reps = [cxb_list[i]]
    
    reps_list.append(reps)

    reps_list2 = []
    
    for reps in reps_list:
      reps = ','.join(reps)
      reps_list2.append(reps)
    
    # Run Cuffnorm
    
    if cuffnorm: 
      
      report_cuff_version('cuffnorm')
      
      dir_cnorm = out_folder + '/cuffnorm/'
      dir_cnorm = new_dir(dir_cnorm)

      cmdArgs = ['cuffnorm'] + basic_options[3:-1]
      cmdArgs.append(dir_cnorm)
      cmdArgs.append('-L')
      cmdArgs.append(conds)
      cmdArgs.append(ofc2)
      cmdArgs += reps_list2
      util.call(cmdArgs,stderr='cuffnorm_stderr.log')
      rm_lines('cuffnorm_stderr.log',util.LOG_FILE_PATH)

    
    # Run Cuffdiff
    
    report_cuff_version('cuffdiff')
    
    dir_cdiff = out_folder + '/cuffdiff/'
    dir_cdiff = new_dir(dir_cdiff)
      
    basic_options[4] = '1' # Cuffdiff should not be run in more than 1 thread to avoid crashing due to insufficient memory
    
    cmdArgs = ['cuffdiff'] + basic_options[:-1]
    cmdArgs.append(dir_cdiff)
    cmdArgs.append('-L')
    cmdArgs.append(conds)
    cmdArgs.append(ofc2)
    cmdArgs += reps_list2
    util.call(cmdArgs,stderr='cuffdiff_stderr.log')
    rm_lines('cuffdiff_stderr.log',util.LOG_FILE_PATH)

    # Run CummeRbund
    
    cmdArgs = ['Rscript','--vanilla', os.environ["cummeRbund"],dir_cdiff]
    util.call(cmdArgs)
    util.info('Plot saved in %s as exploratory_analysis_plots.pdf...' % dir_cdiff)
    
    
       
