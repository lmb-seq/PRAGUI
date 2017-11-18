#!/home/paulafp/applications/anaconda3/bin/python3

import csv
import os
import shutil
import uuid

import cross_fil_util as util


PROG_NAME = 'CAT_FASTQ'
DESCRIPTION = 'Function to concatenate fastq files from different lanes and flowcells.'

def cat_fastq(barcode_csv, fastq_paths_r1,
                  fastq_paths_r2=None, out_top_dir=None, 
                  sub_dir_name=None, file_ext=None):
  
  if not sub_dir_name:
    sub_dir_name = 'strain'
  
  if fastq_paths_r2:
    fastq_paths_r2_list = fastq_paths_r2
  else:
    fastq_paths_r2_list = []
  
  for file_path in [barcode_csv] + fastq_paths_r1 + fastq_paths_r2_list or []:
    is_ok, msg = util.check_regular_file(file_path)
  
    if not is_ok:
      util.critical(msg)
      
  if not file_ext:
    file_ext = util.get_file_ext(fastq_paths_r1[0])
    
  if not out_top_dir:
    file_path = fastq_paths_r1[0]
    out_top_dir = os.path.dirname(file_path)
    
  # # Concatenate FASTQ files 
  
  sample_barcodes = {}
  barcode_samples = {}
  
  # Read CVS
  with open(barcode_csv, 'rU') as file_obj:
    csv_data = csv.reader(file_obj)
    header = next(csv_data)
    
    for seq_run_id, barcode_name, barcode_seq, sample_name in csv_data:
      barcode_name = barcode_name.replace('-','_')
      sample_name = sample_name.replace(' ','_')
      
      if sample_name in barcode_samples:
        util.critical('Multiple samples/strains with name "%s" present in CSV file %s' % (sample_name, barcode_csv))
      
      barcode_samples[sample_name] = barcode_name
      sample_barcodes[barcode_name] = (seq_run_id, sample_name)
  
  # Make subdirs
  for sample_name in barcode_samples:
    dir_name = os.path.join(out_top_dir, sub_dir_name, sample_name)
    util.makedirs(dir_name, exist_ok=True)
  
  # Concatenate FASTQ files with the same barcode for each read 
  # and save results in corresponding strain folder
  # - now makes a symbolic link of only one file for a strain and not gzipped
  
  strain_fastq_paths = {}
  
  for barcode_name in sample_barcodes:
    seq_run_id, sample_name = sample_barcodes[barcode_name]
    file_pattern = '%s*%s*' % (seq_run_id, barcode_name)
    
    if fastq_paths_r2:
      out_file_name_1 = '%s_r_1%s' % (sample_name, file_ext)
      out_file_name_2 = '%s_r_2%s' % (sample_name, file_ext)
      in_fastq_paths_1 = util.match_files(fastq_paths_r1, file_pattern) # Read pairs files already separated
      in_fastq_paths_2 = util.match_files(fastq_paths_r2, file_pattern)
      
      io_paths =  [(in_fastq_paths_1, out_file_name_1),
                   (in_fastq_paths_2, out_file_name_2)] 
  
    else:
      out_file_name  = '%s%s' % (sample_name, file_ext)
      in_fastq_paths = util.match_files(fastq_paths_r1, file_pattern)
      
      io_paths =  [(in_fastq_paths, out_file_name)]
    
    fastq_paths = []
       
    for in_fastq_paths, out_file_name in io_paths:
      if not in_fastq_paths:
        util.critical('No FASTQ read files found for run %s barcode %s' % (seq_run_id, barcode_name))
      
      out_fastq_path = os.path.join(out_top_dir, sub_dir_name, sample_name, out_file_name)

      if os.path.exists(out_fastq_path):
        util.warn('FASTQ file %s already exists and won\'t be overwritten...' % out_fastq_path)
      
      else:
        # Concatenate or sym link
        
        if len(in_fastq_paths) == 1 and not in_fastq_paths[0].endswith('.gz'):
          util.info('Sym linking %s reads to %s' % (barcode_name, out_fastq_path))
          os.symlink(in_fastq_paths[0], out_fastq_path)
        
        else:
          with open(out_fastq_path, 'wb') as out_file_obj:
            util.info('Concatenating %s reads to %s' % (barcode_name, out_fastq_path))
            
            for fastq_path in in_fastq_paths:
              shutil.copyfileobj(util.open_file(fastq_path, 'rb'), out_file_obj) # Accepts GZIP input
 
      fastq_paths.append(out_fastq_path)
    
    strain_fastq_paths[sample_name] = fastq_paths
    
  return strain_fastq_paths
  



if __name__ == '__main__':

  from argparse import ArgumentParser
  
  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             prefix_chars='-', add_help=True)
  
  arg_parse.add_argument('barcode_csv', metavar='BARCODE_CSV_FILE',
                         help='CSV format file containing barcode strain/sample names') 

  arg_parse.add_argument('fastq_paths', nargs='+', metavar='FASTQ_FILES',
                         help='File paths of FASTQ sequence read files (may contain wildcards).')
  
  arg_parse.add_argument('-outdir', metavar='DIR_NAME', default=None,
                         help='Name of directory for output files and sub-directories (need not exist). Defaults to where first FASTQ file is located') 

  arg_parse.add_argument('-pe', nargs=2, metavar='PAIRED_READ_TAGS', default=['r_1','r_2'],
                         help='The subtrings/tags which are the only differences between paired FASTQ file paths. Default: r_1 r_2') 
                         
  arg_parse.add_argument('-se', default=False, action='store_true',
                         help='Input reads are single-end data, otherwise defaults to paired-end.')

  arg_parse.add_argument('-sub_dir_name', default=None, 
                         help='Name of subdirectory created in DIR_NAME to store output files (need not exist). Defaults to "strain".')

  args = vars(arg_parse.parse_args())

  barcode_csv   = args['barcode_csv']
  fastq_paths   = args['fastq_paths']
  pair_tags     = args['pe']
  is_single_end = args['se']
  out_top_dir   = args['outdir']
  sub_dir_name  = args['sub_dir_name']
  
  if len(pair_tags) != 2:
    util.critical('When specified, exactly two paired-end filename tags must be given.')
  
  if is_single_end:
    fastq_paths_r1 = fastq_paths
    fastq_paths_r2 = None
  
  else:
    if not fastq_paths:
      util.critical('No FASTQ files found')

    elif len(fastq_paths) % 2 != 0:
      util.critical('When using paired-end data an even number of FASTQ files must be specified.')
  
    fastq_paths_r1, fastq_paths_r2 = util.pair_fastq_files(fastq_paths, pair_tags)
  
  
  cat_fastq(barcode_csv, fastq_paths_r1,
                  fastq_paths_r2=fastq_paths_r2, out_top_dir=out_top_dir, 
                  sub_dir_name=None)
