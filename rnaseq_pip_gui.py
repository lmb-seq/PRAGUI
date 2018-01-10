import uuid
import os
import sys
from appJar import gui
import rnaseq_pip_util as rnapip

current_path = os.path.realpath(__file__)
current_path = os.path.dirname(current_path) + '/cell_bio_util'

sys.path.append(current_path)
import cell_bio_util as util



def test(samples_csv, genome_fasta, genome_gtf, geneset_gtf=None, analysis_type=['DESeq','Cufflinks'][0], trim_galore=None, 
                       skipfastqc=False, fastqc_args=None, aligner=rnapip.DEFAULT_ALIGNER, is_single_end=False, pair_tags=['r_1','r_2'],
                       star_index=None,star_args=None,num_cpu=8,mapq=20,stranded=False,contrast='condition',levels=None,
                       cuff_opt=None, cuff_gtf=False,cuffnorm=False, python_command=None,q=False,log=False):
  print(locals())


app=gui()


def replace_key(dic,old_key, new_key):
  if old_key in dic:
    dic[new_key]=dic[old_key]
    del(dic[old_key])
    return(dic)

def hide(win):
  # app.hideSubWindow('csv')
  app.destroySubWindow('csv')

def next(button,win):
  if button == 'Cancel':
    app.stop()
  else:
    app.addButtons(win,launch)

def add_pe_tags(btn):
   if btn == 'paired-end':
     app.addLabelEntry('tags')
     
def submit(btn):
  if btn == 'cancel':
    app.stop()
  else:
    args = app.getAllEntries()
    args_copy = args.copy()    
    if app.getOptionBox('Samples File') == 'create':
      if 'samples_csv' not in args.keys():
        app.addLabelNumericEntry('samples_csv',1,1)
        app.setEntryDefault('samples_csv')
        app.setLabelTooltip('samples_csv','Number of samples to be analysed.')
        app.addButton('csv',launch,1,2)
      else:
        if not csv_created:
          app.errorBox('Missing Input','Please provide sample information by pressing the csv button.')
          return()
        else:
          args['samples_csv'] = csv_file
    if app.getOptionBox('Samples File') == 'upload' and 'samples_csv' not in args.keys():
      app.addFileEntry('samples_csv',1,1)
      app.infoBox('Info','Please add path to csv file.')
      return()
    if app.getOptionBox('Library type') == 'paired-end' and 'pair_tags' not in args.keys():
      app.addLabelEntry('pair_tags', 3, 1)
      app.setLabelTooltip('pair_tags','Substrings/tags which are the only differences between paired FASTQ file paths. e.g.: r_1 r_2.')  
    else:
      if 'num_cpu' in args:
        args['num_cpu'] = int(args['num_cpu'])
      for key, val in args_copy.items():
        if val in ['',0]:
          del(args[key])    
      # args['samples_csv'] = csv_file
      analysis_type = app.getOptionBox('analysis_type')
      args['analysis_type'] = analysis_type
      library_type = app.getOptionBox('Library type')
      if library_type == 'single-end':
        args['is_single_end'] = True
      checkBox = app.getAllCheckBoxes()
      if 'qsub' in checkBox:
        del(checkBox['qsub'])
        args.update(checkBox)
        for old_key, new_key in [['num_cpu','cpu'],['pair_tags','pe'],['is_single_end','se']]:
          replace_key(args,old_key,new_key)
        command = 'python3 /lmb/home/paulafp/applications/RNAseq_pipeline//rnaseq_pip_util.py '
        for key, item in args.items():
          k_i = '-%s %s' % (key,item)
          command = command + k_i
          qsubArgs = ['echo',command,'|','qsub','-cwd','-pe','smp','4','-j','y']
          util.call(qsubArgs)
          app.infoBox('Info','Running Pipeline via qsub on the LMB cluster.')
          return()
      else:
        args.update(checkBox)
        rnapip.rnaseq_diff_caller(**args)
      #test(**args)


def launch(win):
  samples = int(app.getEntry('samples_csv'))
  app.startSubWindow(win)
  app.addLabel('samples','samples',0,0)
  app.addLabel('read 1','read 1',0,1)
  app.addLabel('read 2','read 2',0,2)
  app.addLabel('condition','condition',0,3)
  global samplesList
  global read1List
  global read2List
  global conditionList
  samplesList = []
  read1List   = []
  read2List   = []
  conditionList   = []
  for i in range(samples):
    row       = app.getRow()
    sample    = 'sample'    + str(i)
    read1     = 'read1_'    + str(i)
    read2     = 'read2_'    + str(i)
    condition = 'condition' + str(i)
    app.addEntry(sample,row,0)
    app.addFileEntry(read1,row,1)
    app.addFileEntry(read2,row,2)
    app.addEntry(condition,row,3)
    samplesList.append(sample)
    read1List.append(read1)
    read2List.append(read2)
    conditionList.append(condition)
  app.addButton('Create CSV',create_csv)
  app.addButton('CLOSE',hide,row+1,3)
  app.stopSubWindow()
  app.showSubWindow(win)

def create_csv(btn):
  global csv_file
  csv_file = '%s/samples_%s.csv' % (os.getcwd(),str(uuid.uuid4())[:8])
  fileObj = open(csv_file,'w')
  header = 'samples\tread1\tread2\tcondition\n'
  fileObj.write(header)
  for i in range(len(samplesList)):
    s_info = [
    app.getEntry(samplesList[i])   ,
    app.getEntry(read1List[i])     ,
    app.getEntry(read2List[i])     ,
    app.getEntry(conditionList[i]) ]
    line = '\t'.join(s_info) + '\n'
    fileObj.write(line)
  fileObj.close()
  global csv_created
  csv_created = True
  for entry in samplesList + read1List + read2List + conditionList:
    app.clearEntry(entry,callFunction=False)


csv_created = False

app.addLabel('title', 'RNAseq Pipeline')
app.setLabelBg('title','lightblue')
# app.addLabelNumericEntry('samples_csv')
# app.setEntryDefault('samples_csv')
# app.setLabelTooltip('samples_csv','Number of samples to analysed.')
app.addLabelOptionBox('Samples File',['create','upload'])
app.addLabelOptionBox('analysis_type',['DESeq','Cufflinks'])
app.addLabelOptionBox('Library type',['paired-end','single-end'])
app.addLabel('genome_fasta','Genome fasta file')
app.addFileEntry('genome_fasta')
app.addLabel('genome_gtf','Genome gtf file')
app.addFileEntry('genome_gtf')
app.addLabel('geneset_gtf','Geneset file')
app.addFileEntry('geneset_gtf')
app.addLabelEntry('trim_galore')
app.setEntryDefault('trim_galore', 'None')
app.addLabelEntry('fastqc_args')
app.setEntryDefault('fastqc_args', 'None')
app.addLabelEntry('star_args')
app.setEntryDefault('star_args', 'None')
app.addLabelNumericEntry('mapq')
app.setEntryDefault('mapq', 20)
app.addLabelNumericEntry('num_cpu')
app.setEntryDefault('num_cpu', 0)
app.addLabelEntry('contrast')
app.setEntryDefault('contrast', 'condition')
app.addLabelNumericEntry('contrast_levels')
app.setEntryDefault('contrast_levels', 0)
app.addLabelEntry('cuff_opt')
app.setEntryDefault('cuff_opt', 'None')
app.addCheckBox('stranded')
app.addCheckBox('skipfastqc')
app.addCheckBox('cuff_gtf')
app.addCheckBox('cuffnorm')
app.addCheckBox('q')
app.addCheckBox('log')
app.addCheckBox('qsub')

app.addButtons(['submit','cancel'],submit)

# app.addButton('csv',launch)


app.go()
