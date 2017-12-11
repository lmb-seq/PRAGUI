import uuid
import os
from appJar import gui

app=gui()


def hide(win):
  app.hideSubWindow('csv')

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
    if not csv_created:
      app.errorBox('Missing Input','Please provide sample information by pressing the csv button.')
      app.stop()   
    elif app.getOptionBox('Library type') == 'paired-end' and 'pe' not in args.keys():
      app.addLabelEntry('pe', 3, 1)
      app.infoBox('Paired-end Reads','Please provide substrings/tags which are the only differences between paired FASTQ file paths. e.g.: r_1 r_2. And press submit.')  
    else:
      fileObj = open('/home/paulafp/Desktop/test.txt','w')
      global commandArgs
      commandArgs = []
      for key, val in args.items():
        if val not in ['',0] and key not in ['samples_csv','genome_fasta']:
          flag = '-%s' % key
          commandArgs += [flag,str(val)]
      analysis_type = app.getOptionBox('analysis_type')
      commandArgs += ['-analysis_type',analysis_type]
      checkBox = app.getAllCheckBoxes()
      if app.getOptionBox('Library type') == 'single-end':
        commandArgs.append('-se')
      for key, val in checkBox.items():
        if checkBox[key]:
          flag = '-%s' % key
          commandArgs.append(flag)
      #samples_csv   = app.getEntry('samples_csv')
      #fileObj.write('\n' + analysis_type + '\t' + str(samples_csv))
      #stranded = app.getCheckBox('stranded')
      #if stranded:
      #  fileObj.write('\nstranded=True')
      fileObj.write(' '.join(commandArgs))
      fileObj.close()

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
app.addLabelNumericEntry('samples_csv')
app.setEntryDefault('samples_csv')
app.setLabelTooltip('samples_csv','Number of samples to analysed.')
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
app.addLabelNumericEntry('cpu')
app.setEntryDefault('cpu', 0)
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

app.addButtons(['submit','cancel'],submit)

app.addButton('csv',launch)


app.go()
