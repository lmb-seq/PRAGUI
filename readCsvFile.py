#!/lmb/home/paulafp/applications/anaconda3/bin/python3

# Define function to read csv files
# (Function taken from T. J. Stevens and W. Boucher's python book)

import csv
import sys
import numpy

def readCsvFile(filename, converters=None, separator=',', newline='',header=True):

  dataList = []

  if sys.version_info.major >2:
    fileObj = open(filename, 'r', newline=newline)
  else:
    fileObj = open(filename,'rb')

  reader = csv.reader(fileObj, delimiter=separator)

  for n, row in enumerate(reader):
    if header:
      if n > 0: # n= 0 is the header, which we ignore
        for index, datum in enumerate(row):
          if converters:
            print(converters)
            convertFunc = converters[index] # The converters argument should be a list
                                            # of functions (one for each column)
                                            # and is meant to transform the data into the 
                                            # the correct format (e.g. int, float, etc.).
                                            # In case one column should not be converted,
                                            # just use None for that column.
            if convertFunc:
              row[index] = convertFunc(datum)
        
        dataList.append(row)
    else: 
      for index, datum in enumerate(row):
        if converters:
          print(converters)
          convertFunc = converters[index]
          if convertFunc:
            row[index] = convertFunc(datum)
      
      dataList.append(row)
  
  
  fileObj.close()

  return numpy.array(dataList)
