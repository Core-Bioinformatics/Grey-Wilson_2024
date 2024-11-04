import scrublet as scr
import glob
import pandas as pd

filenamesList = glob.glob('Analysis/DoubletDetection/DoubletExpMat/*.csv')
for myfile in filenamesList:
  print(myfile)
  current_sample = myfile.split(sep='/')[6].split(sep='.')[0]
  counts_matrix = pd.read_csv(myfile)
  counts_matrix = counts_matrix.set_index('Unnamed: 0')
  print('Read')
  scrub = scr.Scrublet(counts_matrix)
  doublet_scores, predicted_doublets = scrub.scrub_doublets()
  print(predicted_doublets[:5])
  print(current_sample+counts_matrix.index[:5])
  mydf = pd.DataFrame({'CellID':current_sample+"_"+counts_matrix.index,'Doublet':predicted_doublets})
  mydf.to_csv('Analysis/DoubletDetection/doublets_scrublet.csv',mode='a',index=False, header=False)
  
import doubletdetection
filenamesList = glob.glob('Analysis/DoubletDetection/DoubletExpMat/*.csv')
for myfile in filenamesList:
  print(myfile)
  current_sample = myfile.split(sep='/')[6].split(sep='.')[0]
  counts_matrix = pd.read_csv(myfile)
  counts_matrix = counts_matrix.set_index('Unnamed: 0')
  counts_matrix = counts_matrix.loc[:, (counts_matrix != 0).any(axis=0)]
  print('Read')
  clf = doubletdetection.BoostClassifier()
  labels = clf.fit(counts_matrix).predict()
  scores = clf.doublet_score()
  mydf = pd.DataFrame({'CellID':current_sample+"_"+counts_matrix.index,'Doublet':labels,'Score':scores})
  mydf.to_csv('Analysis/DoubletDetection//doublets_doubletdetection.csv',mode='a',index=False, header=False)
