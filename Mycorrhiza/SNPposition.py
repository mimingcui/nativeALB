from mycorrhiza.dataset import Myco, Structure
from sklearn.model_selection import train_test_split
from mycorrhiza.analysis import CrossValidate
from mycorrhiza.plotting import mixture_plot
from mycorrhiza.settings import const
from sklearn.metrics import mutual_info_score
from sklearn import preprocessing
from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels
from sklearn.model_selection import StratifiedKFold
import seaborn as sn
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from matplotlib import figure
import pickle
from tqdm import tqdm
import random
import pandas as pd
import sys

#================================================================================================================ SNP Selection
if not os.path.exists('results'):
    os.makedirs('results')

vcf = sys.argv[1]
name = sys.argv[2]
strFile = sys.argv[3]
mergeFile = sys.argv[4]
num_selected_snps = int(sys.argv[5])
threshold = 0.6

num_plot_node = 21

mergeDict = {}
if mergeFile != 'None':
	df = pd.read_csv(mergeFile, header = None)
	for pop, mergedPop in zip(df.iloc[:,0], df.iloc[:,1]):
	    mergeDict[pop.strip()] = mergedPop.strip()

rawdata = pd.read_csv(strFile, sep = ' ',  header = None)
if mergeFile != 'None':
    rawdata.iloc[:,1] = rawdata.iloc[:,1].map(mergeDict)
rawdata.to_csv('./{}_collapsed.str'.format(name), sep = ' ', index = None, header = None)
print(rawdata.iloc[:,1])
print(set(rawdata.iloc[:,1]))
popList = list(set(rawdata.iloc[:,1]))
print(popList)
popList.sort()
# print(rawdata.iloc[:,1])
halfdata = rawdata.iloc[::2]
rawdata.iloc[:,1].to_csv('./{}.pop'.format(name), sep = '\t', index = None,  header = None)
rawdata.iloc[:,1].to_csv('./{}.csv'.format(name),  index = None,  header = None)
# print(set(popList))
# print(rawdata)

label = rawdata.iloc[:,:2]
data = rawdata.iloc[:,2:-1]

def spliceData(start, end, data):
    reducedX = data.iloc[:,start:end]
    splicedX = reducedX
    reducedX = pd.concat([data.iloc[:,:2], reducedX], axis=1)
    # print(reducedX)

    reducedX.to_csv('./test.vcf.str', sep=' ', index=False, header=False)

    finalData = Structure(file_path='./test.vcf.str')
    finalData.load()

    cv = CrossValidate(dataset=finalData, out_path='./results')
    result = cv.run(n_partitions=1, n_loci=0, n_splits=5, n_estimators=60, n_cores=1)

    accuracy = result.output_accuracy()
    mixture_plot(cv)
    return splicedX

def modifiedFS(data, label, X_train, y_train, num_features, threshold):
    sortedMIList = []
    
    labelDict = {} # Stores the indices of high MI features according to labels
    sample_index = [i for i in range(data.shape[0])]
    for pop in popList:
        mutualInfos = []
        pop_label = y_train[1].copy()
        pop_label.loc[pop_label.iloc[:] != pop] = 0
        pop_label.loc[pop_label.iloc[:] == pop] = 1
        print(pop_label.sum())
        for i in range(X_train.shape[1]):
            mutualInfos.append(mutual_info_score(pop_label, X_train.iloc[:,i]))
        sortedMI = np.array(mutualInfos).argsort()
        sortedMI = sortedMI[::-1]
        print('This is the sorted MI index for {} :'.format(pop) , sortedMI)
        sortedMIList.append(sortedMI)
        labelDict[pop] = sortedMI
    print(labelDict)
    sortedMI = []
    features = pd.DataFrame()
    MI_indices_pop = pd.DataFrame.from_dict(labelDict)
    [row, column] = MI_indices_pop.shape
    for i in range(row) :
        for j in range(column):
            sortedMI.append(MI_indices_pop.iloc[i, j])
    print(len(sortedMI))
    features = pd.concat([data.iloc[:,sortedMI[0]], data.iloc[:,sortedMI[1]]], axis=1)
    featureIndexes = [sortedMI[0], sortedMI[1]]
    accuracyList = []
    idx = 2
    while len(featureIndexes) < num_features:
        print(len(sortedMI))
        while sortedMI[idx] in featureIndexes:
            idx = idx+1
        print('Num iteration : ', idx)
        cor = False
        for i in range(features.shape[1]):
            if abs(features.iloc[:,i].corr(data.iloc[:,sortedMI[idx]])) > threshold:
                cor = True
                print(cor)
                break
        if cor:
            idx = idx +1
            continue
            
            
        else :
            features = pd.concat([features, data.iloc[:,sortedMI[idx]]], axis =1)
            featureIndexes.append(sortedMI[idx])
            idx = idx +1
        print('Num selected SNPs :',len(featureIndexes))
        
               
    return pd.concat([label, features],axis=1), accuracyList, featureIndexes

skf = StratifiedKFold(n_splits=2)
skf.get_n_splits(data, label[1])
print(data)
print(label[1])  

selected_features, accuracies , featureIndexes = modifiedFS(data, label, data,  label,  num_selected_snps, threshold)
selected_features.to_csv('selected_full.csv', sep=',', index=False, header=False)
selected_features.to_csv('selected_full.str', sep=' ', index=False, header=False)
print(featureIndexes)
with open('fullFeatureIndexes_{}.pkl'.format(name), 'wb') as f:
    pickle.dump(featureIndexes, f)
try:
    os.remove("./*.nex")
    print("network files removed")
except OSError:
    pass

os.system("cut -f 1,2 {} > SNPID.txt".format(vcf))

SNPID = []
with open('SNPID.txt', 'r') as input:
   for line in input:
       if '#' not in line:
            SNPID.append(line.strip())

# with open('SNPID.txt', 'r') as f:
#     SNPID = [line.strip() for line in f]
SNPID = [0]+SNPID
print(SNPID)
print(len(SNPID))
print(data.shape)
# ALB = [0]+ALB
# AGM = [0]+AGM
# print(len(AGM))
# print(len(ALB))
# print(AGM)
# print(ALB)
             
selected_ID = [ SNPID[i] for i in featureIndexes]
   
selected_ID= pd.DataFrame(selected_ID)
             
selected_ID.to_csv('results/{}_selected_ID.csv'.format(name), header = False, index = False)


print(selected_ID.shape)
