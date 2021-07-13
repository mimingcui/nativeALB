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

name = sys.argv[1]
strFile = sys.argv[2]
mergeFile = sys.argv[3]
num_selected_snps = int(sys.argv[4])
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

fullData = Structure(file_path='./{}_collapsed.str'.format(name))
fullData.load()
cv = CrossValidate(dataset=fullData, out_path='./results')
resultFull = cv.run(n_partitions=1, n_loci=0, n_splits=2, n_estimators=60, n_cores=1)
accuracyFull = resultFull.output_accuracy()

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

def modifiedFS(data, label, X_train, X_test, y_train, y_test, train_index, test_index, num_features, threshold):
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
        print(len(featureIndexes))
        
               
    return pd.concat([label, features],axis=1), accuracyList, train_index, test_index, featureIndexes

skf = StratifiedKFold(n_splits=2)
skf.get_n_splits(data, label[1])
print(data)
print(label[1])  
i = 0
for train_index, test_index in skf.split(data, label[1]):
    print("TRAIN:", train_index, "TEST:", test_index)
    X_train, X_test = data.iloc[train_index,:], data.iloc[test_index, :]
    y_train, y_test = label.iloc[train_index, :], label.iloc[test_index, :]
    selected_features, accuracies , train_index, test_index, featureIndexes = modifiedFS(data, label, X_train, X_test, y_train, y_test, train_index, test_index, num_selected_snps, threshold)
    selected_features.to_csv('selected_mixed{}.csv'.format(i), sep=',', index=False, header=False)
    selected_features.to_csv('selected_mixed{}.str'.format(i), sep=' ', index=False, header=False)
    print(featureIndexes)
    with open('featureIndexes_{}_{}.pkl'.format(name, i), 'wb') as f:
        pickle.dump(featureIndexes, f)
    with open('trainIndexes_{}_{}.pkl'.format(name, i), 'wb') as f:
        pickle.dump(train_index, f)
    with open('testIndexes_{}_{}.pkl'.format(name, i), 'wb') as f:
        pickle.dump(test_index, f)
    i = i+1
    try:
        os.remove("./*.nex")
        print("network files removed")
    except OSError:
        pass

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
try:
    os.remove("./*.nex")
    print("network files removed")
except OSError:
    pass

def readIndices(file, fold):
    ranks = []
    for i in range(fold):
        with open(file.format(name, i), 'rb') as f:
            selec = pickle.load(f)
        ranks.append(selec)
    return ranks


train_indices = readIndices('trainIndexes_{}_{}.pkl', 2)
test_indices = readIndices('testIndexes_{}_{}.pkl', 2)
train_index1 = train_indices[0]
train_index2 = train_indices[1]
test_index1 = test_indices[0]
test_index2 = test_indices[1]
def acc_curve(file, test_index, index):
    selected_features = pd.read_csv(file, sep=',', header=None)
    selected_features2 = selected_features.iloc[test_index, :]
    selected_features2.iloc[:,:3+index].to_csv('temp.str', sep=' ', index=False, header=False)
    finalData = Structure(file_path='temp.str')
    finalData.load()
    cv = CrossValidate(dataset=finalData, out_path='./results')
    result = cv.run(n_partitions=1, n_loci=0, n_splits=2, n_estimators=60, n_cores=1)
    accuracy = result.output_accuracy()
    y_pred = result.pred_populations
    y_true = result.real_populations
    label = result.q_populations
    
    mat=confusion_matrix(y_true, y_pred, labels=label)
    rowSum = np.sum(mat, axis=1)
    
    accuList = []
    for i in range (len(rowSum)):
        accuList.append(float(mat[i][i])/rowSum[i])
    return accuList, accuracy
        

def plot_heat(file, index, test_index):
    selected_features = pd.read_csv(file, sep=',', header=None)
    selected_features.iloc[:,:3+index].to_csv('temp.str', sep=' ', index=False, header=False)
    finalData = Structure(file_path='temp.str')
    finalData.load()
    cv = CrossValidate(dataset=finalData, out_path='./results')
    result = cv.run(n_partitions=1, n_loci=0, n_splits=2, n_estimators=60, n_cores=1)
    accuracy = result.output_accuracy()
    y_pred = result.pred_populations
    y_true = result.real_populations
    label = result.q_populations
    print(label)
    mat=confusion_matrix(y_true, y_pred, labels=label)
    df_cm = pd.DataFrame(mat, index = [i for i in label], columns = [i for i in label])
    #print(df_cm)
    plt.figure(figsize = (10,7))

    sn.heatmap(df_cm, annot=True)
    # remove these later
    b, t = plt.ylim() # discover the values for bottom and top
    b += 0.5 # Add 0.5 to the bottom
    t -= 0.5 # Subtract 0.5 from the top
    plt.ylim(b, t) # update the ylim(bottom, top) values
    # plt.show() # ta-da!
    # plt.show(block = False);
    return mat

def get_accuracies(data, test_index,  index):
    realAccuList = []
    realAccuList2 = []
    wholeAccuList = []
    wholeAccuList2 = []
    for i in range (1, index):
        realAccu2 , wholeAccu2 =acc_curve(data, test_index, int(num_selected_snps/(index-1))*i)
        realAccuList2.append(realAccu2)
        wholeAccuList2.append(wholeAccu2)
        print(wholeAccuList2)
    return realAccuList2, wholeAccuList2


realAccuList, wholeAccuList = get_accuracies('selected_mixed0.csv', test_index1, num_plot_node)
realAccuList2, wholeAccuList2 =get_accuracies('selected_mixed1.csv', test_index2, num_plot_node)

def average_eachFeat(realAccuList, realAccuList2):
    result = []
    for a, b in zip(realAccuList, realAccuList2) :
        result.append([sum(x)/2 for x in zip(a,b)])
    return result
    
realAccuSum = average_eachFeat(realAccuList, realAccuList2)
wholeAccuSum = [sum(x)/2 for x in zip(wholeAccuList,wholeAccuList2)]

plt.figure()
plt.plot([i*int(num_selected_snps/(num_plot_node-1)) for i in range(1,num_plot_node)], wholeAccuSum, 'cyan', label='train')
plt.xlabel('Number of SNPs selected')
plt.ylabel('Accuracy')
plt.hlines(accuracyFull, 0, num_selected_snps, colors='k', linestyles='solid', label='Accuracy with all SNPs')
plt.title('{} Prediction Accuracy with different number of SNPs selected'.format(name))
plt.savefig('results/Accuracy when different Number of SNPs are selected')
plt.close()

def plot_eachFeat(realAccuList):
    x = [i for i in range(num_plot_node+1)]
    y = list(map(list, zip(*realAccuList)))
    plt.figure()
    for i in range(len(y)):
        timefilteredForce = plt.plot([i*int(num_selected_snps/(num_plot_node-1)) for i in range(1,num_plot_node)], y[i], label = popList[i])
        timefilteredForce = plt.xlabel('Number of Features')
        timefilteredForce = plt.ylabel('Accuracy')
        plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.)
    plt.savefig('results/each feature accuracy')
    plt.close()

plot_eachFeat(realAccuSum)

def prob_distribution(file, test_index, index, output):
    selected_features = pd.read_csv(file, sep=',', header=None)
    selected_features2 = selected_features.iloc[test_index, :]
    selected_features2.iloc[:,:3+index].to_csv('temp.str', sep=' ', index=False, header=False)
    finalData = Structure(file_path='temp.str')
    finalData.load()
    cv = CrossValidate(dataset=finalData, out_path='./results')
    result = cv.run(n_partitions=1, n_loci=0, n_splits=2, n_estimators=60, n_cores=1)
    q_matrix, q_populations, real_populations, identifiers = mixture_plot(cv)
    q_matrix = pd.DataFrame(data = q_matrix, index =real_populations , columns=q_populations)
    q_matrix.to_csv(output, sep=',', index=set(real_populations), header=q_populations)
    accuracy = result.output_accuracy()
    y_pred = result.pred_populations
    y_true = result.real_populations
    label = result.q_populations
    mat=confusion_matrix(y_true, y_pred, labels=label)
    rowSum = np.sum(mat, axis=1)
    
    accuList = []
    for i in range (len(rowSum)):
        accuList.append(float(mat[i][i])/rowSum[i])
    return accuList, accuracy

prob_distribution('selected_mixed0.csv', test_index1, num_selected_snps, 'results/test1.csv')
prob_distribution('selected_mixed1.csv', test_index2, num_selected_snps, 'results/test2.csv')

print(realAccuSum[10])
y_pos = np.arange(len(popList))
performance = realAccuSum[-1]

plt.figure(figsize=(15,3))
plt.bar(y_pos, performance, align='center', alpha=0.8)
plt.xticks(y_pos, popList)
plt.ylabel('Accuracy')
plt.title('Accuracy of each class with {} SNPs'.format(num_selected_snps))
plt.savefig('results/Accuracy of each population')
plt.close()

matrix0 = plot_heat('selected_mixed0.csv', num_selected_snps, test_index1)
matrix1 = plot_heat('selected_mixed1.csv', num_selected_snps, test_index2)

result = np.add(matrix0, matrix1)
# print repr(result)
df_cm = pd.DataFrame(result, index = [i for i in popList], columns = [i for i in popList])
plt.figure(figsize = (10,7))
sn.heatmap(df_cm, annot=True)
b, t = plt.ylim() # discover the values for bottom and top
b += 0.5 # Add 0.5 to the bottom
t -= 0.5 # Subtract 0.5 from the top
plt.ylim(b, t) # update the ylim(bottom, top) values
plt.savefig('results/Confusion Matrix')
plt.close()

fullData = Structure(file_path='./{}_collapsed.str'.format(name))
fullData.load()
cv = CrossValidate(dataset=fullData, out_path='./results')
resultFull = cv.run(n_partitions=1, n_loci=0, n_splits=2, n_estimators=60, n_cores=1)
mixture_plot(cv)
accuracyFull = resultFull.output_accuracy()
