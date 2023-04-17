def detectExistDirectory(tempDir):
	if os.path.exists(tempDir):
		shutil.rmtree(tempDir)

	os.makedirs(tempDir)
	return None

def generateIterativeFeatureFile(numbers, inputFilename):
	#g = pd.read_csv("out_%d.svm"%(numbers))
	g = "out_%d.svm"%(numbers)
	f = pd.read_csv(inputFilename)
	f.iloc[:,:numbers].to_csv(g,index=0)
	return None

def removeFile(startFilename):
	fileList = os.listdir(os.curdir)
	for eachFile in fileList:
		if re.search(r'^%s'%startFilename, eachFile) == None:
			pass
		else:
			os.remove(eachFile)

def childProcessing(numbers, inputFilename):
	generateIterativeFeatureFile(numbers, inputFilename)
	startFilename = "{0}/out_{1}.svm".format(tempDir, numbers)
	commands = []
	commands.append("python {0}/train_all.py {1}".format(pathPrefix, startFilename)) # 5 cross vaild

	for eachCmd in commands:
		os.system(eachCmd)

	outResult = os.getcwd() + os.path.sep + "outResult.txt"
	df = pd.read_csv('{0}.out'.format(startFilename), sep='\t')
	data = [startFilename]
	acc = df.acc.to_list()
	data.extend(acc)
	ndf = pd.DataFrame(columns=['file','gbn','knn','svm', 'RF'],data=[data])
	#tempLine = open(r'%s.out'%(startFilename)).readlines()[0]
	ndf.to_csv(outResult, header=False, index=False, mode='a', sep='\t')
	print("	***%s - Finished!***"%startFilename)

	removeFile(startFilename)

def mainProcessing(inputFilename, tempDirect, cpuNum):
	os.chdir(tempDirect)
	maxFeatureNum = 3040#obtainFileMaxFeatureNumber(inputFilename)
	featureNumList = list(range(1,maxFeatureNum+1))
	
	num = 0
	threads = []
	outResult = os.getcwd() + os.path.sep + "outResult.txt"
	with open(outResult, 'w') as f:
		f.write('file\tgnb\tknn\tsvm\n')
	for fInde in range(1, 3040):
		featureNum = featureNumList[fInde*1]
		num += 1

		t = threading.Thread(target=childProcessing, args=(featureNum, inputFilename))
		threads.append(t)
		t.start()

		if (num == cpuNum) or (featureNum == maxFeatureNum):
			#for t in threads:
			#	t.start()
			for t in threads:
				t.join()

			num = 0
			threads = []
		else:
			pass
	os.chdir(os.pardir)


import os
import re
import sys
import shutil
import threading
import pandas as pd

pathPrefix = os.getcwd() + os.path.sep
paraInput = sys.argv[1]
inSvmFeatureFile = pathPrefix + paraInput
tempDir = pathPrefix + "temp"
cpuCount = 1

if __name__ == '__main__':
	import sys
	each = sys.argv[1]
	detectExistDirectory(tempDir)
	mainProcessing(inSvmFeatureFile, tempDir, cpuCount)

"""
Useful:
	python3 IFS.py featureSetSortedFile.svm
"""
