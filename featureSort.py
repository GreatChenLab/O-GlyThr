import sklearn
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
import os
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix
import pickle

def gnb(x,y,x_test,y_test):

	from sklearn.naive_bayes import GaussianNB

	gnb = GaussianNB().fit(x,y)
	gnb_prob = gnb.predict(x_test)

	auroc = roc_auc_score(y_test, gnb.predict_proba(x_test)[:,1])

	gnb_file = open('gnb.pkl','wb')
	pickle.dump(gnb,gnb_file)
	gnb_file.close()

	tn, fp, fn, tp = confusion_matrix(y_test,gnb_prob).ravel()
	sn, sp, mcc = performance(tn, fp, fn, tp)
	gnb_acc = accuracy_score(y_test, gnb_prob)
	
	print(gnb_acc, sn, sp, mcc, auroc)

	return gnb_acc, sn, sp, mcc, auroc


def RF(x,y,x_test,y_test, data, fsort=None):
	from sklearn.ensemble import RandomForestClassifier

	rf_params = {'ccp_alpha': 0.001, 'criterion': 'entropy', 'max_depth': 110, 'max_leaf_nodes': 60, 'min_samples_leaf': 5, 'min_samples_split': 5, 'n_estimators': 63}
	rf = RandomForestClassifier(
		n_estimators=rf_params['n_estimators'],
		ccp_alpha=rf_params['ccp_alpha'],
		criterion=rf_params['criterion'],
		max_depth=rf_params['max_depth'],
		max_leaf_nodes=rf_params['max_leaf_nodes'],
		min_samples_leaf=rf_params['min_samples_leaf'],
		min_samples_split=rf_params['min_samples_split'],
		random_state=11).fit(x,y)

	if fsort is not None:
		importances=rf.feature_importances_

		feat_labels = fcol
		out = open('{0}.out'.format(each),'w')
		indices = np.argsort(importances)[::-1]
		print(indices)
		name = ['label']
		for f in range(x.shape[1]):
			name.append(feat_labels[indices[f]])
			out.write("%2d) %-*s %f\n" % (f + 1, 30, feat_labels[indices[f]], importances[indices[f]]))
		out.close()
		ndf = pd.DataFrame(columns=name)
		ndf = pd.concat([ndf, data])
		ndf.to_csv('./featureSort.csv', index=False)
	rf_prob = rf.predict(x_test)
	rf_acc = accuracy_score(y_test, rf_prob)

	auroc = roc_auc_score(y_test, rf.predict_proba(x_test)[:,1])

	rf_file = open('rf.pkl','wb')
	pickle.dump(rf,rf_file)
	rf_file.close()

	tn, fp, fn, tp = confusion_matrix(y_test, rf_prob).ravel()
	sn, sp, mcc = performance(tn, fp, fn, tp)

	print(rf_acc,sn, sp, mcc, auroc)

	return rf_acc,sn, sp, mcc, auroc


def KNN(x,y,x_test,y_test):

	from sklearn.neighbors import KNeighborsClassifier

	knn_params = {'n_neighbors': 7}

	knn = KNeighborsClassifier(n_neighbors=knn_params['n_neighbors']).fit(x,y)
	knn_prob = knn.predict(x_test)
	knn_acc = accuracy_score(y_test, knn_prob)

	auroc = roc_auc_score(y_test, knn.predict_proba(x_test)[:,1])

	knn_file = open('knn.pkl','wb')
	pickle.dump(knn,knn_file)
	knn_file.close()

	tn, fp, fn, tp = confusion_matrix(y_test, knn_prob).ravel()
	sn, sp, mcc = performance(tn, fp, fn, tp)

	print(knn_acc,sn, sp, mcc, auroc)

	return knn_acc,sn, sp, mcc, auroc

def SVM(x,y,x_test,y_test):

	from sklearn.svm import SVC
	from sklearn.model_selection import GridSearchCV

	svm_params = {
		'C': 10, 
		'gamma': 0.0001
	}
	svm = SVC(C=svm_params['C'], gamma=svm_params['gamma'], probability=True).fit(x,y)
	svm_prob = svm.predict(x_test)
	svm_acc = accuracy_score(y_test, svm_prob)
	print(svm_prob)
	auroc = roc_auc_score(y_test, svm.predict_proba(x_test)[:,1])

	svm_file = open('svm.pkl','wb')
	pickle.dump(svm,svm_file)
	svm_file.close()

	tn, fp, fn, tp = confusion_matrix(y_test, svm_prob).ravel()
	sn, sp, mcc = performance(tn, fp, fn, tp)

	print(svm_acc,sn, sp, mcc, auroc)

	return svm_acc,sn, sp, mcc, auroc

def performance(TN, FP, FN, TP):
	SN = TP / (TP + FN)  # recall
	SP = TN / (TN + FP)
	# Precision = TP / (TP + FP)
	# F1 = (2 * TP) / (2 * TP + FP + FN)
	fz = TP * TN - FP * FN
	fm = (TP + FN) * (TP + FP) * (TN + FP) * (TN + FN)
	MCC = fz / pow(fm, 0.5)
	return SN, SP, MCC

if __name__ == '__main__':
	each = 'all_feature.csv'
	data = pd.read_csv('./train/{0}'.format(each))
	test = pd.read_csv('./test/{0}'.format(each))

	fcol = data.columns[1:]
	#global fcol = data.columns[1:] #特征名称列表
	
	x = np.array(data.iloc[:,1:])
	y = np.array(data['label'])
	x_test = np.array(test.iloc[:,1:])
	y_test = np.array(test['label'])

	print(x.shape,y.shape)

	# kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=11)

	#gnb_acc, gsn, gsp, gmcc, gauroc = gnb(x,y,x_test,y_test)
	rf_acc, rsn, rsp, rmcc, rauroc = RF(x,y,x_test,y_test,data, fsort=1)
	#knn_acc, ksn, ksp, kmcc, kauroc = KNN(x,y,x_test,y_test)
	#svm_acc, ssn, ssp, smcc, sauroc = SVM(x,y,x_test,y_test)

	#out = open(f'/home/yhyang/protein/task/out/result_test.out','w')
	#out.write(str(gnb_acc)+'\t'+str(gsn)+'\t'+str(gsp)+'\t'+str(gmcc)+'\t'+str(gauroc)+'\n')
	#out.write(str(rf_acc)+'\t'+str(rsn)+'\t'+str(rsp)+'\t'+str(rmcc)+'\t'+str(rauroc)+'\n')
	#out.write(str(knn_acc)+'\t'+str(ksn)+'\t'+str(ksp)+'\t'+str(kmcc)+'\t'+str(kauroc)+'\n')
	#out.write(str(svm_acc)+'\t'+str(ssn)+'\t'+str(ssp)+'\t'+str(smcc)+'\t'+str(sauroc)+'\n')
	#out.close()