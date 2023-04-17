# O-GlyThr is a tool for prediction of O-glycosites linked to the threonine residues in Homo sapiens, which identifies a “T” as a positive O-glycosite if its probability is 0.5 or higher. 

## You can use our training data to train a new model according to the following steps:
1.	Create two file folders named “train” and “test”.
2.	Run “python all_feature_extraciton.py pos_train.fa neg_train.fa”, and copy the result “all_feature.csv” into the “train” folder. 
3.	Run “python all_feature_extraciton.py pos_test.fa neg_test.fa” , and copy the result “all_feature.csv” into the “test” folder.
4.	Run “python featureSort.py”, and you can obtain the result file: “featureSort.csv”.
5.	Run “python IFS.py featureSort.csv”

## You can also perform O-glycosites prediction with our provided model (rf.pkl):
Run “python test_pred.py test.fa output.txt”. 
All the needed files are in the file folder named “testOGlyThr”. 

Note: the code “all_feature_extraciton.py” was learned from the following work. Lv H, Zhang Y, Wang JS, Yuan SS, Sun ZJ, Dao FY, Guan ZX, Lin H, Deng KJ. iRice-MS: An integrated XGBoost model for detecting multitype post-translational modification sites in rice. Brief Bioinform. 2022 Jan 17;23(1): bbab486.
