#-*- coding: utf-8 -*-

"""
Created on Mon Jan 6 09:20:06 2020 

@author: XiongYao
"""

import os, sys, logging
import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import pickle
from userConfig import AAS3Drf
from sklearn.preprocessing import Imputer, StandardScaler


def _FillNone(method, train, test):
    if method == "mean":
        imp = Imputer(missing_values="NaN", strategy="mean", axis=0)
        imp.fit(train)
        train_new = imp.transform(train)
        test_new = imp.transform(test)
    return train_new, test_new
        

def _Scaler(method, train, test):
    if method == "standardization":
        scaler = StandardScaler()
        train_new = scaler.fit_transform(train)
        test_new = scaler.transform(test)
    return train_new, test_new


def AAS3DrfPred(infile, outfile, work_path):
    logging.info("Predicting...")
    
    trainfile = os.path.join(AAS3Drf, "PredMod", "preData", "TrainDataset")
    testfile = os.path.join(AAS3Drf, "Data", "USER_DATA.tsv")

    train = pd.read_table(trainfile, header=0, na_values="-", sep="\t")
    title = train.columns.values.tolist()
    fea_order = title[title.index("Functional_Effect")+1:]
    x_train = train.ix[:,fea_order]
    
    test = pd.read_table(testfile, header=0, na_values="-", sep="\t")
    x_test = test.ix[:,fea_order]

    imp_train, imp_test = _FillNone("mean", x_train, x_test)
    scale_train, scale_test = _Scaler("standardization", imp_train, imp_test)

    model_file = os.path.join(AAS3Drf, "PredMod", "RF_Model")
    open_model = open(model_file, "r")
    rf_model = pickle.load(open_model)

    y_pred, y_pred_prob = rf_model.predict(scale_test), rf_model.predict_proba(scale_test)
    test.insert(4,"Result",y_pred)
    test.insert(4,"Score",y_pred_prob[:,1])
    test["Result"] = test["Result"].replace((1,0),("Disease","Neutral"))

    if outfile is None:
        name = os.path.splitext(infile)[0]
        outfile = name+"_AAS3D-RF.out"
    
    test.to_csv(os.path.join(work_path, outfile), sep="\t", na_rep="-", index=False)
