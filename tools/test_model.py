#!/usr/local/bin/python
# Author: Eli Draizen
# Date: 16-3-2014
# File: classify.py
 
#Standard Libraries
import argparse
import json
import os

#Required Libraries
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, matthews_corrcoef
from matplotlib.backends.backend_pdf import PdfPages

sns.set(style="white", context="talk")

#Thresholds decided by manualyy looking at the headers in the hmmsearchout. Useful to compare.
original_thresholds = {
    "H2A": {
        "H2A.X":275, 
        "H2A.Z":230, 
        "H2A.B":115, 
        "H2A.L":145, 
        "H2A.M":80, 
        "macroH2A":270, 
        "canonicalH2A":265
    }
}

def get_model_scores(model_output):
    """Get the bit score for each hit/domain in a hmmersearch result

    Parameters:
    -----------
    model_output : str or File-like object
        Path to hmmersearch output file

    Return:
    -------
    A list of all bitscores
    """
    return [hsp.bitscore for query in SearchIO.parse(model_output, "hmmer3-text") \
        for hit in query for hsp in hit] 

def test_model(model_name, postive_file, negative_file, measure="SPC"):
    """Test the model by calcuating

    Returns:
    A dictionary with containg the AUCROC and Threshold. An image is also saved
    with the ROC and score histograms. 
    """
    postive_scores = get_model_scores(postive_file)
    negative_scores = get_model_scores(negative_file)
    all_scores = postive_scores+negative_scores

    y_true = [1]*len(postive_scores) + [0]*len(negative_scores)
    y_score = np.array(all_scores)

    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)

    threshold, values = calcualte_threshold(postive_scores, negative_scores, measure=measure, thresholds=reversed(thresholds))

    pp = PdfPages("{}_model_evaluation.pdf".format(model_name))

    sns.set(style="darkgrid")
    f, axes = plt.subplots(1, 3)
    trans = f.transFigure.inverted()
    colors = sns.color_palette("Set2", 7)

    sns.kdeplot(postive_scores, shade=True, color=sns.xkcd_rgb["denim blue"], label="Scores for postive examples", ax=axes[i, 0])
    sns.kdeplot(negative_scores, shade=True, color=sns.xkcd_rgb["pale red"], label="Scores for negative examples",  ax=axes[i, 0])
    if i==len(scores)-1: axes[1,0].set_xlabel("Bit score")
    if i==2: axes[i,0].set_ylabel("Density")
    if i==0: 
        axes[i,0].legend(loc="upper left")
        axes[i,0].set_title("Kernel Density of Scores")
    else:
        axes[i,0].legend_.remove()
    axes[i,1].set_xlim([0, 1.0])
    axes[i,1].set_ylim([0.0, 1.05])

    
    axes[1,1].plot(fpr,tpr, color=colors[0], lw=3., label="ROC (AUC: {})".format(roc_auc))
    if i==len(scores)-1: axes[1,1].set_xlabel("False Positive Rate")
    if i==2: axes[1,1].set_ylabel("True Positive Rate")
    axes[1,1].legend(loc="lower right")
    axes[1,1].set_xlim([-0.05, 1.0])
    axes[1,1].set_ylim([0.0, 1.05])
    if i==0: 
        axes[1,1].set_title("ROC")

    ax.axvline(theshold)
 
    for i, (measure, values) in enumerate(values.iteritems()):
        axes[i,2].plot(thresholds, values, label=measure, linewidth=2, color=colors[i])

    if i==0: 
        axes[i,2].legend()
        axes[i,2].set_title("Coosing Cutoff")
    if i==2: axes[i,2].set_ylabel("Rate")
    if i==len(scores)-1:axes[i,2].set_xlabel("Threshold")

    fig.suptitle("{} Model Evaluation".format(model_name), fontsize=20)

    pp.savefig()
    pp.close()

    return {"roc_auc":roc_auc, "threshold":threshold}

def calcualte_threshold(positives, negatives, measure="SPC", measure_threshold=0.95, thresholds=None):
    """Plot the TPR the FPR vs threshold values
 
    Input:
    postives - list of scores of postive runs
    negatives - list of scores of negative runs
    measure - choose coffectiong by 95% Specificity ("SPC"), or matthews_corrcoef ("MCC")
    """
    assert measure in ["TPR", "FPR", "SPC", "MCC", "PPV", "NPV", "FDR", "ACC"]
    y_true = [1]*len(positives)+[0]*len(negatives)
    values = {name:[] for name in ["TPR", "FPR", "SPC", "MCC", "PPV", "NPV", "FDR", "ACC"]}
    saveTheshold = None
    saveValue = 1.0
    thresholds = list(thresholds or map(lambda i: i/10., xrange(1,10000)))
    for threshold in thresholds:
        #print threshold
        TN = sum([1 for score in negatives if score < threshold])
        FP = sum([1 for score in negatives if score >= threshold])
        TP = sum([1 for score in positives if score >= threshold])
        FN = sum([1 for score in positives if score < threshold])
        #print "FP", FP
        values["FPR"].append(float(FP)/(FP+TN))
        values["TPR"].append(float(TP)/(TP+FN))
        values["SPC"].append(float(TN)/(FP+TN))

        y_pred = [int(score >= threshold) for scores in (positives, negatives) for score in scores]
        values["MCC"].append(matthews_corrcoef(y_true, y_pred))
        values["PPV"].append(float(TP)/(TP+FP) if TP+FP>0 else 0.0)
        values["MCC"].append(float(TN)/(TN+FN) if TN+FN>0 else 0.0)
        values["FDR"].append(float(FP)/(TP+FP) if TP+FP>0 else 0.0)
        values["ACC"].append(float(TP+TN)/(len(positives)+len(negatives)))
        
        if round(values[measure][-1]*20)/20 == measure_threshold and \
            abs(measure_threshold-values[measure][-1])<abs(measure_threshold-saveValue):
            saveTheshold = threshold
            saveValue = values[measure][-1]

    return saveTheshold, values


