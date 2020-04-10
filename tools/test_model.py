#!/usr/local/bin/python
# Author: Eli Draizen
# Date: 16-3-2014
# File: classify.py

import matplotlib
matplotlib.use('Agg')
 
#Standard Libraries
import argparse
import json
import os
import logging

#BioPython
from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq

#Required Libraries
import numpy as np
import seaborn as sns

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, matthews_corrcoef
from matplotlib.backends.backend_pdf import PdfPages

from scipy.interpolate import interp1d

sns.set(style="white", context="talk")

#Thresholds decided by manualy looking at the headers in the hmmsearch output. Useful to compare.
original_thresholds = {
    "H2A": {
        "H2A.X":275, 
        "H2A.Z":230, 
        "H2A.B":115, 
        "H2A.L":145, 
        "H2A.M":80, 
        "macroH2A":270, 
        "canonical_H2A":265
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

def test_model(model_name, save_dir, postive_file, negative_file, measure="SPC", measure_threshold=0.95):
    """Test the model by calcuating

    Returns:
    A dictionary with containg the AUCROC and Threshold. An image is also saved
    with the ROC and score histograms. 
    """
    logging.info(model_name)
    postive_scores = get_model_scores(postive_file)
    negative_scores = get_model_scores(negative_file)
    all_scores = postive_scores+negative_scores
    # print all_scores

    if len(negative_scores) == 0:
        return {"roc_auc":0, "threshold":min(postive_scores)}

    y_true = [1]*len(postive_scores) + [0]*len(negative_scores)
    y_score = np.array(all_scores)

    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)

    best_threshold, thresholds, values = calcualte_threshold(
        postive_scores, 
        negative_scores, 
        measure=measure,
        measure_threshold=measure_threshold, 
        thresholds=reversed(thresholds))


    pp = PdfPages(os.path.join(save_dir, "{}_model_evaluation.pdf".format(model_name)))

    sns.set(style="darkgrid")
    f, axes = plt.subplots(3)
    trans = f.transFigure.inverted()
    colors = sns.color_palette("Set2", 7)

    sns.kdeplot(np.array(postive_scores), shade=True, color=sns.xkcd_rgb["denim blue"], label="Scores for postive examples", ax=axes[0])
    sns.kdeplot(np.array(negative_scores), shade=True, color=sns.xkcd_rgb["pale red"], label="Scores for negative examples",  ax=axes[0])
    axes[0].set_xlabel("Bit score")
    axes[0].set_ylabel("Density")
    axes[0].legend(loc="upper left")
    #axes[0].set_title("Kernel Density of Scores")
    axes[1].set_xlim([0, 1.0])
    axes[1].set_ylim([0.0, 1.05])

    
    axes[1].plot(fpr,tpr, color=colors[0], lw=3., label="ROC (AUC: {})".format(roc_auc))
    axes[1].set_xlabel("False Positive Rate")
    axes[1].set_ylabel("True Positive Rate")
    axes[1].legend(loc="lower right")
    axes[1].set_xlim([-0.05, 1.0])
    axes[1].set_ylim([0.0, 1.05])
    #axes[1].set_title("ROC")
 
    for i, (measure, values) in enumerate(values.items()):
        label = "SPC: (>={})".format(best_threshold) if measure=="SPC" else measure
        axes[2].plot(list(thresholds), values, label=label, linewidth=2, color=colors[i])
    axes[2].axvline(best_threshold)

    axes[2].legend()
    #axes[2].set_title("Coosing Cutoff")
    axes[2].set_ylabel("Rate")
    axes[2].set_xlabel("Threshold")

    f.suptitle("{} Model Evaluation".format(model_name), fontsize=20)

    pp.savefig()
    pp.close()

    return {"roc_auc":roc_auc, "threshold":best_threshold}

def calcualte_threshold(positives, negatives, measure="SPC", measure_threshold=0.95, thresholds=None, attempt=0):
    """Plot the TPR the FPR vs threshold values
 
    Input:
    postives - list of scores of postive runs
    negatives - list of scores of negative runs
    measure - choose coffectiong by 95% Specificity ("SPC"), or matthews_corrcoef ("MCC")
    """
    assert measure in ["TPR", "FPR", "SPC", "MCC", "PPV", "FDR", "ACC"]
    y_true = [1]*len(positives)+[0]*len(negatives)
    values = {name:[] for name in ["TPR", "FPR", "SPC", "MCC", "PPV", "FDR", "ACC"]}
    saveThreshold = None
    saveValue = 1.0
    thresholds = list(sorted(thresholds or [i/10. for i in range(1,10000)]))

    for threshold in thresholds:
        TN = sum([1 for score in negatives if score < threshold])
        FP = sum([1 for score in negatives if score >= threshold])
        TP = sum([1 for score in positives if score >= threshold])
        FN = sum([1 for score in positives if score < threshold])

        values["FPR"].append(float(FP)/(FP+TN))
        values["TPR"].append(float(TP)/(TP+FN))
        values["SPC"].append(float(TN)/(FP+TN))

        y_pred = [int(score >= threshold) for scores in (positives, negatives) for score in scores]
        values["MCC"].append(matthews_corrcoef(y_true, y_pred))
        values["PPV"].append(float(TP)/(TP+FP) if TP+FP>0 else 0.0)
        values["FDR"].append(float(FP)/(TP+FP) if TP+FP>0 else 0.0)
        values["ACC"].append(float(TP+TN)/(len(positives)+len(negatives)))
        
    specificity_curve_inverse = interp1d(values[measure], thresholds)
    saveThreshold = specificity_curve_inverse(measure_threshold) #modified by Alexey 0.95 => measure_threshold
                
    return saveThreshold, thresholds, values


