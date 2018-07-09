# BORKOV
This is my undergraduate thesis. Automatic Segmentation of Multiple Sclerosis Lesions using KNN and MARKOV Random Fields

# Introduction
BORKOV is an AI application that automates the segmentation process of multiple sclerosis lesions in MR images using KNN and MRF.
The results are probabilistic. 
We require a certain level of confidence input on how much you want to restrict the probability of each node to be to be considered a lesion.
We had a lot of constraints working in this such as we don't have access to a wide array of datasets and specifically time. 
We came to a conclusion that this unsupervised method is for enhancement due to our findings having high false positive ratios but it is lessened with higher confidence.

# Thesis Paper
It's included in the repo.

# How To Use
Kindly read our thesis paper where we found our datasets. You will need them. You will get everything there as well as what to run and how to use.

# Technology
We mainly used SimpleITK for the MR images support.
As well as, scikit for KNN methods.
MRF was done from scratch.
Tktinker for the simple GUI.
