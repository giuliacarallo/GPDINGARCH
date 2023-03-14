# GPDINGARCH

This contains codes for the paper "Generalized Poisson Difference Autoregressive Processes" by Carallo, G., Casarin, R. and Robert, C. P.
available at: https://arxiv.org/abs/2002.04470

AdaptiveGPDINGARCH.m and CyberEstimation.m are the main code files for the estimation of two applications: Schiphol and Cyber threats datasets. Datasets are available in the .mat files Schiphol. mat and CyberAttacks2017_2018.mat

SeqForecastAdap.m is the main code for the forecasting exercise, can be adapted for all models and applications studied in the paper. The other files are the functions necessary for the code to run:

alpha0Sampl.m philambdaAdap.m varphiAdap.m xySampl2New.m are the sampling functions;

gpbranching.m is the function to sample from the Generalized Poisson distibution;

lkl.m is the function for evaluating the likelihood of the Generalized Poisson Difference distribution;

mut.m is the function for generating the garch process;

gibbsadap.m and mut_forec2.m are the functions for approximating forecasting distributions;

 
