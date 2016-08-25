% part 1 of /Users/alex/Documents/research/proton/code/calibration/motion/icra17_figures.m
% generates figures for the ICRA 2017 paper

% calibrates using the datasets taken on 2/26/16, 7/26/16, and 8/11/16
% loads the dataset taken on 8/15/16 and does machine learning

addpath(genpath('RANSAC-Toolbox'))
addpath('libsvm/matlab')

DATADIR = '/Volumes/shared/Projects/Proton Pack/Data';

