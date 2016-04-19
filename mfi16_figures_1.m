% part 1 of /Users/alex/Documents/research/proton/code/calibration/motion/mfi16_figures.m
% generates tbl:regression for the MFI 2016 paper

% calibrates using the datasets taken on 2/26/16 and 2/23/16
% loads the datasets taken on 2/19/16 and 3/10/16 and does machine learning

addpath(genpath('RANSAC-Toolbox'))
addpath('libsvm/matlab')

DATADIR = '../../nri/data';

