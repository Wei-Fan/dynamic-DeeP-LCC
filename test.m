clc; close all; clear;
addpath('_fcn');
addpath('_data');

% Random setup for OVM
load(['_data/hdv_ovm_random_largescale.mat']);

% Loop through the fields of the input struct
alpha = hdv_parameter.alpha;
beta = hdv_parameter.beta;
s_st = hdv_parameter.s_st;
s_go = hdv_parameter.s_go;
s_star = hdv_parameter.s_star;