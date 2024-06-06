close all; clear all; clc
warning off;
addpath(genpath('function'));
load('bbcsport.mat');
lambda =0.01;
beta = 0.5;  
alpha = 2;    
idx = 1;  

[S,restL]=TMVC(X,lambda, beta, alpha ,gt);
idx = idx + 1;     
disp(restL)
