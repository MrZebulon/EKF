clc;
close all;
clear all;
%% MEKF instance initialization
model = SimplifiedModel();
instance = MEKF(zeros(10, 1), diag(1e-9 * ones(1, 10)), model);