function R=evaluate_kriging_model_ccdbr(x,omega,varargin)
load kriging_model_ccdbr.mat
Y=predictor(x',dmodel);
R{1}=20*log10(max(Y,0.001));
