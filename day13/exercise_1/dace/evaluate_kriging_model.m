function R=evaluate_kriging_model(x,omega,varargin)
load kriging_model.mat
Y=predictor(x',dmodel);
R{1}=Y;
