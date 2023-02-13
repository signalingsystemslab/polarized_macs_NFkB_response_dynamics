function [] = calculate_features(ids, expt_names)
% ---------------------------------------------------------------------
% calculates all feature values for a collection of experiments used for
% feature based analysis
%
% INPUT
%   ids = list of n numerical id numbers (from spreadsheet) corresponding 
%   to experiments of interest [###, ###, ...] or names of excel
%   spreadsheets with time series data
%
%   expt_names = list of n strings naming each experiment in ids ["...",
%   "...", ...]
% OUTPUT
% generates excel spreadsheet for each experiment in current directory
% ---------------------------------------------------------------------
for i = 1:length(ids)
    metrics = nfkbmetrics(ids(i));
    metrics = struct2table(metrics);
    metrics = splitvars(metrics);
    writetable(metrics, strcat(expt_names(i), "_metrics.xlsx"))
end
end