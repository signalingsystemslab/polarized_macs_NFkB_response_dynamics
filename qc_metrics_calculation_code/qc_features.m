function [] = qc_features(ids, expt_names)
% ---------------------------------------------------------------------
% calculates the signaling codon values (Adelaja et al. 2021) that were
% used for quality control for a collection of experiments
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

tot_data = table();
num_traj = zeros(length(ids), 1);
for i = 1:length(ids)
    [qc_metrics] = nfkbQCmetrics(ids(i));
    qc_metrics = struct2table(qc_metrics);
    num_traj(i) = size(qc_metrics, 1);
    tot_data = [tot_data; qc_metrics];
end

nanzscore = @(x)(x-nanmean(x, 1))./nanstd(x, 0, 1);

Duration = mean([nanzscore(tot_data.dur_t), nanzscore(tot_data.n_pks)], 2);
EarlyVsLate = -nanzscore(tot_data.eVl);
OscVsNon = nanzscore(tot_data.oVn);
PeakAmplitude = mean([nanzscore(tot_data.max_val), nanzscore(tot_data.pk2pk),nanzscore(tot_data.pk1_amp)], 2);
Speed = mean([nanzscore(tot_data.max_pk1_spd), -nanzscore(tot_data.pk1_t), nanzscore(tot_data.deriv2)], 2);
TotalActivity = nanzscore(tot_data.tot_act);

final_data = table(Duration, EarlyVsLate, OscVsNon, PeakAmplitude, Speed, TotalActivity);

row = 1;
for i = 1:length(ids)
    sc = final_data(row:row+num_traj(i)-1, :);
    writetable(sc, strcat(expt_names(i), ".xlsx"))
    row = row + num_traj(i);
end
end