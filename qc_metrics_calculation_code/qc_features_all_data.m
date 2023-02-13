function [] = qc_features_all_data()
% ---------------------------------------------------------------------
% calculates the signaling codon values (Adelaja et al. 2021)
% for a collection of experiments that were used for quality control

% OUTPUT
% generates excel spreadsheet for each experiment in current directory
% ---------------------------------------------------------------------
path = "C:\Users\apeks\Box\Hoffmann\Spring 2021\time_series_data_files\";
stimuli = ["TNF", "R84", "PIC", "P3K", "CpG", "FLA", "FSL", "LPS", "UST"];
polarization = ["", "ib", "ig", "i0", "i3", "i4"];
rep = ["1", "2"];
tot_data = table();
num_traj = zeros(((length(stimuli)-1)*length(polarization)*2)+2, 1);
count = 1;
for i = stimuli
    if i == "UST"
        pol_queque = [""];
    else
        pol_queque = polarization;
    end
    for j = pol_queque
        for k = rep
            i+j+k
            [qc_metrics] = nfkbQCmetrics(path+i+j+k);
            qc_metrics = struct2table(qc_metrics);
            size(qc_metrics)
            tot_data = [tot_data; qc_metrics];
            num_traj(count) = size(qc_metrics, 1);
            count = count + 1;
        end
    end
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
count = 1;
for i = stimuli
    if i == "UST"
        pol_queque = [""];
    else
        pol_queque = polarization;
    end
    for j = pol_queque
        for k = rep
            cond_metrics = final_data(row:row+num_traj(count)-1, :);
            writetable(cond_metrics, i + j + k +"_QCmetrics.xlsx")
            row = row + num_traj(count);
            count = count + 1;
        end
    end
end
end