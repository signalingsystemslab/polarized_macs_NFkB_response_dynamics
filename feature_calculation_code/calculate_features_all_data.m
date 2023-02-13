function [] = calculate_features_all_data()
% ---------------------------------------------------------------------
% calculates all feature values for a collection of experiments used for
% feature based analysis

% OUTPUT
% generates excel spreadsheet for each experiment in current directory
% ---------------------------------------------------------------------
path = "C:\Users\apeks\Box\Hoffmann\Spring 2021\time_series_data_files\";
stimuli = ["TNF", "R84", "PIC", "P3K", "CpG", "FLA", "FSL", "LPS", "UST"];
polarization = ["", "ib", "ig", "i0", "i3", "i4"];
rep = ["1", "2"];
for i = stimuli
    if i == "UST"
        pol_queque = [""];
    else
        pol_queque = polarization;
    end
    for j = pol_queque
        for k = rep
            i+j+k
            metrics = nfkbmetrics(path+i+j+k);
            metrics = struct2table(metrics);
            metrics = splitvars(metrics);
            writetable(metrics, i + j + k + "_metrics.xlsx")
        end
    end
end
end