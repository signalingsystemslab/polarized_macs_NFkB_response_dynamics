function [metrics, aux, graph, info, measure] = nfkbmetrics(expt)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Adpated from Adelaja A, Taylor B, Sheu KM, Liu Y, Luecke S, Hoffmann A. 
% Six distinct NFKB signaling codons convey discrete information to distinguish stimuli 
% and enable appropriate macrophage responses. Immunity. 2021;54(5):916-930.e7. 
% doi:10.1016/j.immuni.2021.04.011
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% metrics = nfkbmetrics(id)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NFKBMETRICS uses the filter_nfkb function to filter and preprocess NFkB trajectories,
% then calculates many related metrics regarding activation
% 
% INPUTS (required):
% expt            experiment ID (from Google Spreadsheet specified in
% "locations.mat") or excel filename of time-series data
%
% OUTPUT: 
% metrics       structure with all output metric fields
% aux           Extra data (e.g. fourier information (FFT, power, frequencies), thresholds used in envelope/duration)
% graph, info, measure         main structure outputs from filter_nfkb if computed
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%% INITIALIZATION. Load and process data.
MinLifetime = 97; 
endFrame = 100;
Delay = 3;

if isnumeric(expt)
    [graph, info, measure] = filter_nfkb(expt,'MinLifetime',MinLifetime,'ConvectionShift',1,...
        'Delay',Delay, 'TrimFrame', endFrame);
    graph.var = graph.var(:,1:min(endFrame, size(graph.var,2)));
    graph.var_nfkb_no_base_ded = graph.var_nfkb_no_base_ded(:,1:min(endFrame, size(graph.var_nfkb_no_base_ded,2)));
    graph.t = graph.t(1:min(size(graph.var,2),endFrame));
    baseline = info.baseline;
    % 1) basic time series. Interpolate over "normal" interval (12 frames per hr) if required
    % - Taylor B
    t = min(graph.t):1/12:max(graph.t);
    if length(t)~=length(graph.t)
        metrics.time_series = nan(size(graph.var,1),length(t));
        for i = 1:size(graph.var,1)
            metrics.time_series(i,:) = interp1(graph.t,graph.var(i,:),t);
        end
        time_series_no_base_ded = nan(size(graph.var_nfkb_no_base_ded,1),length(t));
        for i = 1:size(graph.var_nfkb_no_base_var_ded,1)
            time_series_no_base_ded(i,:) = interp1(graph.t,graph.var_nfkb_no_base_ded(i,:),t);
        end
    else
        metrics.time_series = graph.var;
        time_series_no_base_ded = graph.var_nfkb_no_base_ded;
    end
    FramesPerHour = 12;
    baseline_stdv = nanstd(metrics.time_series(:,1:Delay),0,2);
    metrics.time_series = metrics.time_series(:, Delay:end);
    time_series_no_base_ded = time_series_no_base_ded(:, Delay:end);
    t= t(Delay:end);
else
    metrics.time_series = xlsread(expt, 'time_series');
    time_series_no_base_ded = xlsread(expt, 'time_series_no_bsl_ded');
    baseline = xlsread(expt, 'baseline');
    baseline_stdv = xlsread(expt, 'bsl_stdv');
    t= xlsread(expt, 'time_h');
    FramesPerHour = 12;
end

%% BASIC METRICS
% - Taylor B
% 2) integrated activity, trapezoidal approximation
integrals = nan(size(metrics.time_series));
nan_removed = metrics.time_series;
nan_removed(isnan(nan_removed)) = 0;
for i = 1:size(integrals,1)
    integrals(i,:) = cumtrapz(t,nan_removed(i,:));
end

% 3) differentiated activity - use central finite difference
smoothed = medfilt1(metrics.time_series,3,[],2);
derivatives = (smoothed(:,3:end) - smoothed(:,1:end-2))/(1/6);
metrics.derivatives = derivatives(:, 1:25);%keep only first two hours
%% TRIM to a common length
endFrame = endFrame-Delay+1;
try
    metrics.time_series = metrics.time_series(:,1:endFrame);
    integrals = integrals(:,1:endFrame);
    derivatives = derivatives(:,1:(endFrame-2));
    smoothed = smoothed(:,1:endFrame);
    t= t(1:endFrame);
catch
    disp(['Note: vectors too short to cap @ ',num2str(endFrame),' frames'])
end

%% INTEGRAL WINDOWS
% Integrals within one-hour windows (0-1, 1-2, 2-3) and three hour windows
% (0-3, 1-4, etc) of activity - adapted from Taylor B
max_hr = floor(max(t));
metrics.intwin1 = nan(size(metrics.time_series,1),max_hr);
metrics.intwin3 = nan(size(metrics.time_series,1),max_hr-2);
for i = 1:(max_hr)
    win = t>=(i-1) & t<=(i);
    metrics.intwin1(:,i) = trapz(t(win),nan_removed(:,win),2);
    if i<= (max_hr-2)
        win = t>=(i-1) & t<=(i+2);
        metrics.intwin3(:,i) = trapz(t(win),nan_removed(:,win),2);
    end
end

%added difference between hour 0-1 integral and hour 1-2 integral
metrics.phase_diff1 = metrics.intwin1(:, 1)-metrics.intwin1(:, 2);

%added difference between hour 0-3 integral and hour 3-6 integral
metrics.phase_diff3 = metrics.intwin3(:, 1)-metrics.intwin3(:, 4);

%added half hour windows as well
max_HalfHr = floor(max(t*2));
metrics.intwin0_5 = nan(size(metrics.time_series,1),max_HalfHr);
for i = 1:(max_HalfHr)
    win = t>=((i-1)/2) & t<=(i/2);
    metrics.intwin0_5(:,i) = trapz(t(win),nan_removed(:,win),2);
end

%% MAX/MIN METRICS
% - adapted from Taylor B
metrics.max_amplitude = nanmax(metrics.time_series(:,1:end),[],2);
metrics.min_amplitude = nanmin(metrics.time_series(:,1:end),[],2);
metrics.range = metrics.max_amplitude-metrics.min_amplitude;

%adding calculations for half max in early activity
[metrics.maxAmp_early, metrics.time2Max, metrics.timeUp2halfMax, metrics.timeDown2halfMax] = early_activity(metrics.time_series);          

function [max_amp, time2Max, timeUp2halfMax, timeDown2halfMax] = early_activity(time_series)
    [max_amp, time2Max] = nanmax(time_series(:,2:25),[],2);
    time2Max = time2Max + 1; %add 1 to account for starting consideration after zero timepoint
    half_max = ((max_amp-time_series(:, 1))/2) + time_series(:, 1);
    cross_halfMax = time_series-half_max>=0;
    [~, timeUp2halfMax] = nanmax(cross_halfMax, [], 2);
    timeDown2halfMax = nan(size(time2Max));
    for index=1:length(time2Max)
        cross_halfMax = time_series(index, time2Max(index):end)<=half_max(index);
        [~, timeDown2halfMax(index)] = nanmax(cross_halfMax);
    end
    timeUp2halfMax = (timeUp2halfMax - 1)/12;
    timeDown2halfMax = (timeDown2halfMax - 1)/12;
    time2Max = (time2Max - 1)/12;
    % if max_amp below time zero value --> timeUp/down will be zero
end

% METRICS OF DURATION RELATIVE TO MAX AMP
% added aditional relative duration metric
aux.peak_thresholds = metrics.max_amplitude*[0.5, 0.7, 0.9];

% Number of frames above a given threshold
metrics.peak_duration = zeros(size(metrics.time_series,1),size(aux.peak_thresholds, 2));
for i = 1:size(metrics.time_series,1)
    for j = 1:size(aux.peak_thresholds, 2)
        metrics.peak_duration(i,j) = nansum(smoothed(i, 1:end)>aux.peak_thresholds(i, j))/FramesPerHour;
    end
end

%% ACTIVITY METRICS
% Determine which cells are responders and compute off-times - Luecke S
Wliml = 1; %first/lower time point of window to check for activity
Wlimu = 48; %last/upper time point of window to check for activity, ie check in the first 4 hours after stimulation
blockLengthThresh = 5; %number of consecutive frames cell needs to pass activity threshold to be considered a responder
OnThresh = 3; %value set as the standard deviation theshold above baseline reqd to be considered a "responder" (also called sigmaThresh)
smoothed_by_sigma = smoothed./baseline_stdv;
[metrics.responder_index, responders_fraction, metrics.off_times] = get_activity_metrics(smoothed_by_sigma, Wliml, Wlimu, OnThresh, blockLengthThresh);
metrics.off_times = metrics.off_times/FramesPerHour;
metrics.off_times(metrics.off_times<0) = 0;

%% METRICS OF DURATION
% Envelope width: maximum consecutive time above a threshold (envelope must
% begin within 1st 6 hrs) - Taylor B
smoothed2 = medfilt1(metrics.time_series,5,[],2);
lowerThresh = 0;
upperThresh = 2.25; %down from 7.1
aux.thresholds = linspace(lowerThresh,upperThresh,4); %down from 25 thresholds
metrics.envelope = zeros(size(metrics.time_series,1),4);
for j = 1:length(aux.thresholds)
    thresholded = smoothed2(:, 1:end)>aux.thresholds(j);
    for i = 1:size(thresholded,1)
        curr = 1;
        idx_start = 1;
        while (curr<size(thresholded,2)) && (idx_start< (6*FramesPerHour))
            idx_start = find(thresholded(i,curr:end)==1,1,'first')+curr-1;
            if ~isempty(idx_start)
                idx_stop = find(thresholded(i,idx_start:end)==0,1,'first')+idx_start-1;
                if isempty(idx_stop)
                    idx_stop = find(~isnan(thresholded(i,:)),1,'last');
                end
                if (idx_stop-idx_start) > metrics.envelope(i,j)
                    metrics.envelope(i,j) = (idx_stop-idx_start);
                end
                curr = idx_stop;
            else
                break
            end
        end
    end
end
metrics.envelope = metrics.envelope/FramesPerHour;

% Number of frames above a given threshold
metrics.duration = zeros(size(metrics.time_series,1),length(aux.thresholds));
for i = 1:length(aux.thresholds)
    metrics.duration(:,i) = nansum(smoothed(:, 1:end)>aux.thresholds(i),2)/FramesPerHour;
end

%% METRICS OF OSCILLATION
% Calculate fourier distribution (via FFT) & power - Taylor B
Fs = 1/300;
depth = max(metrics.off_times)*FramesPerHour;
NFFT = 2^nextpow2(depth); % Next power of 2 from chosen depth
aux.fft = zeros(size(metrics.time_series,1),NFFT/2+1);
aux.freq = Fs/2*linspace(0,1,NFFT/2+1);
aux.power = zeros(size(aux.fft));
off_pad = 12; % Signal time added to trajectory in  FFT calculation (keeps transients from being recorded as osc.)

for i = 1:size(metrics.time_series,1)
    if(metrics.off_times(i)>0)
        y = metrics.time_series(i,1:(depth));
        off_frame = min([length(y), metrics.off_times(i)*FramesPerHour+1+off_pad]); % (Pad w/ 1 extra hr of content)
        y(off_frame:end) = nan;
        y(isnan(y)) = [];
        y = y-nanmean(y);
        if ~isempty(y)
            Y = fft(y,NFFT)/length(y);
            aux.fft(i,:) = abs(Y(1:NFFT/2+1));
            aux.power(i,:) = abs(Y(1:NFFT/2+1).^2);
        end
    end
end

% Find the point of peak (secondary) power
metrics.peakfreq = nan(size(aux.power,1),1);
for i =1:size(metrics.time_series,1)
    [pks,locs] = globalpeaks(aux.power(i,:),2);
    % Ensure we're not getting a totally spurious peak
    if min(pks) < (0.1*max(pks))
        locs(pks==min(pks)) = [];
    end
    if length(locs)>1
        idx = max(locs(1:2));
        metrics.peakfreq(i) = 3600*aux.freq(idx);
    elseif ~isempty(locs)
         metrics.peakfreq(i) = 3600*aux.freq(max([locs,3]));
    else
        metrics.peakfreq(i) = 3600*aux.freq(1);
    end
end

%% METRICS OF AMPLITUDE AND TIMING
% modified from Taylor B to extract only positive value peaks
% 1st peak only time/amplitude/prominence/width

pk_feats = {'pk1_amp', 'pk1_time', 'pk1_width', 'pk1_prom'};
%pk_feats = {'pk1_amp', 'pk1_time', 'pk1_width', 'pk1_prom', 'pk1_height'};

for j=1:length(pk_feats)
    metrics.(pk_feats{j}) = nan(size(time_series_no_base_ded,1),1);
end
for j = 1:size(metrics.pk1_time,1)    
    [pks, locs, width, prom, heights] = global_positive_peaks(time_series_no_base_ded(j,1:min([90,MinLifetime])),baseline(j),5);
    % Supress any peaks that are within 6 frames of each other.
    [locs, order] = sort(locs,'ascend');
    pks = pks(order); width = width(order); prom = prom(order); heights = heights(order);
    while min(diff(locs))<6
        tmplst = find(diff(locs)==min(diff(locs)),1,'first');
        tmplst = tmplst + (pks(tmplst)>=pks(tmplst+1));
        pks(tmplst) = []; locs(tmplst) = []; width(tmplst) = []; prom(tmplst) = []; heights(tmplst) = [];
    end
    %remove any peaks occuring prior to fifteen minutes or after 4 hours
    pks(locs< 3 | locs>48) = [];
    width(locs< 3 | locs>48) = [];
    prom(locs< 3 | locs>48) = [];
    heights(locs< 3 | locs>48) = [];
    locs(locs< 3 | locs>48) = [];
    if ~isempty(locs)
        metrics.pk1_time(j) = locs(1);
        metrics.pk1_amp(j) = pks(1);
        metrics.pk1_width(j) = width(1);
        metrics.pk1_prom(j) = prom(1);
        %metrics.pk1_height(j) = heights(1);
    end
    %if length(locs)>1
    %    metrics.pk2_time(j) = locs(2);
    %    metrics.pk2_amp(j) = pks(2);
    %    metrics.pk2_width(j) = width(2);
    %    metrics.pk2_prom(j) = prom(2);
    %    metrics.pk2_height(j) = heights(2);
    %end
end

% %look for smallest value in the hour following pk1 time
% valleys = nan(length(metrics.pk1_time),1);
% for i=1:length(metrics.pk1_time)
%     if isnan(metrics.pk1_time(i))
%         valleys(i) = nan;
%     else
%         valleys(i) = nanmin(metrics.time_series(i, metrics.pk1_time(i):metrics.pk1_time(i)+12));
%     end
% end
% 
% metrics.pk1_range = metrics.pk1_amp - valleys;

%convert pk1 time to hours
FramesPerHour = 12;
metrics.pk1_time = (metrics.pk1_time-1)/FramesPerHour;
%metrics.pk2_time = (metrics.pk2_time-1)/FramesPerHour;

function [peaks, locs, HMW, prom, heights] = global_positive_peaks(vect, baseline, num_peaks)
% seeks to find the "dominant" peaks in an input vector - output will be sorted
% from most-to-least dominant.
[all_peaks, all_locs, all_HMW, all_prom] = findpeaks(vect, 'WidthReference', 'halfheight');

all_peaks(((all_locs==1)) | (all_locs==length(vect))) = [];
all_HMW(((all_locs==1)) | (all_locs==length(vect))) = [];
all_prom(((all_locs==1)) | (all_locs==length(vect))) = [];
all_locs(((all_locs==1)) | (all_locs==length(vect))) = [];

getclosest = @(idx,vect) vect(find(abs(vect-idx)==min(abs(vect-idx)),1,'first'));
peaks = [];
locs  = [];
HMW = [];
prom = [];
heights = [];

while length(peaks) < num_peaks
    % Eliminate peaks that have been found already
    tmp_peaks = all_peaks;
    tmp_locs = all_locs;
    tmp_HMW = all_HMW;
    tmp_prom = all_prom;
    tmp_peaks(ismember(tmp_locs,locs)) = [];
    tmp_locs(ismember(tmp_locs,locs)) = [];
    tmp_HMW(ismember(tmp_locs,locs)) = [];
    tmp_prom(ismember(tmp_locs,locs)) = [];
    if isempty(tmp_peaks)
      break
    end
    
    % For each candidate, identify nearest peaks - maximize difference btw candidate and two nearest troughs.  
    diffs = zeros(size(tmp_peaks));
    loc_compare = [1 locs length(vect)];
    for l = 1:length(tmp_locs)
        tmp = loc_compare; tmp(tmp>=tmp_locs(l)) = inf;
        trough1 = min(vect(getclosest(tmp_locs(l),tmp):tmp_locs(l)));
        tmp = loc_compare; tmp(tmp<=tmp_locs(l)) = inf;
        trough2 = min(vect(tmp_locs(l):getclosest(tmp_locs(l),tmp)));
        diffs(l) = tmp_peaks(l) - max([trough1, trough2]);
    end
    
    all_HMW(tmp_peaks-baseline<0) = [];
    all_prom(tmp_peaks-baseline<0) = [];
    all_locs(tmp_peaks-baseline<0) = [];
    diffs(tmp_peaks-baseline<0) = [];
    all_peaks(tmp_peaks-baseline<0) = [];
    
    peaks = [peaks, tmp_peaks(find(diffs==max(diffs),1,'first'))];
    locs = [locs, tmp_locs(find(diffs==max(diffs),1,'first'))];
    HMW = [HMW, tmp_HMW(find(diffs==max(diffs),1,'first'))];
    prom = [prom, tmp_prom(find(diffs==max(diffs),1,'first'))];
    heights = [heights, max(diffs)];
end
peaks = peaks - baseline;
end
%% ADDITIONAL PEAK METRICS
% modified from Adelaja A

smoothed_deriv = zeros(size(derivatives));
for j = 1:size(smoothed_deriv, 1)
    smoothed_deriv(j, :) = smooth(derivatives(j, :), 'lowess');
end
FramesPerHour = 12;
pk1_frame = metrics.pk1_time * FramesPerHour + 1;
max_pk1_speed = nan(size(pk1_frame)); max_pk1_speed_frame = nan(size(pk1_frame));
for i=1:length(pk1_frame)
    if ~isnan(pk1_frame(i))
        [max_pk1_speed(i), max_pk1_speed_frame(i)] = nanmax(smoothed_deriv(i,1:pk1_frame(i)),[],2);
    end
end

metrics.max_pk1_speed = max_pk1_speed;
metrics.max_pk1_speed_frame = max_pk1_speed_frame;
%% SIGNAL STATISTICS
% modified from Adelaja A
% Calculates signal statistics for NFkB trajectories

Data = metrics.time_series(:, 1:end);
Data = fillmissing(Data,'linear', 2, 'EndValues','extrap');
smoothData = zeros(size(Data));
for i = 1:size(smoothData, 1)
    smoothData(i, :) = smooth(Data(i, :), "sgolay");
end

Fs = 12; %frames per hour
freq_range = [0.33 1];

%psd and power (scales PSD by equiv noise bandwidth of window)
[psd,fq]=pwelch(smoothData',size(Data,2),10,256,Fs,'one-sided','psd');
%metrics.fq=fq';
fq=fq';
%metrics.psd=transpose(psd./sum(psd,1));
psd=transpose(psd./sum(psd,1));

%[pwr,~]=pwelch(smoothData',size(Data,2),10,256,Fs,'one-sided','power');
%metrics.power=transpose(pwr./sum(pwr,1));
%power=transpose(pwr./sum(pwr,1));

%oscpower aka bandpower
psd = transpose(psd) ; fq = transpose(fq);% 
bp = bandpower(psd,fq,freq_range, 'psd')';
metrics.oscpower =bp;
%% INTEGRAL BASED METRICS
% modified from Adelaja A
% Calcuation of integral metrics based only off nonnegative activity
% contribtuions
    integrals_pos = zeros(size(integrals));
for i = 1:size(integrals_pos,1)
    if integrals(i,1) > 0
        integrals_pos(i,1) = integrals(i,1);
    end
    for j = 2:size(integrals_pos,2)
        if integrals(i,j)>= integrals(i,j-1)
         integrals_pos(i,j) =  integrals_pos(i,j-1) + (integrals(i,j) - integrals(i,j-1));
        elseif integrals(i,j)< integrals(i,j-1)
         integrals_pos(i,j) =  integrals_pos(i,j-1);
        end
    end
end

metrics.quarter_activity = nanmax(integrals_pos(:, 1:round(end/4)),[],2);
metrics.half_activity = nanmax(integrals_pos(:, 1:round(end/2)),[],2);
metrics.total_activity = nanmax(integrals_pos,[],2);

FramesPerHour = 12;
quarterMaxIntegral = metrics.total_activity/4;
halfMaxIntegral = metrics.total_activity/2;
threequarterMaxIntegral = metrics.total_activity*3/4;

distances = abs(integrals_pos- quarterMaxIntegral);
[~, idx] = nanmin(distances,[],2);
idx(idx==1) = NaN;
metrics.time2QuarterMaxIntegral = (idx-1)/FramesPerHour;

distances = abs(integrals_pos- halfMaxIntegral);
[~, idx] = nanmin(distances,[],2);
idx(idx==1) = NaN;
metrics.time2HalfMaxIntegral = (idx-1)/FramesPerHour;

distances = abs(integrals_pos- threequarterMaxIntegral);
[~, idx] = nanmin(distances,[],2);
idx(idx==1) = NaN;
metrics.time2ThreeQuarterMaxIntegral = (idx-1)/FramesPerHour;

%integrals_pos = integrals_pos(:, 2:end); %removes first entry which is zero
%integrals_pos = integrals_pos(:, 3:3:end);% take every third element (15 minute cummulative integrals)
%metrics.integrals_pos =integrals_pos; 
%% FOLD CHANGE METRICS
% modified from Adelaja A
fold_change = time_series_no_base_ded(:,1:end)./baseline;
metrics.max_fold_change = nanmax(fold_change,[],2);
end