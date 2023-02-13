function [qc_metrics, metrics, aux, graph, info, measure] = nfkbQCmetrics(expt)
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
% qc_metrics    structure with only metrics needed for QC feature
%               calculations (line 517)
% metrics       structure with all other previous output metric fields
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
    metrics.baseline = info.baseline;
    % 1) basic time series. Interpolate over "normal" interval (12 frames per hr) if required
    % - Taylor B
    t = min(graph.t):1/12:max(graph.t);
    if length(t)~=length(graph.t)
        metrics.time_series = nan(size(graph.var,1),length(t));
        for i = 1:size(graph.var,1)
            metrics.time_series(i,:) = interp1(graph.t,graph.var(i,:),t);
        end
        metrics.time_series_no_base_ded = nan(size(graph.var_nfkb_no_base_ded,1),length(t));
        for i = 1:size(graph.var_nfkb_no_base_var_ded,1)
            metrics.time_series_no_base_ded(i,:) = interp1(graph.t,graph.var_nfkb_no_base_ded(i,:),t);
        end
    else
        metrics.time_series = graph.var;
        metrics.time_series_no_base_ded = graph.var_nfkb_no_base_ded;
    end
    FramesPerHour = 12;
    metrics.baseline_stdv = nanstd(metrics.time_series(:,1:Delay),0,2);
    metrics.time_series = metrics.time_series(:, Delay:end);
    metrics.time_series_no_base_ded = metrics.time_series_no_base_ded(:, Delay:end);
    metrics.t = t(Delay:end);
else
    metrics.time_series = xlsread(expt, 'time_series');
    metrics.time_series_no_base_ded = xlsread(expt, 'time_series_no_bsl_ded');
    metrics.baseline = xlsread(expt, 'baseline');
    metrics.baseline_stdv = xlsread(expt, 'bsl_stdv');
    metrics.t = xlsread(expt, 'time_h');
    FramesPerHour = 12;
end

%% BASIC METRICS
% - Taylor B
% 2) integrated activity, trapezoidal approximation
metrics.integrals = nan(size(metrics.time_series));
nan_removed = metrics.time_series;
nan_removed(isnan(nan_removed)) = 0;
for i = 1:size(metrics.integrals,1)
    metrics.integrals(i,:) = cumtrapz(metrics.t,nan_removed(i,:));
end

% 3) differentiated activity - use central finite difference
smoothed = medfilt1(metrics.time_series,3,[],2);
metrics.derivatives = (smoothed(:,3:end) - smoothed(:,1:end-2))/(1/6);

%% TRIM to a common length
endFrame = endFrame-Delay+1;
try
    metrics.time_series = metrics.time_series(:,1:endFrame);
    metrics.integrals = metrics.integrals(:,1:endFrame);
    metrics.derivatives = metrics.derivatives(:,1:(endFrame-2));
    smoothed = smoothed(:,1:endFrame);
    metrics.t = metrics.t(1:endFrame);
catch
    disp(['Note: vectors too short to cap @ ',num2str(endFrame),' frames'])
end

%% INTEGRAL WINDOWS
% Integrals within one-hour windows (0-1, 1-2, 2-3) and three hour windows
% (0-3, 1-4, etc) of activity - Taylor B
max_hr = floor(max(metrics.t));
metrics.intwin1 = nan(size(metrics.time_series,1),max_hr);
metrics.intwin3 = nan(size(metrics.time_series,1),max_hr-2);
for i = 1:(max_hr)
    win = metrics.t>=(i-1) & metrics.t<(i);
    metrics.intwin1(:,i) = trapz(metrics.t(win),metrics.time_series(:,win),2);
    if i<= (max_hr-2)
        win = metrics.t>=(i-1) & metrics.t<(i+2);
        metrics.intwin3(:,i) = trapz(metrics.t(win),metrics.time_series(:,win),2);
    end
end

%% MAX/MIN/MED METRICS
% - Taylor B
metrics.max_amplitude = nanmax(metrics.time_series(:,1:end),[],2);
metrics.max_integral = nanmax(metrics.integrals(:,1:end),[],2);
metrics.max_derivative = nanmax(metrics.derivatives(:,1:end),[],2);
metrics.min_derivative = nanmin(metrics.derivatives(:,1:end),[],2);
metrics.median_derivative = nanmedian(metrics.derivatives(:,1:end),2);
metrics.mean_derivative = nanmean(metrics.derivatives(:,1:end),2);           

%% ACTIVITY METRICS
% Determine which cells are responders and compute off-times - Luecke S
Wliml = 1; %first/lower time point of window to check for activity
Wlimu = 48; %last/upper time point of window to check for activity, ie check in the first 4 hours after stimulation
blockLengthThresh = 5; %number of consecutive frames cell needs to pass activity threshold to be considered a responder
OnThresh = 3; %value set as the standard deviation theshold above baseline reqd to be considered a "responder" (also called sigmaThresh)
smoothed_by_sigma = smoothed./metrics.baseline_stdv;
[metrics.responder_index, metrics.responders_fraction, metrics.off_times] = get_activity_metrics(smoothed_by_sigma, Wliml, Wlimu, OnThresh, blockLengthThresh);
metrics.off_times = metrics.off_times/FramesPerHour;
metrics.off_times(metrics.off_times<0) = 0;

%% METRICS OF DURATION
% Envelope width: maximum consecutive time above a threshold (envelope must
% begin within 1st 6 hrs) - Taylor B
smoothed2 = medfilt1(metrics.time_series,5,[],2);
lowerThresh = 0;
upperThresh = 7.1;
aux.thresholds = linspace(lowerThresh,upperThresh,25);
metrics.envelope = zeros(size(metrics.time_series,1),25);
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

% Find total oscillatory content of particular cells (using thresholds from 0.35 to 0.7 hrs^(-1))
freq_thresh = aux.freq( (aux.freq >= (0.35/3600)) & (aux.freq <= (0.7/3600)));
metrics.oscfrac = nan(size(aux.power,1),length(freq_thresh));
for j = 1:length(freq_thresh)
    for i =1:size(metrics.time_series,1)
        metrics.oscfrac(i,j) = nansum(aux.power(i,aux.freq >= freq_thresh(j))) /nansum(aux.power(i,:));
        if isnan(metrics.oscfrac(i,j))
            metrics.oscfrac(i,j) = 0;
        end
    end
end

%% METRICS OF AMPLITUDE AND TIMING
% modified from Taylor B to extract only positive value peaks
% 1st + 2nd peak time/amplitude/prominence/width/height

pk_feats = {'pk1_amp', 'pk1_time', 'pk1_width', 'pk1_prom', 'pk1_height',...
        'pk2_amp', 'pk2_time', 'pk2_width', 'pk2_prom', 'pk2_height'};
for j=1:length(pk_feats)
    metrics.(pk_feats{j}) = nan(size(metrics.time_series_no_base_ded,1),1);
end
for j = 1:size(metrics.pk1_time,1)    
    [pks, locs, width, prom, heights] = global_positive_peaks(metrics.time_series_no_base_ded(j,1:min([90,MinLifetime])),metrics.baseline(j),5);
    % Supress any peaks that are within 6 frames of each other.
    [locs, order] = sort(locs,'ascend');
    pks = pks(order); width = width(order); prom = prom(order); heights = heights(order);
    while min(diff(locs))<6
        tmplst = find(diff(locs)==min(diff(locs)),1,'first');
        tmplst = tmplst + (pks(tmplst)>=pks(tmplst+1));
        pks(tmplst) = []; locs(tmplst) = []; width(tmplst) = []; prom(tmplst) = []; heights(tmplst) = [];
    end
    pks(locs< 3) = [];
    locs(locs< 3) = [];
    width(locs< 3) = [];
    prom(locs< 3) = [];
    heights(locs< 3) = [];
    if ~isempty(locs)
        metrics.pk1_time(j) = locs(1);
        metrics.pk1_amp(j) = pks(1);
        metrics.pk1_width(j) = width(1);
        metrics.pk1_prom(j) = prom(1);
        metrics.pk1_height(j) = heights(1);
    end
    if length(locs)>1
        metrics.pk2_time(j) = locs(2);
        metrics.pk2_amp(j) = pks(2);
        metrics.pk2_width(j) = width(2);
        metrics.pk2_prom(j) = prom(2);
        metrics.pk2_height(j) = heights(2);
    end
end
FramesPerHour = 12;
metrics.pk1_time = (metrics.pk1_time-1)/FramesPerHour;
metrics.pk2_time = (metrics.pk2_time-1)/FramesPerHour;

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
%% MISC PEAK METRICS
% - Adelaja A

metrics.pk2_ratio           = (metrics.pk2_amp)./(metrics.pk1_amp);
metrics.pk2_ratio_prom      = (metrics.pk2_prom)./(metrics.pk1_prom);

smoothed = zeros(size(metrics.derivatives));
for j = 1:size(smoothed, 1)
    smoothed(j, :) = smooth(metrics.derivatives(j, :), 'lowess');
end
FramesPerHour = 12;
pk1_frame = metrics.pk1_time * FramesPerHour + 1;
max_pk1_speed = nan(size(pk1_frame)); max_pk1_speed_frame = nan(size(pk1_frame));
for i=1:length(pk1_frame)
    if ~isnan(pk1_frame(i))
        [max_pk1_speed(i), max_pk1_speed_frame(i)] = nanmax(smoothed(i,1:pk1_frame(i)),[],2);
    end
end

metrics.max_pk1_speed = max_pk1_speed;

%% PEAK STATISTICS
% Calculates peak statistics - modified from Adelaja A

min_pks = 2;
max_pk_diff= 35;
FramesPerHour = 12;

[metrics.peak_times,metrics.peak_amps, metrics.valley_times, metrics.valley_amps]=nfkbpeaks(metrics.time_series_no_base_ded(:, 1:end), metrics.baseline, 'BeginFrame',3,'MinHeight',0.75,'MinDist',6);
ipt=diff(metrics.peak_times,1,2);

ipt(ipt>max_pk_diff)=nan;
tot_pks = sum(~isnan(ipt),2)+1;
metrics.kept = true([size(ipt,1),1]);
ipt(tot_pks<min_pks,:) = [];
metrics.kept(tot_pks<min_pks)=false;
ipt = ipt.*60/FramesPerHour;%convert to minutes

metrics.mean_ipt    = nanmean(ipt,2);
metrics.median_ipt  = nanmedian(ipt,2);
metrics.std_ipt     = nanstd(ipt,[],2);
metrics.max_ipt     = nanmax(ipt,[],2);
metrics.min_ipt     = nanmin(ipt,[],2);
metrics.cv_ipt      = metrics.std_ipt./metrics.mean_ipt;
metrics.ipt         = ipt;
metrics.num_peaks   = tot_pks;

metrics.mean_peak_amp   = nanmean(metrics.peak_amps, 2); 
metrics.median_peak_amp = nanmedian(metrics.peak_amps, 2); 
metrics.std_peak_amp    = nanstd(metrics.peak_amps, [],2); 
metrics.cv_peak_amp     = metrics.std_peak_amp./metrics.mean_peak_amp; 

metrics.peak2trough = metrics.peak_amps(:,1:end-1)-metrics.valley_amps;
minVals = zeros(size(metrics.time_series_no_base_ded,1),1);
for row=1:size(minVals,1)
    pk_frame =metrics.peak_times(row, 1);
    if isnan(pk_frame)
        continue;
    end
    minVals(row) = min(metrics.time_series(row, 1:pk_frame));
end
troughs = [minVals, metrics.valley_amps];

metrics.trough2peak         = troughs-metrics.peak_amps;

metrics.max_peak2trough     = nanmax(metrics.peak2trough, [],2);
metrics.max_trough2peak     = nanmax(metrics.trough2peak, [],2);

metrics.mean_peak2trough    = nanmean(metrics.peak2trough, 2);
metrics.mean_trough2peak    = nanmean(metrics.trough2peak,2);

metrics.median_peak2trough  = nanmedian(metrics.peak2trough,2); 
metrics.median_trough2peak  = nanmedian(metrics.trough2peak,2);

metrics.min_peak2trough     = nanmin(metrics.peak2trough,[],2); 
metrics.min_trough2peak     = nanmin(metrics.trough2peak, [],2);

metrics.std_trough2peak     = nanstd(metrics.trough2peak,[],2);
metrics.std_peak2trough     = nanstd(metrics.peak2trough,[],2); 
metrics.cv_trough2peak      = metrics.std_trough2peak./metrics.mean_trough2peak;
metrics.cv_peak2trough      = metrics.std_peak2trough./metrics.mean_peak2trough; 

%% SIGNAL STATISTICS
% Calculates signal statistics for NFkB trajectories - Adelaja A

Data = metrics.time_series(:, 1:end);
Data = fillmissing(Data,'linear', 2, 'EndValues','extrap');
smoothData = zeros(size(Data));
for i = 1:size(smoothData, 1)
    smoothData(i, :) = smooth(Data(i, :), "sgolay");
end

Fs = 12; %frames per hour
freq_range = [0.33 1];

metrics.powerbw         = powerbw(smoothData',Fs, freq_range)';
metrics.medfreq         = medfreq(smoothData', Fs,freq_range)';
metrics.meanfreq        = meanfreq(smoothData', Fs, freq_range)'; 
metrics.peak2rms        = peak2rms(smoothData,2);     
metrics.rms             = rms(smoothData')';
metrics.peak2peak       = peak2peak(smoothData,2);
metrics.mean_movmad     = mean( movmad(smoothData,3,2),2);
metrics.mean_movstd     = mean( movstd(smoothData,3,[],2),2);
metrics.mean_movvar     = mean( movvar(smoothData,3,[],2),2);

%psd and power (scales PSD by equiv noise bandwidth of window)
[psd,fq]=pwelch(smoothData',size(Data,2),10,256,Fs,'one-sided','psd');
metrics.fq=fq';
metrics.psd=transpose(psd./sum(psd,1));

[pwr,~]=pwelch(smoothData',size(Data,2),10,256,Fs,'one-sided','power');
metrics.power=transpose(pwr./sum(pwr,1));

%oscpower aka bandpower
psd = transpose(metrics.psd) ; fq = transpose(metrics.fq);% 
bp = bandpower(psd,fq,freq_range, 'psd')';
metrics.oscpower =bp;

%oscillation frequency--find peaks within the frequency range
bandfilter= @(x) x<= max(freq_range) & x>= min(freq_range);normalize =@(x) x/sum(x);
ix =bandfilter(metrics.fq);
peakFun =@(a) arrayfun(@(j) findpeaks(a(j,:), metrics.fq(ix),...
                'SortStr', 'descend', 'MinPeakProminence', 0.0055), 1:size(a,1), 'UniformOutput',false);
[peaks,locs] = peakFun(metrics.psd(:,ix)) ; %peaks = psd, locs = frequency
freq =zeros(size(peaks)); 
    for j = 1:numel(peaks)
        if numel(peaks{j}) > 1
        % more than one peak within the range, take weighted
        % sum of frequency 
            wgts = normalize(peaks{j}); 
            freq(j) = sum(locs{j}.*wgts); 
        elseif ~isempty(peaks{j})
            freq(j) = locs{j}; 
        end
    end
metrics.oscfreq = freq';

%oscillatory bandwidth       
metrics.oscbandwidth     = obw(smoothData',Fs)';

%max entropy
max_entropy= zeros(size(smoothData,1),1);    
time_pts = max_entropy;
for j =1:numel(max_entropy)
  [sig_entropy, tp]= pentropy(smoothData(j,:), Fs/3600,'FrequencyLimit',freq_range./3600);   
   [max_entropy(j), ix] = max(sig_entropy);
   time_pts(j) = tp(ix)/3600;
end
metrics.max_entropy = max_entropy; 
metrics.noise_est =  wnoisest(Data)';

%% INTEGRAL BASED METRICS
% modified from Adelaja A
% Calcuation of integral metrics based only off nonnegative activity
% contribtuions
    integrals_pos = zeros(size(metrics.integrals));
for i = 1:size(integrals_pos,1)
    if metrics.integrals(i,1) > 0
        integrals_pos(i,1) = metrics.integrals(i,1);
    end
    for j = 2:size(integrals_pos,2)
        if metrics.integrals(i,j)>= metrics.integrals(i,j-1)
         integrals_pos(i,j) =  integrals_pos(i,j-1) + (metrics.integrals(i,j) - metrics.integrals(i,j-1));
        elseif metrics.integrals(i,j)< metrics.integrals(i,j-1)
         integrals_pos(i,j) =  integrals_pos(i,j-1);
        end
    end
end

metrics.integrals_pos =integrals_pos; 
metrics.max_pos_integral = nanmax(integrals_pos,[],2);

endFrame = min(96, size(metrics.integrals,2));
halfMaxIntegral = nanmax(metrics.integrals_pos(:,1:endFrame),[],2)/2;
distances = abs(metrics.integrals_pos(:,1:endFrame)- halfMaxIntegral);
[~, idx] = nanmin(distances,[],2);
idx(idx==1) = NaN;
FramesPerHour = 12;
metrics.time2HalfMaxIntegral = (idx-1)/FramesPerHour;

%% FOLD CHANGE METRICS
% modified from Adelaja A
metrics.fold_change = metrics.time_series_no_base_ded(:,1:end)./metrics.baseline;
metrics.max_fold_change = nanmax(metrics.fold_change,[],2);
metrics.max_value = nanmax(metrics.time_series_no_base_ded - metrics.baseline, [], 2);

%% OSCILATORY CATEGORY
% Calculates the percentage of trajectories that fall into different
% oscillation categories- Adelaja A

cutoff_fq = 0.42; % freq cutoff for oscillatory & non-oscillatory
cats = zeros(size(metrics.peakfreq)); 
cats(metrics.off_times==0) = 1;
cats((metrics.off_times>0)&(metrics.peakfreq<cutoff_fq)) = 3;
cats((metrics.peakfreq>=cutoff_fq))=2;
varNames = {'off', 'osc', 'non_osc'};
metrics.osc_cats = categorical(cats,1:3, varNames); 

%% COLLECT METRICS NEEDED FOR QC FEATURES CALCULATIONS
% Duration
qc_metrics.dur_t = metrics.duration(:, 2);
qc_metrics.n_pks = metrics.num_peaks;
%EarlyVsLate
qc_metrics.eVl = metrics.time2HalfMaxIntegral;
%OscVsNon
qc_metrics.oVn = metrics.oscpower;
%PeakAmplitude
qc_metrics.max_val = metrics.max_value;
qc_metrics.pk2pk = metrics.peak2peak;
qc_metrics.pk1_amp = metrics.pk1_amp;
%Speed 
qc_metrics.max_pk1_spd = metrics.max_pk1_speed;
qc_metrics.pk1_t = metrics.pk1_time;
qc_metrics.deriv2 = metrics.derivatives(:, 2);
%TotalActivity
qc_metrics.tot_act = metrics.max_pos_integral;
end