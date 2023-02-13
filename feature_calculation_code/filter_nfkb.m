function [graph, info, measure] = filter_nfkb(id,varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Adpated from Adelaja A, Taylor B, Sheu KM, Liu Y, Luecke S, Hoffmann A. 
% Six distinct NFKB signaling codons convey discrete information to distinguish stimuli 
% and enable appropriate macrophage responses. Immunity. 2021;54(5):916-930.e7. 
% doi:10.1016/j.immuni.2021.04.011
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%% Create input parser object, add required params from function input
p = inputParser;
% Required: ID input
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||isstruct(x)||exist(x,'file'),...
    'ID input must be spreadsheet ID or full file path');
addRequired(p,'id',valid_id);

%optional parameters
valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
    'Convection correction parameter must be single integer >= 0');
addParameter(p,'ConvectionShift',1, valid_conv);
addParameter(p,'MinLifetime',97, @isnumeric);
addParameter(p,'StartThresh',2, @isnumeric);
addParameter(p, 'GraphLimits', [-0.25, 4]);
addParameter(p,'Delay',3); %% shift this number timepoint to t=0
addParameter (p,'TrimFrame',100);

% Parse parameters, assign to variables
parse(p,id, varargin{:})
MinLifetime = p.Results.MinLifetime;
start_thresh = p.Results.StartThresh;
max_shift = p.Results.ConvectionShift; % Max allowable frame shift in XY-specific correction
startFrame = p.Results.Delay;
endFrame = p.Results.TrimFrame;
graph_limits = p.Results.GraphLimits; % Min/max used in graphing
%% Load data
[measure, info] = loadID(id); 
info.ImageExpr = info.parameters.nfkbdimModule.ImageExpr;
info.graph_limits = graph_limits;
%% Filtering
robuststd = @(distr, cutoff) nanstd(distr(distr < (nanmedian(distr)+cutoff*nanstd(distr))));

% Filtering, part 1: cell fate and cytoplasmic intensity
droprows = [];
droprows = [droprows, sum(isnan(measure.NFkBdimNuclear(:,1:4)),2)>2]; % Cells existing @ expt start
droprows = [droprows, sum(isnan(measure.NFkBdimNuclear(:,1:MinLifetime)),2)>3]; % Long-lived cells
droprows = [droprows, sum(measure.NFkBdimCytoplasm(:,1:4)==0,2)>0]; % Very dim cells

% Baseline correction
nfkb = measure.NFkBdimNuclear(:,:);
nfkb = nfkb/mean(info.parameters.adj_distr(2,:));
nfkb_baseline = nanmean([nfkb(:,1:startFrame)],2);
nfkb = nfkb-nfkb_baseline;

% Filtering, part 2: eliminate outlier cells (based on mean value)
nfkb_lvl = reshape(nfkb(max(droprows,[],2) == 0,:),[1 numel(nfkb(max(droprows,[],2) == 0,:))]);
droprows =  [droprows, (nanmean(abs(nfkb-nanmean(nfkb_lvl)),2)./nanstd(nfkb_lvl))>=3];
droprows =  [droprows, (nanmean(abs(nfkb-nanmean(nfkb_lvl)),2)./nanstd(nfkb_lvl))>=1.7];

% Filtering, part 3: nuclear stain intensity and starting NFkB value
keep = max(droprows,[],2) == 0;
start_lvl = nanmin(nfkb(keep,1:3),[],2);
nuc_lvl = nanmedian(measure.MeanIntensityNuc(keep,1:31),2);
nuc_thresh = nanmedian(nuc_lvl)+2.5*robuststd(nuc_lvl(:),2);
area_thresh = 90;

droprows =  [droprows, prctile(nfkb(:,1:8),18.75,2) > start_thresh];
droprows =  [droprows, nanmedian(measure.MeanIntensityNuc(:,1:31),2) > nuc_thresh];
droprows =  [droprows, nanmedian(measure.Area,2) < area_thresh];

info.keep = max(droprows,[],2) == 0;
nfkb = nfkb(info.keep,:);
nfkb_baseline = nfkb_baseline(info.keep,:);

info.actualKeepInfo = info.keep;
info.baseline = nfkb_baseline;
%% Initialize outputs, do final corrections
graph.celldata = info.CellData(info.actualKeepInfo,:);

% Correct for XY positions that activate late
[graph.var, shift_xy] = alignTrajectories(nfkb, graph.celldata, 60, max_shift);
graph.var_nfkb_no_base_ded = graph.var + info.baseline;
[~,graph.order] = sort(nansum(graph.var(:,1:min([size(graph.var,2),endFrame])),2),'descend');

graph.shift = shift_xy;

graph.t = ((-startFrame+1)/info.parameters.FramesPerHour):(1/info.parameters.FramesPerHour):48;
graph.t = graph.t(1:endFrame);
graph.opt = maketicks(graph.t,graph_limits,0);
graph.opt.Name = 'NFkB Activation'; 

