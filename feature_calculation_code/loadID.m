function [measure, info, AllMeasurements] = loadID(id, verbose)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Adpated from Adelaja A, Taylor B, Sheu KM, Liu Y, Luecke S, Hoffmann A. 
% Six distinct NFKB signaling codons convey discrete information to distinguish stimuli 
% and enable appropriate macrophage responses. Immunity. 2021;54(5):916-930.e7. 
% doi:10.1016/j.immuni.2021.04.011
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [measure, info, AllMeasurements] = loadID(id, options)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% LOADID pulls results from an experimental set using a "Scope Runs" Google Spreadsheet - specified in "locations.mat"
%    - choose a set by its ID number
%
% INPUTS
% id          ID# of sets get data from (or AllMeasurements.mat file location, or AllMeasurements object)
%
% OUTPUTS:
% measure          full measurement information struct
% info             general information about experiment and tracking
% AllMeasurements  originally-saved output file, augmented w/ measurement-specific information
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if nargin<2
    verbose =1 ;
end


tic
home_folder = mfilename('fullpath'); % Load locations (for images and output data)
slash_idx = strfind(home_folder,filesep);
load([home_folder(1:slash_idx(end)), 'locations.mat'],'-mat')

if ischar(id) || isnumeric(id) % Load file if a location or row index of a spreadsheet entry 
    % Find/load AllMeasurements.mat - a full file path can be specfied, or an
    % ID corresponding to an entry on the ScopeRuns spreadsheet.
    if ~exist(num2str(id), 'file') && isnumeric(id)
        data = readScopeRuns(locations.spreadsheet, id);
        info.name = [data.save_folder{1}];
        load([locations.data,filesep,data.save_dir{1},filesep,info.name,filesep,'AllMeasurements.mat'])
        info.savename = [locations.data,filesep,data.save_dir{1},filesep,info.name,filesep,'AllMeasurements.mat'];

    elseif exist(num2str(id), 'file')
        id = namecheck(id);
        load(id)
        info.savename = id;
    else
        error(['Specified file/index (''id'') is invalid'])
    end
elseif isstruct(id)
    AllMeasurements = id;
    info.savename = [locations.data,AllMeasurements.parameters.SaveDirectory,filesep,'AllMeasurements.mat'];
else
    error(['loadID accepts an "AllMeasurements" structure, or a file location/spreadsheet row index.'])
end

info.locations = locations;


% Parse AllMeasurements
info.CellData = AllMeasurements.CellData;
info.fields = fieldnames(AllMeasurements);
info.ImageDirectory = [locations.scope, AllMeasurements.parameters.ImagePath];
measure = struct;
for i = 1:length(info.fields)
    if ~strcmpi(info.fields{i},'parameters') && ~strcmpi(info.fields{i},'CellData')
        measure.(info.fields{i}) = AllMeasurements.(info.fields{i}); 
    end
end
info.fields = fieldnames(measure);

p = AllMeasurements.parameters;

info.parameters = p;

toc1 = toc;
if verbose
    disp(['Loaded "', info.savename, '" in ', num2str(round(toc1*100)/100),' sec']) 
end
