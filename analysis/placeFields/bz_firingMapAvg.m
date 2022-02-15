function [firingMaps] = bz_firingMapAvg(positions,spikes,varargin)

% USAGE
% [firingMaps] = bz_firingMapAvg(positions,spikes,varargin)
% Calculates averaged firing map for a set of linear postions 
%
% INPUTS
%
%   spikes    - buzcode format .cellinfo. struct with the following fields
%               .times 
%   positions - [t x y ] or [t x] position matrix or
%               cell with several of these matrices (for different conditions)
%      or
%   behavior  - buzcode format behavior struct - 
%   <options>      optional list of property-value pairs (see table below)
% ===================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'			smoothing size in bins (0 = no smoothing, default = 2)
%     'nBins'			number of bins (default = 50)
%     'speedThresh'		speed threshold to compute firing rate
%     'minTime'			minimum time spent in each bin (in s, default = 0)
%     'mode'			interpolate' to interpolate missing points (< minTime),
%                   	or 'discard' to discard them (default)
%     'maxDistance'		maximal distance for interpolation (default = 5)
%     'maxGap'			z values recorded during time gaps between successive (x,y)
%                   	samples exceeding this threshold (e.g. undetects) will not
%                	    be interpolated; also, such long gaps in (x,y) sampling
%                 	    will be clipped to 'maxGap' to compute the occupancy map
%                 	    (default = 0.100 s)
%     'orderKalmanVel'	order of Kalman Velocity Filter (default 2)
%     'saveMat'   		- logical (default: false) that saves firingMaps file
%     'CellInspector'  	- logical (default: false) that creates an otuput
%                   	compatible with CellInspector

%
%
% OUTPUT
%
%   firingMaps - cellinfo struct with the following fields
%                .rateMaps              gaussian filtered rates
%                .rateMaps_unsmooth     raw rate data
%                .rateMaps_box          box filtered rates
%                .countMaps             raw spike count data
%                .occuMaps              position occupancy data
%                .cmBin                 cm/bins ratio
%
% Antonio FR, 10/2019

%% parse inputs
p=inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'smooth_buz',2,@isnumeric);
addParameter(p,'smooth_tint',7,@isnumeric);
addParameter(p,'speedThresh',0.05,@isnumeric);
addParameter(p,'nBins',50,@isnumeric);
addParameter(p,'maxGap',0.1,@isnumeric);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'CellInspector',false,@islogical);
addParameter(p,'mode','discard',@isstr);
addParameter(p,'maxDistance',5,@isnumeric);
addParameter(p,'orderKalmanVel',2,@isnumeric);
addParameter(p,'speedFilter',false,@islogical);
addParameter(p,'cmBin',2.5,@isnumeric);
addParameter(p,'periodicAnalysis',false,@islogical);
addParameter(p,'numRand',1000,@isnumeric);
addParameter(p,'spikeShuffling',true,@islogical);
addParameter(p,'buzAnalysis',false,@islogical);
addParameter(p,'tintAnalysis',true,@islogical);

parse(p,varargin{:});
smooth_buz = p.Results.smooth_buz;
smooth_tint = p.Results.smooth_tint;
speedThresh = p.Results.speedThresh;
nBins = p.Results.nBins;
maxGap = p.Results.maxGap;
minTime = p.Results.minTime;
saveMat = p.Results.saveMat;
CellInspector = p.Results.CellInspector;
mode = p.Results.mode;
maxDistance = p.Results.maxDistance;
order = p.Results.orderKalmanVel;
speedFilter = p.Results.speedFilter;
cmBin = p.Results.cmBin;
periodicAnalysis = p.Results.periodicAnalysis;
numRand = p.Results.numRand;
basepath = p.Results.basepath;
spikeShuffling = p.Results.spikeShuffling;
buzAnalysis = p.Results.buzAnalysis;
tintAnalysis = p.Results.tintAnalysis;

%% In case firingMapsAvg already exists
if ~isempty(dir([basepath filesep '*firingMapsAvg.cellinfo.mat']))
    disp('Firing Maps already detected! Loading file.');
    file = dir([basepath filesep '*firingMapsAvg.cellinfo.mat']);
    load(file.name);
    return
end
% We need to obtain the number of bins
tracking = getSessionTracking();
if isfield(tracking,'apparatus')
    if size(tracking.apparatus,2) == 1
        numApparatus = 1;
        apparatus = tracking.apparatus{1};
        apparatus_name = apparatus.name;
        xLim = apparatus.boundingbox.xmax - apparatus.boundingbox.xmin;
        yLim = apparatus.boundingbox.ymax - apparatus.boundingbox.ymin;  
    else
        numApparatus = size(tracking.apparatus,2);
        for i=1:numApparatus
            apparatus{i} = tracking.apparatus{i};
            apparatus_name{i} = apparatus{i}.name;
            xLim{i} = apparatus{i}.boundingbox.xmax - apparatus{i}.boundingbox.xmin;
            yLim{i} = apparatus{i}.boundingbox.ymax - apparatus{i}.boundingbox.ymin;
        end
    end
end

if numApparatus == 1
    clear nBins
    nBins{1} = round(xLim/cmBin);
else
    clear nBins
    for i = 1:numApparatus
        nBins{i} = round(xLim{i}/cmBin);
    end
end

if isstruct(positions)
    positions = positions.maps;
end

% number of conditions
  if iscell(positions)
     conditions = length(positions); 
  elseif isvector(positions)
     conditions = 1;
  end
        
  %%% TODO: conditions label
%   if conditions == 3
%       nBins(1) = numberBins{1};
%       nBins(2) = numberBins{1};
%       nBins(3) = numberBins{2};
%   end
  

% for i=1:length(positions)   
%     positions{i} = positions{i}(tracking.events.subSessionsMask == i,:);  
% end
%% Calculate
% Erase positions below speed threshold
for iCond = 1:size(positions,2)
    % Compute speed
    post = positions{iCond}(:,1);
    % - 1D 
    if size(positions{iCond},2)==2
        posx = positions{iCond}(:,2);
        [~,~,~,vx,vy,~,~] = KalmanVel(posx,posx*0,post,order);
    elseif size(positions{iCond},2)==3
        posx = positions{iCond}(:,2);
        posy = positions{iCond}(:,3);
        [~,~,~,vx,vy,~,~] = KalmanVel(posx,posy,post,order);
    else
        warning('This is not a linear nor a 2D space!');
    end
    % Absolute speed
    v = sqrt(vx.^2+vy.^2);    
    if length(v) < length(posx)
        dif = length(posx) - length(v);
        s = v;
        sx = vx;
        sy = vy;
        clear v
        clear vx
        clear vy
        
        v = zeros(1,length(posx));
        vx = zeros(1,length(posx));
        vy = zeros(1,length(posx));
        
        v = [nan(dif,1); s];
        vx = [nan(dif,1); sx];
        vy = [nan(dif,1); sy];        
    end
%     tracking.speed.v = v/100; % cm/s
%     tracking.speed.vx = vx/100; % cm/s
%     tracking.speed.vy = vy/100; % cm/s
%     [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); sessionInfo.rates.lfp = 1250;  
%     basepath = sessionInfo.session.path;
%     save([basepath filesep sessionInfo.FileName '.Tracking.Behavior.mat'],'tracking');
    % Compute timestamps where speed is under threshold
    if speedFilter == 1
        positions{iCond}(v<speedThresh,:) = [];
    end
end

% Need to take into account that probably we have 3 conditions ( linear
% Track has 2) and we need to change nBins and apparatus
islinearTrack = false;
for i = 1:length(tracking.apparatus)
    if strcmp(tracking.apparatus{i}.name,'Linear Track  N-S')
        islinearTrack = true;
    end
end
    
if islinearTrack && conditions ~= length(nBins)
    bins = nBins;
    clear nBins
    nBins{1} = bins{1};
    nBins{2} = bins{1};
    nBins{3} = bins{2};  
    
    tracking.apparatus{3} = tracking.apparatus{2};
    tracking.apparatus{2} = tracking.apparatus{1};
    
    tracking.pixelsmetre{3} = tracking.pixelsmetre{2};
    tracking.pixelsmetre{2} = tracking.pixelsmetre{1};
    
    tracking.zone{3} = tracking.zone{2};
    tracking.zone{2} = tracking.zone{1};
    
    tracking.folders{3} = tracking.folders{2};
    tracking.folders{2} = tracking.folders{1};
end


%% get firign rate maps & map stats
for unit = 1:length(spikes.times)
    for c = 1:conditions
        if size(positions{c},2) == 2 
            % Only timestamps and linear position ( LinearTrack)
            map{unit}{c} = Map(positions{c},spikes.times{unit},'smooth',smooth_buz,'minTime',minTime,...
                'nBins',nBins{1},'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
            stats{unit}{c} = MapStats(map{unit}{c},spikes.times{unit},'nBins',nBins{c},'verbose','off');

        elseif size(positions{c},2) == 3
            % Timestamps + x and y positions (Open Field)
            if buzAnalysis
                map{unit}{c} = Map(positions{c},spikes.times{unit},'smooth',smooth_buz,'minTime',minTime,...
                    'nBins',nBins{c},'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
            elseif tintAnalysis
%                 map_tint{unit}{c} = firingMapsBuild(positions{c},spikes.times{unit},'smooth',smooth,'minTime',minTime,...
%                     'nBins',nBins{c},'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
                [n,bin] = histc(spikes.times{unit},positions{c}(:,1));
                bin = bin(bin >0);
                bin = unique(bin);
                map{unit}{c} = firingMapsBuild_tint(positions{c},bin,'boundingbox',tracking.apparatus{c}.boundingbox,'pixels_metre',tracking.pixelsmetre{c},...
                    'tracking',tracking,'smooth',smooth_tint,'minTime',minTime,...
                    'nBins',nBins{c},'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);                
            end
            
            stats{unit}{c} = MapStats(map{unit}{c},spikes.times{unit},'nBins',nBins{c},'verbose','off');

            % shuffling of spikes times
            if spikeShuffling
                if buzAnalysis
                    disp(['Shuffling unit:', num2str(unit)])
                    shuffling{unit}{c} = bz_SpikeShuffling(positions{c},spikes.times{unit},'smooth',smooth_buz,'minTime',minTime,...
                        'nBins',nBins{c},'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance,'numRand',numRand);
                elseif tintAnalysis
                    disp(['Shuffling unit:', num2str(unit), ' out of ', num2str(length(spikes.times))])
                    [shuffling{unit}{c},shuffling_maps{unit}{c}] = bz_SpikeShuffling(positions{c},spikes.times{unit},'smooth',smooth_tint,'minTime',minTime,...
                        'nBins',nBins{c},'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance,'numRand',numRand);
                end

            else
                shuffling = []; shuffling_maps = [];
            end
            % Periodic Firing
            if periodicAnalysis
                periodic{unit}{c} = bz_PeriodicPower(map{unit}{c},shuffling_maps{unit}{c},'random',true,'nBins',nBins{c});
            else
                periodic{unit}{c} = [];
            end 

        end      
    end
end
% cmBin = (max(positions{1}(:,2))-min(positions{1}(:,2)))/nBins;
%%% TODO: pass rest of inputs to Map

%% restructure into cell info data type

% inherit required fields from spikes cellinfo struct
firingMaps.UID = spikes.UID;
try firingMaps.sessionName = spikes.sessionName;
catch
    firingMaps.sessionName = spikes.basename;
end
try
firingMaps.region = spikes.region; 
catch
   %warning('spikes.region is missing') 
end

if buzAnalysis
    firingMaps.params.analysis = 'buzcode';
    firingMaps.params.smooth = smooth_buz;
elseif tintAnalysis
    firingMaps.params.analysis = 'tint';
    firingMaps.params.smooth = smooth_tint;
end

% firingMaps.params.smooth_buz = smooth_buz;
firingMaps.params.smooth_tint = smooth_tint;
firingMaps.params.minTime = minTime;
firingMaps.params.nBins = nBins;
firingMaps.params.maxGap = maxGap;
firingMaps.params.mode = mode;
firingMaps.params.maxDistance = maxDistance;
firingMaps.params.cmBin = cmBin;
firingMaps.params.numRand = numRand;

for unit = 1:length(spikes.times)
    for c = 1:conditions
    firingMaps.rateMaps{unit,1}{c} = map{unit}{c}.z;
    firingMaps.countMaps{unit,1}{c} = map{unit}{c}.count;
    firingMaps.occupancy{unit,1}{c} = map{unit}{c}.time;
    end
end

% Save stats
firingMaps.stats = stats;

% Save shuffling
firingMaps.shuffling = shuffling;
firingMaps_shuffling = shuffling_maps;


% Save periodic
if periodicAnalysis
    firingMaps.periodic = periodic;
end


for unit = 1:length(spikes.times)
    for c = 1:conditions
        firingMaps.rateMapsUnSmooth{unit,1}{c} =map{unit}{c}.zUnSmooth;
        firingMaps.countMapsUnSmooth{unit,1}{c} = map{unit}{c}.countUnSmooth;
        firingMaps.occupancyUnSmooth{unit,1}{c} = map{unit}{c}.timeUnSmooth;
    end
end


if saveMat
    try
        save([firingMaps.sessionName '.firingMapsAvg.cellinfo.mat'],'firingMaps'); 
    catch
        save([firingMaps.sessionName '.firingMapsAvg.cellinfo.mat'],'firingMaps','-v7.3'); 
    end
    
    try
        save([firingMaps.sessionName '.firingMapsAvg_shuffling.cellinfo.mat'],'firingMaps_shuffling');
    catch
        save([firingMaps.sessionName '.firingMapsAvg_shuffling.cellinfo.mat'],'firingMaps_shuffling','-v7.3');
    end
end


end
