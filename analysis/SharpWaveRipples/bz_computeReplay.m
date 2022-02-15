function replay = bz_computeReplay(spikes,ripples,varargin)
% Computes replay by plotting spikes time-locked to ripples

% USAGE:
%   replay = bz_computeReplay(spikes,ripples,varargin)
% 
% INPUT:
%   spikes - 
%   ripples - 
%
%
%
%
%
%
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'win',0.2,@isnumeric);
addParameter(p,'nBin',0.01,@isnumeric);

parse(p,varargin{:})

win = p.Results.win;
nBin = p.Results.nBin;
basepath = p.Results.basepath;

cd(basepath);
sessionInfo = bz_getSessionInfo(basepath,'noPrompts',true);
if ~isempty(dir([sessionInfo.FileName '.MergePoints.events.mat']))
    file = dir([sessionInfo.FileName '.MergePoints.events.mat']);
    load(file.name);
end

if ~exist('spikes','var')
    spikes = loadSpikes();
end

if ~exist('ripples','var')
    if exist(dir([sessionInfo.FileName '.ripples.Subsession.events.mat']))
        file = dir([sessionInfo.FileName '.ripples.SubSession.events.mat'])
        load(file.name)
    end
end

if ~isempty(dir([sessionInfo.FileName '.channelinfo.ripples.mat']))
    file = dir([sessionInfo.FileName '.channelinfo.ripples.mat']);
    load(file.name);
end

tracking = getSessionTracking();

%Only computing ripples for HomeCage Recordings
% Let's create first a cell with all the timestamps for each cell for each
% of the ripples
for i = 1:length(MergePoints.foldernames)
    if ~any(ismember(tracking.folders,MergePoints.foldernames{i}))
        disp(['Computing Replay in folder: ', MergePoints.foldernames{i}]);
        for j = 1:length(ripples{i}.peaks)
            ts = [ripples{i}.peaks-win/2 ripples{i}.peaks+win/2];
            for k = 1:spikes.numcells
                [status,interval,index] = InIntervals(spikes.times{k},ts);
                status_ = find(status == 1);
                interval_ = interval(status);
                for o = 1:length(status_)
                    spkdata{k}(o,:) = spikes.times{k}(status_(o)) - ripples{i}.peaks(interval_(o));
                end         
            end
            plotRasterAbad(spkdata);
        end
    end
end




end

