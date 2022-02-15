function [meanFr] = bz_meanFr(spikes,varargin)
% Computes meanFr over all the tracking folders




%% Default Parameters
p = inputParser();
% addParameter(p,'spikes',[],@bz_isCellInfo);
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'isLinearTrack',false,@islogical);

parse(p,varargin{:})
% spikes = p.Results.spikes;
basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
isLinearTrack = p.Results.isLinearTrack;

% In case meanFr already exists
if ~isempty(dir([basepath filesep '*meanFr.cellinfo.mat']))
    disp('meanFr already detected. Loading file.');
    file = dir([basepath filesep '*meanFr.cellinfo.mat']);
    load(file.name);
    return
end



%% Get Session Info
sessionInfo = bz_getSessionInfo(basepath);

%% Get Session
session = loadSession();
%% Get tracking

try
    tracking = getSessionTracking();
catch
    disp('No exist tracking... Error on computing mean Fr');
end


for jj=1:spikes.numcells
    for ii=1:length(tracking.folders)
        timestamps = tracking.events.subSessions(ii,:);
        spk = InIntervals(spikes.times{jj},timestamps);
        spk_times = spikes.times{jj}(spk);
        meanFr{jj}{ii} = length(spk_times) / (timestamps(2)-timestamps(1));
    end
end


%% Load firingMaps to include meanFr inside

if ~isempty(dir([basepath filesep sessionInfo.FileName, '*.firingMapsAvg.cellinfo.mat']))
    disp('Loading firing Maps...')
    file = dir([basepath filesep sessionInfo.FileName, '*.firingMapsAvg.cellinfo.mat']);
    load(file.name)
end

if isLinearTrack
    for i=1:length(firingMaps.stats)
        for j=1:length(firingMaps.stats{i})
            if size(firingMaps.rateMaps{i}{j},1) == 1
                % Linear Track
                firingMaps.stats{i}{j}.meanFr = meanFr{i}{1};
            else
                firingMaps.stats{i}{j}.meanFr = meanFr{i}{2};
            end
        end
    end
else
   for i=1:length(firingMaps.stats)
        for j=1:length(firingMaps.stats{i})
            firingMaps.stats{i}{j}.meanFr = meanFr{i}{j};
        end
   end 
end

if saveMat
    try
        save([basepath filesep, sessionInfo.FileName '.firingMapsAvg.cellinfo.mat'],'firingMaps')
        save([basepath filesep, sessionInfo.FileName '.meanFr.cellinfo.mat'],'meanFr')
    catch
        save([basepath filesep, sessionInfo.FileName '.firingMapsAvg.cellinfo.mat'],'firingMaps','-v7.3')
        save([basepath filesep, sessionInfo.FileName '.meanFr.cellinfo.mat'],'meanFr')
    end        
end
