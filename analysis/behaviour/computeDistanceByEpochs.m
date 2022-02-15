function distanceByEpochs = computeDistanceByEpochs(tracking,varargin)
% Computes the distance run by animal in epochs defined by the user (
% default 5 min)
%   
%   USAGE
%       distanceByEpochs = computeDistanceByEpochs(tracking,varargin)
%
%   INPUTS
%       tracking - buzcode tracking struct
%       
%
%
% Pablo Abad 2021.

%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'epoch',5,@isnumeric); % in minutes
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveplt',true,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
epoch = p.Results.epoch;
saveMat = p.Results.saveMat;
saveplt = p.Results.saveplt;
%% Load SessionInfo
[sessionInfo] = bz_getSessionInfo(pwd,'noPrompts',true);

% Need to find number of tracking files
numFiles = length(tracking.folders);

if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
    load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
end

[status,interval,index] = InIntervals(tracking.timestamps,MergePoints.timestamps);

a = linspace(0,tracking.timestamps(end),tracking.timestamps(end)/epoch/60);

count = 1;
for i=1:length(a)-1
    timestamp = InIntervals(tracking.timestamps,[a(i) a(i+1)]);
    if ~isempty(find(timestamp == 1))
        count = count+1;
        distance = computeDistance(tracking.position.x(timestamp),tracking.position.y(timestamp));
        distance = cumsum(distance(~isnan(distance)));
        distanceByEpochs(i) = max(distance);
    end
end

figure,
plot(a(1:end-1)/60,distanceByEpochs)
if saveplt
    saveas(gcf,['SummaryFigures\',sessionInfo.FileName,'.distanceByEpochs.png'])
end

if saveMat
    save([basepath filesep sessionInfo.FileName,'.distanceByEpochs.Behavior.mat'],'distanceByEpochs','-v7.3');
end  
end

