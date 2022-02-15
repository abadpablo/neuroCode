% Prueba Replay Bayesian
basepath = pwd;
[sessionInfo] = bz_getSessionInfo(pwd,'noPrompts',true);
% Load spikes
spikes = loadSpikes();
% Load Ripples
% if ~isempty(dir([sessionInfo.FileName '.ripples.events.mat']))
%     file = dir([sessionInfo.FileName '.ripples.events.mat']);
%     load(file.name)
% end
if ~isempty(dir([sessionInfo.FileName '.ripples.SubSession.events.mat']))
    file = dir([sessionInfo.FileName '.ripples.SubSession.events.mat']);
    load(file.name)
end
ripples = ripples{3};

% Load firingMaps
if ~isempty(dir([sessionInfo.FileName '.firingMapsAvg.cellinfo.mat']))
    file = dir([sessionInfo.FileName '.firingMapsAvg.cellinfo.mat']);
    load(file.name)
end
template = zeros(spikes.numcells,length(firingMaps.rateMaps{1}{1}));
for i = 1:spikes.numcells
    template(i,:) = firingMaps.rateMaps{i}{1};
end

figure,
imagesc(template)
colormap(jet(15))

include = [];

[replayScores] = replay_Bayesian(spikes,ripples,template,include);
