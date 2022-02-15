function [meanCrossCorr] = computeMeanCrossCorrelationRegions(spikes,varargin)
% Computes mean cross correlation between all cell pairs for different
% timestamps ( epochs in a task ) or between subsessions or whatever
% between regions
% 
% USAGE
%   meanCrossCorr = computeMeanCrossCorrelation(spikes,varargin)
%
% INPUT
%   spikes - spike time cellinfo struct
%   intervals - intervals for which to compute meanCrossCorrr
%   timestamps - timestamps defining start and end of subsessions
%   binsize - binsize for ccg computation
%   duration - duration for ccg computation
%   
%
% OUTPUT
%   
%   Struct with different fields

%% Defaults and Params
p = inputParser;
addParameter(p,'intervals',[],@isnumeric);
addParameter(p,'timestamps',[],@isnumeric);
addParameter(p,'binsize',0.005,@isnumeric);
addParameter(p,'duration',0.6,@isnumeric);
addParameter(p,'foldername',[],@isstr);
addParameter(p,'showFig',true,@islogical);
addParameter(p,'saveFig',true,@islogical);

parse(p,varargin{:})

intervals = p.Results.intervals;
timestamps = p.Results.timestamps;
binsize = p.Results.binsize;
duration = p.Results.duration;
foldername = p.Results.foldername;
showFig = p.Results.showFig;
saveFig = p.Results.saveFig;

win = [-duration/2 duration/2];
sessionInfo = bz_getSessionInfo();
for i = 1:length(sessionInfo.AnatGrps)
    region{i} = sessionInfo.AnatGrps(i).region;
end
region = unique(region);

if ~isempty(timestamps)
    for i=1:length(spikes.times)
        [status,interval,index] = InIntervals(spikes.times{i},timestamps);
        spikes.times{i} = spikes.times{i}(index>0);
    end
    % Let's compute only for the same area first
    for i = 1:length(region)
        reg = region{i};
        counter = 0;
        for j = 1:length(spikes.times)
            if strcmpi(sessionInfo.AnatGrps(spikes.shankID(j)).region,reg)
                counter = counter+1;
                spks{counter} = spikes.times{j};
            end
        end
        [ccg,t] = CCG(spks,[],'binSize',binsize,'duration',duration);
        indCell = [1:size(ccg,2)];
        for k = 1:size(ccg,2)      
            zmean(:,:,k) = mean(zscore(squeeze(ccg(:,k,indCell(indCell~=k)))',[],2));
            mean_crosscorr = mean(zmean,3);
            mean_crosscorr = mean_crosscorr - min(mean_crosscorr); 
            mean_crosscorr = mean_crosscorr/max(mean_crosscorr) * (max(indCell)-1) * std(mean_crosscorr);          
        end
        if showFig
            figure;
            set(gcf,'Position',get(0,'ScreenSize'))
            plot(t,mean_crosscorr)
            xlim([win(1) win(2)])
%             if saveFig
%                 saveas(gcf,['SummaryFigures\meanCrossCorr_',foldername,'.png']); 
%             end
        end
        clear spks
    end
            
    for i = 1:length(spikes.times)
        for j = 1:length(region)
            % First region
            if strcmpi(sessionInfo.AnatGrps(spikes.shankID(i)).region,region{j})
                spks{i} = spikes.times{i};
            end
        end
    end
    
    [ccg, t] = CCG(spikes.times,[],'binSize',binsize,'duration',duration);
    indCell = [1:size(ccg,2)];

    for i=1:size(ccg,2)      
        zmean(:,:,i) = mean(zscore(squeeze(ccg(:,i,indCell(indCell~=i)))',[],2));
        mean_crosscorr = mean(zmean,3);
        mean_crosscorr = mean_crosscorr - min(mean_crosscorr); 
        mean_crosscorr = mean_crosscorr/max(mean_crosscorr) * (max(indCell)-1) * std(mean_crosscorr);          
    end
    if showFig
        figure;
        set(gcf,'Position',get(0,'ScreenSize'))
        plot(t,mean_crosscorr)
        xlim([win(1) win(2)])
        if saveFig
            saveas(gcf,['SummaryFigures\meanCrossCorr_',foldername,'.png']); 
        end
    end
    
    % OUTPUT   
    meanCrossCorr.foldername = foldername;
    meanCrossCorr.mean_crosscorr = mean_crosscorr;   
else
    [ccg, t] = CCG(spikes.times,[],'binSize',binsize,'duration',duration);
    indCell = [1:size(ccg,2)];

    for i=1:size(spikes.UID,2)      
        zmean(:,:,i) = mean(zscore(squeeze(ccg(:,i,indCell(indCell~=i)))',[],2));
        mean_crosscorr = mean(zmean,3);
        mean_crosscorr = mean_crosscorr - min(mean_crosscorr); 
        mean_crosscorr = mean_crosscorr/max(mean_crosscorr) * (max(indCell)-1) * std(mean_crosscorr);          
    end
    if showFig
        figure;
        set(gcf,'Position',get(0,'ScreenSize'))
        plot(t,mean_crosscorr)
        xlim([win(1) win(2)])
        if saveFig
            saveas(gcf,'SummaryFigures\meanCrossCorr.png')
        end
    end
    % OUTPUT   
    meanCrossCorr.mean_crosscorr = mean_crosscorr;
end

end


