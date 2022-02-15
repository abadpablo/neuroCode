function bz_PlotRippleStats_abad(maps,data,stats,varargin)

%bz_PlotRippleStats - Plot descriptive stats for ripples (100~200Hz oscillations).
%
%  USAGE
%
%    bz_PlotRippleStats(maps,data,stats,<options>)
%
%    Use the ouputs of <a href="matlab:help RippleStats">RippleStats</a> as input parameters. Using cell arrays of
%    repeated measures will yield pairwise statistics.
%
%    maps           ripple instantaneous amplitude, frequency, and phase maps
%    data           frequency and amplitude at peak, durations
%    stats          autocorrelogram and correlations
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'frequency'   sampling rate (in Hz) (default = 1250Hz)
%     'durations'   durations before and after ripple peak (in s)
%                   (default = [-0.5 0.5])
%    =========================================================================
%
%  SEE
%
%    See also FindRipples, RippleStats, SaveRippleEvents.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Script adapted by Pablo Abad in order to change and add features to the
% figures

% Default values
p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'samplingRate',1250,@isnumeric);
addParameter(p,'durations',[-50 50]/1000,@isnumeric);
addParameter(p,'nCorrBins',1000,@isnmueric);
addParameter(p,'showFig',true,@islogical);
addParameter(p,'saveFig',true,@islogical);
addParameter(p,'timestamps',[],@isnumeric);
addParameter(p,'foldername',[],@isstr);
addParameter(p,'ripple',[],@isstruct);

parse(p,varargin{:});

basepath = p.Results.basepath;
samplingRate = p.Results.samplingRate;
durations = p.Results.durations;
nCorrBins = p.Results.nCorrBins;
showFig = p.Results.showFig;
saveFig = p.Results.saveFig;
timestamps = p.Results.timestamps;
foldername = p.Results.foldername;
ripple = p.Results.ripple;

sessionInfo = bz_getSessionInfo(basepath);
session = loadSession();
if ~isempty(dir([basepath filesep sessionInfo.FileName '*.MergePoints.events.mat']))
    disp('Loading MergePoints...')
    file = dir([basepath filesep sessionInfo.FileName '*.MergePoints.events.mat']);
    load(file.name);
end



if isa(maps,'cell')
else % remove the completetely, eventually...
	m{1} = maps;
	d{1} = data;
	s{1} = stats;
    maps = m;
    data = d;
    stats = s;
end
nElectrodes = length(maps);

%  nBins = floor(samplingRate*diff(durations)/2)*2+1; % must be odd
nBins = 	size(maps{1}.ripples,2);
nHalfCenterBins = 3;
%  centerBin = durations(1)/sum(durations)*(nBins-1)+1;
centerBin = ceil(nBins/2);
centerBins = centerBin-nHalfCenterBins:centerBin+nHalfCenterBins;

% Let's break the ripples to see in which paradigm they occur
if ~isempty(ripple)
    ts = ripple.peaks;
    [status,a,b] = InIntervals(ts,MergePoints.timestamps);
    for i=1:max(a)
        c = find(a == i);
        if ~isempty(c)
            endparadigm(i) = c(end);
        else
            endparadigm(i) = NaN;
        end
    end
    dx = diff(durations);
    x = durations(1):dx/size(maps{1}.ripples,2):durations(2);
    endparadigm_vector = ones(max(a),length(x));
    for i=1:length(endparadigm)
        endparadigm_vector(i,:) = endparadigm_vector(i,:)*endparadigm(i);
    end
end
% Stats for individual electrodes
for electrode = 1:nElectrodes

	% Unsorted color plots
	f = figure;
	set(f,'name',['Ripples (unsorted) - ' int2str(electrode)]);
	dx = diff(durations);
	x = durations(1):dx/size(maps{electrode}.ripples,2):durations(2);
	subplot(2,2,1);PlotColorMap(maps{electrode}.ripples,1,'bar','on',...
        'cutoffs',[min(maps{electrode}.ripples(:)) max(maps{electrode}.ripples(:))],...
        'x',x);xlabel('Ripples');
    if exist('endparadigm_vector','var')
        hold on;
        plot(x,endparadigm_vector,'k')
    end
	subplot(2,2,2);PlotColorMap(maps{electrode}.phase,1,'bar','on','type',...
        'circular','x',x);xlabel('Ripple Phase');
    if exist('endparadigm_vector','var')
        hold on;
        plot(x,endparadigm_vector,'k')
    end
	subplot(2,2,3);PlotColorMap(maps{electrode}.frequency,1,'bar','on',...
        'cutoffs',[100 250],'x',x);xlabel('Ripple Frequency');
    if exist('endparadigm_vector','var')
        hold on;
        plot(x,endparadigm_vector,'k')
    end
	subplot(2,2,4);PlotColorMap(maps{electrode}.amplitude,1,'bar','on',...
        'x',x);xlabel('Ripple Amplitude'); 
    if exist('endparadigm_vector','var')
        hold on;
        plot(x,endparadigm_vector,'k')
    end
    
    if ~isempty(foldername)
        saveas(gcf,['SummaryFigures\Ripples Stats1_',foldername,'.png']);
    else
        saveas(gcf,'SummaryFigures\Ripples Stats1.png');
    end

	% Ripples stats: ripples, peak frequency vs amplitude, autocorrelogram, peak frequency vs duration, peak amplitude vs duration
	f = figure;set(f,'name',['Ripple Stats - ' int2str(electrode)]);

	subplot(2,2,1);a = gca;hold on;
	plot(((1:nBins)'-ceil(nBins/2))/nBins*diff(durations),maps{electrode}.ripples','b');

	subplot(2,2,2);
	b = bar(stats{electrode}.acg.t,stats{electrode}.acg.data);set(b,'FaceColor',[0 0 0]);xlabel('Autocorrelogram');
%  	b = bar(((0:nCorrBins)-nCorrBins/2)/1000,stats{electrode}.acg.data);xlim([-nCorrBins nCorrBins]/2000);set(b,'FaceColor',[0 0 0]);xlabel('Autocorrelogram');

	subplot(2,3,4);a = gca;
	PlotDistribution2(data{electrode}.peakAmplitude,data{electrode}.peakFrequency,'nbins',1000); %,'smooth',5
	axes(a);xlabel(['r=' num2str(stats{electrode}.amplitudeFrequency.rho(1,2)) ' p=' num2str(stats{electrode}.amplitudeFrequency.p(1,2))]);ylabel('Frequency vs Amplitude');

	subplot(2,3,5);a = gca;
	PlotDistribution2(data{electrode}.duration,data{electrode}.peakFrequency,'nbins',1000); %,'smooth',5
	axes(a);xlabel(['r=' num2str(stats{electrode}.durationFrequency.rho(1,2)) ' p=' num2str(stats{electrode}.durationFrequency.p(1,2))]);ylabel('Frequency vs Duration');
    line([prctile(data{electrode}.duration,99.5) prctile(data{electrode}.duration,99.5)],[100 200],'color','k')

	subplot(2,3,6);a = gca;
	PlotDistribution2(data{electrode}.duration,data{electrode}.peakAmplitude,'nbins',1000); %,'smooth',5
	axes(a);xlabel(['r=' num2str(stats{electrode}.durationAmplitude.rho(1,2)) ' p=' num2str(stats{electrode}.durationAmplitude.p(1,2))]);ylabel('Amplitude vs Duration');
    
    if ~isempty(foldername)
        saveas(gcf,['SummaryFigures\Ripples Stats2_',foldername,'.png']);
    else
        saveas(gcf,'SummaryFigures\Ripples Stats2.png');
    end
end



% Pairwise stats
% deal with this later...
% if nElectrodes == 2,
% 
% 		% Match ripple pairs
% 		ripplePairs = MatchPairs(ripples{1}(:,2),ripples{2}(:,2),'error',0.01); % within 10 ms
% 		% Ripples on LFP 1 only
% 		only1 = isnan(ripplePairs(:,2));
% 		% Ripples on LFP 2 only
% 		only2 = isnan(ripplePairs(:,1));
% 		% Ripples on both electrodes
% 		both = ~only1 & ~only2;
% 		% Peak time on LFP 1 when there are ripples on both electrodes
% 		ripplePairTimes = ripplePairs(both,1);
% 		% For each ripple on LFP 1, whether there are ripples on both electrodes
% 		ripplesOn1AlsoOn2 = ~isnan(ripplePairs(~only2,2));
% 		% For each ripple on LFP 2, whether there are ripples on both electrodes
% 		ripplesOn2AlsoOn1 = ~isnan(ripplePairs(~only1,1));
% 		% For each ripple on LFP 1, whether there are no ripples on LFP 2
% 		ripplesOn1NotOn2 = isnan(ripplePairs(~only2,2));
% 		% For each ripple on LFP 2, whether there are no ripples on LFP 1
% 		ripplesOn2NotOn1 = isnan(ripplePairs(~only1,1));
% 
% %  		% Compute phase shift map
% %  		size(maps{1}.phase)
% %  		[p1,i1] = Sync(maps{1}.phase,ripplePairTimes,'durations',durations);
% %  		[p2,i2] = Sync(maps{2}.phase,ripplePairTimes,'durations',durations);
% %  		shift(:,1) = p1(:,1);
% %  		shift(:,2) = p1(:,2)-p2(:,2);
% %  		phaseShiftMap = SyncMap(shift,i1,'durations',durations,'nbins',nBins,'smooth',0);
% %  		phaseShiftMap(phaseShiftMap<-pi) = phaseShiftMap(phaseShiftMap<-pi) + 2*pi;
% %  		phaseShiftMap(phaseShiftMap>pi) = phaseShiftMap(phaseShiftMap>pi) - 2*pi;
% %
% %  		% Determine average phase shift using center bins
% %  		meanCenterPhaseShift = mean(phaseShiftMap(:,centerBins),2);
% %  		[unused,order] = sortrows(meanCenterPhaseShift);
% 
% 		% Plot phase difference, frequency and amplitude (1 vs 2)
% 		f = figure;set(f,'name','Pairwise Stats');
% 
% %  		subplot(2,2,1);PlotColorMap(phaseShiftMap,1,'bar','on','type','circular');xlabel('Phase Shift');
% %  		hold on;plot([1 1]*size(phaseShiftMap,2)/2,ylim,'w');
% %
% %  		subplot(2,2,2);PlotColorMap(phaseShiftMap(order,:),1,'bar','on','type','circular');xlabel('Phase Shift');
% %  		hold on;plot([1 1]*size(phaseShiftMap,2)/2,ylim,'w');
% 
% 		subplot(2,3,4);a = gca;
% 		[rho,p] = corrcoef(data{1}.peakFrequency(ripplesOn1AlsoOn2),data{2}.peakFrequency(ripplesOn2AlsoOn1));
% 		PlotDistribution2(data{1}.peakFrequency(ripplesOn1AlsoOn2),data{2}.peakFrequency(ripplesOn2AlsoOn1),'nbins',1000,'smooth',50);
% 		axes(a);xlabel(['r=' num2str(rho(1,2)) ' p=' num2str(p(1,2))]);ylabel('Frequency (1 vs 2)');
% 
% 		subplot(2,3,5);a = gca;
% 		[rho,p] = corrcoef(data{1}.peakAmplitude(ripplesOn1AlsoOn2),data{2}.peakAmplitude(ripplesOn2AlsoOn1));
% 		PlotDistribution2(data{1}.peakAmplitude(ripplesOn1AlsoOn2),data{2}.peakAmplitude(ripplesOn2AlsoOn1),'nbins',1000,'smooth',50);
% 		axes(a);xlabel(['r=' num2str(rho(1,2)) ' p=' num2str(p(1,2))]);ylabel('Amplitude (1 vs 2)');
% 
% 		subplot(2,3,6);a = gca;
% 		rippleDurations1 = data{1}.duration(ripplesOn1AlsoOn2);
% 		rippleDurations2 = data{2}.duration(ripplesOn2AlsoOn1);
% 		discard = rippleDurations1>0.05|rippleDurations2>0.05;
% 		rippleDurations1 = rippleDurations1(~discard);
% 		rippleDurations2 = rippleDurations2(~discard);
% 		[rho,p] = corrcoef(rippleDurations1,rippleDurations2);
% 		PlotDistribution2(rippleDurations1,rippleDurations2,'nbins',1000,'smooth',50);
% 		axes(a);xlabel(['r=' num2str(rho(1,2)) ' p=' num2str(p(1,2))]);ylabel('Duration (1 vs 2)');
% 		% Plot info for unilateral vs bilateral ripples
% 		f = figure;set(f,'name','Unilateral vs Bilateral Stats');
% 
% 		subplot(2,2,1);
% 		S_boxplot([data{1}.peakFrequency(ripplesOn1NotOn2);data{1}.peakFrequency(ripplesOn1AlsoOn2);data{2}.peakFrequency(ripplesOn2NotOn1);data{2}.peakFrequency(ripplesOn2AlsoOn1)],[ones(sum(ripplesOn1NotOn2),1);2*ones(sum(ripplesOn1AlsoOn2),1);3*ones(sum(ripplesOn2NotOn1),1);4*ones(sum(ripplesOn2AlsoOn1),1)]);
% 		p1 = S_ranksum(data{1}.peakFrequency(ripplesOn1NotOn2),data{1}.peakFrequency(ripplesOn1AlsoOn2));
% 		p2 = S_ranksum(data{2}.peakFrequency(ripplesOn2NotOn1),data{2}.peakFrequency(ripplesOn2AlsoOn1));
% 		xlabel(['p=' num2str(p1) ' N=' num2str(sum(ripplesOn1NotOn2)) ',' num2str(sum(ripplesOn1AlsoOn2)) ' p=' num2str(p2) ' N=' num2str(sum(ripplesOn2NotOn1)) ',' num2str(sum(ripplesOn2AlsoOn1))]);
% 		ylabel('Frequency (1-2,1+2,2-1,2+1)');
% 
% 		subplot(2,2,2);
% 		S_boxplot([data{1}.peakAmplitude(ripplesOn1NotOn2);data{1}.peakAmplitude(ripplesOn1AlsoOn2);data{2}.peakAmplitude(ripplesOn2NotOn1);data{2}.peakAmplitude(ripplesOn2AlsoOn1)],[ones(sum(ripplesOn1NotOn2),1);2*ones(sum(ripplesOn1AlsoOn2),1);3*ones(sum(ripplesOn2NotOn1),1);4*ones(sum(ripplesOn2AlsoOn1),1)]);
% 		p1 = S_ranksum(data{1}.peakAmplitude(ripplesOn1NotOn2),data{1}.peakAmplitude(ripplesOn1AlsoOn2));
% 		p2 = S_ranksum(data{2}.peakAmplitude(ripplesOn2NotOn1),data{2}.peakAmplitude(ripplesOn2AlsoOn1));
% 		xlabel(['p=' num2str(p1) ' N=' num2str(sum(ripplesOn1NotOn2)) ',' num2str(sum(ripplesOn1AlsoOn2)) ' p=' num2str(p2) ' N=' num2str(sum(ripplesOn2NotOn1)) ',' num2str(sum(ripplesOn2AlsoOn1))]);
% 		ylabel('Amplitude (1-2,1+2,2-1,2+1)');
% 
% 		subplot(2,2,3);
% 		S_boxplot([data{1}.duration(ripplesOn1NotOn2);data{1}.duration(ripplesOn1AlsoOn2);data{2}.duration(ripplesOn2NotOn1);data{2}.duration(ripplesOn2AlsoOn1)],[ones(sum(ripplesOn1NotOn2),1);2*ones(sum(ripplesOn1AlsoOn2),1);3*ones(sum(ripplesOn2NotOn1),1);4*ones(sum(ripplesOn2AlsoOn1),1)]);
% 		p1 = S_ranksum(data{1}.duration(ripplesOn1NotOn2),data{1}.duration(ripplesOn1AlsoOn2));
% 		p2 = S_ranksum(data{2}.duration(ripplesOn2NotOn1),data{2}.duration(ripplesOn2AlsoOn1));
% 		xlabel(['p=' num2str(p1) ' N=' num2str(sum(ripplesOn1NotOn2)) ',' num2str(sum(ripplesOn1AlsoOn2)) ' p=' num2str(p2) ' N=' num2str(sum(ripplesOn2NotOn1)) ',' num2str(sum(ripplesOn2AlsoOn1))]);
% 		ylabel('Duration (1-2,1+2,2-1,2+1)');
% end
