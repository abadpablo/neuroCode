function [outputArg1,outputArg2] = bz_analyzeSession(varargin)
% Runs a list of descriptive analysis
%
% USAGE:
%   bz_analyzeSession(varargin)
%
% INPUT:
%   basepath              By default pwd
%   listOfAnalysis        Summary analysis that will be computed(as a cell
%                           with strings). Default 'all'. Possible values:
%                           'spikes','analogPulses','digitalPulses','ripples',
%                           'tMazeBehaviour','linearMazeBehaviour','OpenFieldBehaviour',
%                           'YMazeBehaviour','thetaModulation'
% exclude                 Cell with strings with the list of analysis to
%                         exluce
%
% (specific analysis options)
%  exlcludeShanks          Default []
%  getWaveformsFromDat     From 'spikes' summary, default false
%  analogChannelsList      Array of Channel to perform 'analogPulses' psth
%                           and csd. Default 'all'
%  digitalChannelList      Array of channel to perform 'digitalPulses' psth
%                           and csd. Default 'all'.
% analyzeSubSessions       Default false. Indicate if the analysis will be
%                           performed in individual subsessions instead of the whole session
%
% Pablo Abad 2021. Adapted from Manu Valero (computeSessionSummary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'listOfAnalysis','all',@iscellstr);
addParameter(p,'exclude',[],@iscellstr);
addParameter(p,'excludeShanks',[],@isnumeric);
addParameter(p,'getWaveformsFromDat',true,@islogical);
addParameter(p,'analogChannelsList','all',@isnumeric);
addParameter(p,'digitalChannelsList','all',@isnumeric);
addParameter(p,'analyzeSubSessions',false,@islogical);
addParameter(p,'showWaveforms',true,@islogical);
addParameter(p,'forceReloadRipples',false,@islogical);
addParameter(p,'diffLFPs',false,@islogical);
addParameter(p,'cellClassification',false,@islogical);
addParameter(p,'thetaFreq',[4 12],@isnumeric);
addParameter(p,'sgFreq',[30 65], @isnumeric);
addParameter(p,'hgFreq',[66 130],@isnumeric);
addParameter(p,'hfoFreq',[150 185], @isnumeric);
addParameter(p,'showFig',true,@isnumeric);
addParameter(p,'pathExcel','C:\Projects\GLUN3',@isstr);
addParameter(p,'nameExcel',[],@isstr);
addParameter(p,'selectedRippleChannel',[],@isnumeric);
addParameter(p,'selectedSWChannel',[],@isnumeric);
addParameter(p,'pass',5,@isnumeric);

parse(p,varargin{:});

basepath = p.Results.basepath;
listOfAnalysis = p.Results.listOfAnalysis;
exclude = p.Results.exclude;
excludeShanks = p.Results.excludeShanks;
getWaveformsFromDat = p.Results.getWaveformsFromDat;
analogChannelsList = p.Results.analogChannelsList;
digitalChannelsList = p.Results.digitalChannelsList;
analyzeSubSessions = p.Results.analyzeSubSessions;
showWaveforms = p.Results.showWaveforms;
forceReloadRipples = p.Results.forceReloadRipples;
diffLFPs = p.Results.diffLFPs;
thetaFreq = p.Results.thetaFreq;
sgFreq = p.Results.sgFreq;
hgFreq = p.Results.hgFreq;
hfoFreq = p.Results.hfoFreq;
showFig = p.Results.showFig;
pathExcel = p.Results.pathExcel;
nameExcel = p.Results.nameExcel;
selectedRippleChannel = p.Results.selectedRippleChannel;
selectedSWChannel = p.Results.selectedSWChannel;
pass = p.Results.pass*60;

prevPath = pwd;
cd(basepath);

if ischar(listOfAnalysis) && strcmpi(listOfAnalysis,'all')
    listOfAnalysis = {'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','newProtocol','excel','plotPlaceFields'};

end

if ~isempty(exclude)
    listOfAnalysis(ismember(listOfAnalysis,exclude)) = [];
end

session = loadSession;
excludeChannels = [];

for ii=1:length(excludeShanks)
    excludeChannels = [excludeChannels session.extracellular.electrodeGroups.channels{excludeShanks(ii)}];
end

mkdir('SummaryFigures'); % create folder
mkdir('PhaseModulationFig')
mkdir('lfpAnalysisFigures')
close all

%% 1 - Spikes Summary

if any(ismember(listOfAnalysis,'spikes'))
    try
          
        disp('Spike-waveform, ACG and cluster location...');
        spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat,'forceReload',false,'showWaveforms',showWaveforms);

        % plot spikes summary
        disp('Plotting spikes summary...');
        if exist('chanMap.mat','file')
            load('chanMap.mat','xcoords','ycoords','kcoords');
            xcoords = xcoords - min(xcoords); xcoords = xcoords/max(xcoords);
            ycoords = ycoords - min(ycoords); ycoords = ycoords/max(ycoords); 
            xcoords(excludeChannels) = [];
            ycoords(excludeChannels) = [];
            kcoords(excludeChannels) = [];
            
        else
            xcoords = NaN;
            ycoords = NaN;
            kcoords = NaN;
        end
        ccg=CCG(spikes.times,[],'binSize',0.001,'duration',0.08);
        dur = 0.08;
        xt = linspace(-dur/2*1000,dur/2*1000,size(ccg,1));
        figure
        set(gcf,'Position',[100 -100 2500 1200])
        for jj = 1:size(spikes.UID,2)
            fprintf(' **Spikes from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
            subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
            area(xt,ccg(:,jj,jj),'LineStyle','none');
            try 
                xt2 = linspace(dur/3*1000+5,dur/2*1000+5+80,length(spikes.filtWaveform{jj})); % waveform
                spk = spikes.filtWaveform{jj} - min(spikes.filtWaveform{jj}); spk = spk/max(spk) * max(ccg(:,jj,jj));
                hold on
                    plot(xt2,spk,'color','.8 .2 .2');

                plot((xcoords*30) + dur/2*1000+5+60, ycoords*max(ccg(:,jj,jj)),'.','color',[.8 .8 .8],'MarkerSize',5); % plotting xml
                plot((xcoords(spikes.maxWaveformCh1(jj))*30) + dur/2*1000+5+60,...
                    ycoords(spikes.maxWaveformCh1(jj))*max(ccg(:,jj,jj)),'.','color',[.1 .1 .1],'MarkerSize',10); % plotting xml
            end
            title(num2str(jj),'FontWeight','normal','FontSize',10);
            if jj == 1
                ylabel('Counts/ norm');
            elseif jj == size(spikes.UID,2)
                xlabel('Time (100 ms /1.5ms)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,'SummaryFigures\spikes AutoCorr.png');

        win = [-0.3 0.3];
        disp('Plotting CCG...');
        figure;
        set(gcf,'Position',[100 -100 2500 1200])
        [allCcg, t] = CCG(spikes.times,[],'binSize',0.005,'duration',0.6);
        indCell = [1:size(allCcg,2)];
        for jj = 1:size(spikes.UID,2)
            fprintf(' **CCG from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
            subplot(7,ceil(size(spikes.UID,2)/7),jj);
            cc = zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2); % get crosscorr
            imagesc(t,1:max(indCell)-1,cc)
            set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
            hold on
            zmean = mean(zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2));
            zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
            plot(t, zmean,'k','LineWidth',2);
            xlim([win(1) win(2)]); ylim([0 max(indCell)-1]);
            title(num2str(jj),'FontWeight','normal','FontSize',10);

            if jj == 1
                ylabel('Cell');
            elseif jj == size(spikes.UID,2)
                xlabel('Time (s)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,'SummaryFigures\spikes CrossCorr.png'); 
    catch
        warning('Error on Spike-waveform, autocorrelogram and cluster location! ');
    end
end


%% Compute CrossCorrelation between different areas and between same areas
% if any(ismember(listOfAnalysis,'meanCrossCorr'))
%     spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat);
%     sessionInfo = bz_getSessionInfo(pwd,'noPrompts',true);
%     if analyzeSubSessions
%         try
%             if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
%                 file = dir([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')]);
%                 load(file.name);
%             end
%             for ii=1:size(MergePoints.foldernames,2)
%                 timestamps = MergePoints.timestamps(ii,:);
%                 foldername = MergePoints.foldernames{ii};
%                 meanCrossCorrRegions.(foldername) = computeMeanCrossCorrelationRegions(spikes,'timestamps',timestamps,'foldername',foldername);
%                 
% %                 meanCrossCorr.(foldername) = computeMeanCrossCorrelation(spikes,'timestamps',timestamps,'foldername',foldername);
%             end
%             % Plot all subsessions together
%             figure,
%             set(gcf,'Position',get(0,'ScreenSize'))
% 
%             for i=1:size(MergePoints.foldernames,2)
%                 foldername = MergePoints.foldernames{i};
%                 plot(meanCrossCorrRegions.(foldername).mean_crosscorr)
%                 hold on
%                 time = strsplit(foldername,'_');
%                 time_{i} = time{end};
%             end
%             legend(time_)
%             saveas(gcf,'SummaryFigures\meanCrossCorr_AllSubSessions.png')
% 
%             % Save Output
%             save(fullfile(basepath, [sessionInfo.session.name, '.meanCrossCorr.SubSession.cellinfo.mat']),'meanCrossCorr')
% 
%         catch
%             disp('Not possible to run Mean CrossCorr for SubSessions...')
%         end
%     end
% end

%% 2 - Compute CrossCorrelation in different epochs and for subSessions
if any(ismember(listOfAnalysis,'meanCrossCorr'))
    spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat);
    sessionInfo = bz_getSessionInfo(pwd,'noPrompts',true);
    if analyzeSubSessions
        try
            if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
                file = dir([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')]);
                load(file.name);
            end
            for ii=1:size(MergePoints.foldernames,2)
                timestamps = MergePoints.timestamps(ii,:);
                foldername = MergePoints.foldernames{ii};
                meanCrossCorr.(foldername) = computeMeanCrossCorrelation(spikes,'timestamps',timestamps,'foldername',foldername);
            end
            % Plot all subsessions together
            figure,
            set(gcf,'Position',get(0,'ScreenSize'))

            for i=1:size(MergePoints.foldernames,2)
                foldername = MergePoints.foldernames{i};
                plot(meanCrossCorr.(foldername).mean_crosscorr)
                hold on
                time = strsplit(foldername,'_');
                time_{i} = time{end};
            end
            legend(time_)
            saveas(gcf,'SummaryFigures\meanCrossCorr_AllSubSessions.png')

            % Save Output
            save(fullfile(basepath, [sessionInfo.session.name, '.meanCrossCorr.SubSession.cellinfo.mat']),'meanCrossCorr')

        catch
            disp('Not possible to run Mean CrossCorr for SubSessions...')
        end
    end
    % Compute for all Recording anyways
    meanCrossCorrWS = computeMeanCrossCorrelation(spikes);
    save(fullfile(basepath, [sessionInfo.session.name, '.meanCrossCorr.cellinfo.mat']),'meanCrossCorrWS')
end

%% Compute rateRemapping

%% Digital Pulses
if any(ismember(listOfAnalysis,'digitalPulses'))
    try
        disp('Computing digital Pulses')
        sessionInfo = bz_getSessionInfo(basepath);
        digitalIn = getDigitalIn('all');
        spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat,'forceReload',false,'showWaveforms',showWaveforms);
        if ischar(digitalChannelsList) && strcmpi(digitalChannelsList,'all')
            digitalChannelsList = 1:length(digitalIn.timestampsOn);
            digitalChannelsList(1:2) = [];
        end
        
        for mm = 1:length(digitalChannelsList)
            fprintf('Stimulus %3.i of %3.i \n',mm, length(digitalChannelsList)); %\n
            st = digitalIn.timestampsOn{digitalChannelsList(mm)};
            % CSD
            shanks = sessionInfo.AnatGrps;
            shanks(excludeShanks) = [];
            figure
            set(gcf,'Position',[100 100 1400 600])
            for jj = 1:length(shanks)
                lfp = bz_GetLFP(shanks(jj).Channels,'noPrompts', true);
                twin = 0.02;
                [csd,lfpAvg] = bz_eventCSD(lfp,st,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
                taxis = linspace(-twin,twin,size(csd.data,1));
                cmax = max(max(csd.data)); 
                subplot(1,size(sessionInfo.AnatGrps,2),jj);
                contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('DigitalCh #',num2str(digitalChannelsList(mm))),'FontWeight','normal'); 
                colormap jet; try caxis([-cmax cmax]); end
                hold on
                for kk = 1:size(lfpAvg.data,2)
                    plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
                end
            end
            saveas(gcf,['SummaryFigures\digitalPulsesCSD_ch',num2str(digitalChannelsList(mm)), '.png']);

            % PSTH
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            
            if ~isempty(st)
                win = [-0.1 0.5];
                if length(st) > 1000
                    st = randsample(st, 1000);
                    st = sort(st);
                end

                disp('Plotting spikes raster and psth...');
                spikeResponse = [];
                [stccg, t] = CCG([spikes.times st'],[],'binSize',0.005,'duration',1);

                for jj = 1:size(spikes.UID,2)
                    fprintf(' **Pulses from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
                    rast_x = []; rast_y = [];
                    for kk = 1:length(st)
                        temp_rast = spikes.times{jj} - st(kk);
                        temp_rast = temp_rast(temp_rast>win(1) & temp_rast<win(2));
                        rast_x = [rast_x temp_rast'];
                        rast_y = [rast_y kk*ones(size(temp_rast))'];
                    end

                    spikeResponse = [spikeResponse; zscore(squeeze(stccg(:,end,jj)))'];
                    subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                    plot(rast_x, rast_y,'.','MarkerSize',1)
                    hold on
                    plot(t(t>win(1) & t<win(2)), stccg(t>win(1) & t<win(2),end,jj) * kk/max(stccg(:,end,jj))/2,'k','LineWidth',2);
                    xlim([win(1) win(2)]); ylim([0 kk]);
                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                    if jj == 1
                        ylabel('Trial');
                    elseif jj == size(spikes.UID,2)
                        xlabel('Time (s)');
                    else
                        set(gca,'YTick',[],'XTick',[]);
                    end
                end
            end
            saveas(gcf,['SummaryFigures\digitalPulsesRaster_ch',num2str(digitalChannelsList(mm)) ,'.png']); 
            
            figure
            if ~isempty(st)
                imagesc([t(1) t(end)],[1 size(spikeResponse,1)], spikeResponse); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([-.1 t(end)]);
            end
            saveas(gcf,['SummaryFigures\digitalPulsesPsth_',num2str(digitalChannelsList(mm)) ,'ch.png']); 
        end  
    catch
        disp('It is not possible to run digitalPulses')
    end
end
%% 2 - Ripples CSD and PSTH
if any(ismember(listOfAnalysis,'ripples'))
    disp('Ripples CSD and PSTH...');
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts',true);
%     lfp = bz_GetLFP(sessionInfo.channels,'noPrompts',true);
    spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat,'forceReload',false,'showWaveforms',showWaveforms);
    if isempty(dir('*ripples.events.mat')) || forceReloadRipples
        if analyzeSubSessions  
            try
                if ~isempty(selectedRippleChannel)
                    rippleChannel = computeRippleChannel('discardShanks',excludeShanks,'saveMat',true,'selectedRippleChannel',selectedRippleChannel,'selectedSWChannel',selectedSWChannel,'force',forceReloadRipples);
                else
                    rippleChannel = computeRippleChannel('discardShanks',excludeShanks,'saveMat',true);
                end
                if rippleChannel.Ripple_Channel == rippleChannel.Sharpwave_Channel
                    rippleChannel.Sharpwave_Channel = rippleChannel.Sharpwave_Channel +1;
                    rippleChannels = rippleChannel;
                    save([sessionInfo.FileName,'.channelInfo.ripples.mat'],'rippleChannels')
                    clear rippleChannels
                end
                shanks = sessionInfo.AnatGrps;
                shanks(excludeShanks) = [];
                passband = [130 200];
                lfp = bz_GetLFP(sessionInfo.channels,'noPrompts',true);
                filtered = bz_Filter(lfp,'channels',rippleChannel.Ripple_Channel,'filter','butter','passband',passband,'order',3);
                ripples_WS = bz_FindRipples(basepath,rippleChannel.Ripple_Channel);
                [maps_WS,data_WS,stats_WS] = bz_RippleStats(filtered.data,lfp.timestamps,ripples_WS);
                ripples_WS.maps = maps_WS;
                ripples_WS.data = data_WS;
                ripples_WS.stats = stats_WS;
                bz_PlotRippleStats_abad(maps_WS,data_WS,stats_WS,'ripple',ripples_WS)
                
                
                disp('Ripples CSD and PSTH for SubSession...')
            if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
                load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
                count = 1;
                for ii=1:size(MergePoints.foldernames,2)
                    timestamps = MergePoints.timestamps(ii,:);
                    foldername = MergePoints.foldernames{ii};
                    lfp = bz_GetLFP(sessionInfo.channels,'restrict',timestamps);        
                    rippleChannels{ii} = computeRippleChannel('discardShanks',excludeShanks,'timestamps_subSession',timestamps,'foldername',foldername,'saveMat',false);
                    
                    if ~isempty(rippleChannel.Ripple_Channel)
%                         ripples{ii} = bz_FindRipples(basepath,rippleChannels{ii}.Ripple_Channel,'timestamps_subSession',timestamps,'foldername',foldername);
                        ripples{ii} = bz_FindRipples(basepath,rippleChannel.Ripple_Channel,'timestamps_subSession',timestamps,'foldername',foldername);
                        if ~isempty(ripples{ii}.timestamps) && length(ripples{ii}.timestamps) >= 5
                            shanks = sessionInfo.AnatGrps;
                            shanks(excludeShanks) = [];
                            passband = [130 200];
    %                         filtered = bz_Filter(lfp,'channels',rippleChannels{ii}.Ripple_Channel,'filter','butter','passband',passband,'order',3);   
                            filtered = bz_Filter(lfp,'channels',rippleChannel.Ripple_Channel,'filter','butter','passband',passband,'order',3);
%                             ripples{ii} = bz_FindRipples(basepath,rippleChannel.Ripple_Channel,'timestamps_subSession',timestamps,'foldername',foldername);
                            [maps,data,stats] = bz_RippleStats(filtered.data,lfp.timestamps,ripples{ii});
                            ripples{ii}.maps = maps;
                            ripples{ii}.data = data;
                            ripples{ii}.stats = stats;
                            ripples{ii}.foldername = foldername;
                            bz_PlotRippleStats_abad(maps,data,stats,'foldername',foldername)

                            % CSD
                            twin = 0.1;
                            evs = ripples{ii}.peaks;
                            figure
                            set(gcf,'Position',[100 100 1400 600])
                            for jj = 1:size(shanks,2)
                                lfp = bz_GetLFP(shanks(jj).Channels,'noPrompts', true);
                                [csd,lfpAvg] = bz_eventCSD(lfp,evs,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
                                taxis = linspace(-twin,twin,size(csd.data,1));
                                cmax = max(max(csd.data)); 
                                subplot(1,size(shanks,2),jj);
                                contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                                set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('RIPPLES, Shank #',num2str(jj)),'FontWeight','normal'); 
                                colormap jet; caxis([-cmax cmax]);
                                hold on
                                for kk = 1:size(lfpAvg.data,2)
                                    plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
                                end
                            end

                            saveas(gcf,['SummaryFigures\',strcat('ripplesCSD.',foldername),'.png']);  

                            % PSTH
                            st = ripples{ii}.peaks;
                            spikeResponse = [];
                            win = [-0.2 0.2];
                            figure
                            set(gcf,'Position',[100 -100 2500 1200])

                            for jj = 1:size(spikes.UID,2)
                                fprintf(' **Ripple from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
                                rast_x = []; rast_y = [];
                                for kk = 1:length(st)
                                    temp_rast = spikes.times{jj} - st(kk);
                                    temp_rast = temp_rast(temp_rast>win(1) & temp_rast<win(2));
                                    rast_x = [rast_x temp_rast'];
                                    rast_y = [rast_y kk*ones(size(temp_rast))'];
                                end
                                [stccg, t] = CCG({spikes.times{jj} st},[],'binSize',0.005,'duration',1);
                                spikeResponse = [spikeResponse; zscore(squeeze(stccg(:,end,1:end-1)))'];
                                subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                                plot(rast_x, rast_y,'.','MarkerSize',1)
                                hold on
                                plot(t(t>win(1) & t<win(2)), stccg(t>win(1) & t<win(2),2,1) * kk/max(stccg(:,2,1))/2,'k','LineWidth',2);
                                xlim([win(1) win(2)]); ylim([0 kk]);
                                title(num2str(jj),'FontWeight','normal','FontSize',10);

                                if jj == 1
                                    ylabel('Trial');
                                elseif jj == size(spikes.UID,2)
                                    xlabel('Time (s)');
                                else
                                    set(gca,'YTick',[],'XTick',[]);
                                end
                            end
                            saveas(gcf,['SummaryFigures\',strcat('ripplesRaster.',foldername),'.png']); 

                            figure
                            imagesc([t(1) t(end)],[1 size(spikeResponse,2)], spikeResponse); caxis([-3 3]); colormap(jet);
                            xlim([-.2 .2]); set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells');
                            saveas(gcf,['SummaryFigures\',strcat('ripplesPsth.',foldername),'.png']); title('Ripples');
                        else
                            ripples{ii}.maps = [];
                            ripples{ii}.data = [];
                            ripples{ii}.stats = [];
                            ripples{ii}.foldername = foldername;
                        end
                    end
                end
                save(fullfile(basepath, [sessionInfo.session.name, '.ripples.SubSession.events.mat']),'ripples')
                % Maybe saving also the rippleChannels the same way
                save(fullfile(basepath,[sessionInfo.session.name,'.channelInfo.SubSession.ripples.mat']),'rippleChannels')
            end
            catch
                warning('Error on CSD and PSTH from ripples SubSessions !');
            end            
        else  
            try
                disp('Ripples CSD and PSTH for whole Session ...')
                
                rippleChannels = computeRippleChannel('discardShanks',excludeShanks);
                ripples = bz_FindRipples(basepath,rippleChannels.Ripple_Channel,'saveMat',true,'show','');

                shanks = sessionInfo.AnatGrps;
                shanks(excludeShanks) = [];
                passband = [130 200];
                lfp = bz_GetLFP(sessionInfo.channels,'noPrompts',true);
                filtered = bz_Filter(lfp,'channels',rippleChannels.Ripple_Channel,'filter','butter','passband',passband,'order',3);
                [maps,data,stats] = bz_RippleStats(filtered.data,lfp.timestamps,ripples);
                ripples.maps = maps;
                ripples.data = data;
                ripples.stats = stats;
                
                % CSD
                twin = 0.1;
                evs = ripples.peaks;
                figure
                set(gcf,'Position',[100 100 1400 600])
                for jj = 1:size(shanks,2)
                    lfp = bz_GetLFP(shanks(jj).Channels,'noPrompts', true);
                    [csd,lfpAvg] = bz_eventCSD(lfp,evs,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
                    taxis = linspace(-twin,twin,size(csd.data,1));
                    cmax = max(max(csd.data)); 
                    subplot(1,size(shanks,2),jj);
                    contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                    set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('RIPPLES, Shank #',num2str(jj)),'FontWeight','normal'); 
                    colormap jet; caxis([-cmax cmax]);
                    hold on
                    for kk = 1:size(lfpAvg.data,2)
                        plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
                    end
                end
                    
                saveas(gcf,['SummaryFigures\ripplesCSD.png']);  
                    
                    % PSTH
                    st = ripples.peaks;
                    spikeResponse = [];
                    win = [-0.2 0.2];
                    figure
                    set(gcf,'Position',[100 -100 2500 1200])
                    
                    for jj = 1:size(spikes.UID,2)
                        fprintf(' **Ripple from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
                        rast_x = []; rast_y = [];
                        for kk = 1:length(st)
                            temp_rast = spikes.times{jj} - st(kk);
                            temp_rast = temp_rast(temp_rast>win(1) & temp_rast<win(2));
                            rast_x = [rast_x temp_rast'];
                            rast_y = [rast_y kk*ones(size(temp_rast))'];
                        end
                        [stccg, t] = CCG({spikes.times{jj} st},[],'binSize',0.005,'duration',1);
                        spikeResponse = [spikeResponse; zscore(squeeze(stccg(:,end,1:end-1)))'];
                        subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                        plot(rast_x, rast_y,'.','MarkerSize',1)
                        hold on
                        plot(t(t>win(1) & t<win(2)), stccg(t>win(1) & t<win(2),2,1) * kk/max(stccg(:,2,1))/2,'k','LineWidth',2);
                        xlim([win(1) win(2)]); ylim([0 kk]);
                        title(num2str(jj),'FontWeight','normal','FontSize',10);

                        if jj == 1
                            ylabel('Trial');
                        elseif jj == size(spikes.UID,2)
                            xlabel('Time (s)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gcf,['SummaryFigures\ripplesRaster.png']); 
                    
                    figure
                    imagesc([t(1) t(end)],[1 size(spikeResponse,2)], spikeResponse); caxis([-3 3]); colormap(jet);
                    xlim([-.2 .2]); set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells');
                    saveas(gcf,['SummaryFigures\ripplesPsth.png']); title('Ripples');
                    save(fullfile(basepath, [sessionInfo.session.name '.ripples.events.mat']),'ripples')
            catch
                warning('Error on CSD and PSTH from ripples Whole Session !');
            end   
        end
    end
end


%% Replay and Reactivation

% if any(ismember(listOfAnalysis,'replay'))
%     disp('Computing Replay and Reactivation...code to improve !')
%     sessionInfo = bz_getSessionInfo();
%     spikes = loadSpikes();
%     if ~isempty(dir([sessionInfo.FileName '.ripples.SubSession.events.mat']))
%         file = dir([sessionInfo.FileName,'.ripples.SubSession.events.mat']);
%         load(file.name);
%     end
%     
%     replay = bz_computeReplay(spikes,ripples);
%     
% end
%% POWER SPECTRUM PROFILE

if any(ismember(listOfAnalysis,'powerSpectrumProfile'))
    if analyzeSubSessions
        try
            [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts',true);
            if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
                load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
                count = 1;
                for i = 1:size(MergePoints.foldernames,2)
                    timestamps = MergePoints.timestamps(i,:);
                    foldername = MergePoints.foldernames{i};
                    disp(['Theta, gamma and HFO powerSpectrum for folder: ', foldername])
                    % Theta Profile
                    powerProfile_theta{i} = bz_PowerSpectrumProfile([thetaFreq(1) thetaFreq(2)], 'useparfor',false,'timestampsSubSession',timestamps,'foldername',foldername,'channels', sessionInfo.channels, 'showfig',true,'saveMat',false)
                    % Low Gamma Profile                  
                    powerProfile_sg{i} = bz_PowerSpectrumProfile([sgFreq(1) sgFreq(2)], 'useparfor',false,'timestampsSubSession',timestamps,'foldername',foldername,'channels',sessionInfo.channels,'showfig',true,'saveMat',false);
                    % HFO Profile
                    powerProfile_hg{i} = bz_PowerSpectrumProfile([hgFreq(1) hgFreq(2)], 'useparfor',false,'timestampsSubSession',timestamps,'foldername',foldername,'channels',sessionInfo.channels,'showfig',true,'saveMat',false);
                    % get channels of interest % max theta power above pyr layer
                    [~, a] = max(powerProfile_theta{i}.mean);
                    region{i}.CA1sp = powerProfile_theta{i}.channels(a);
                    if ~isempty('*.channelInfo.ripples.mat')
                        file = dir('*.channelInfo.ripples.mat');
                        load(file.name)
                    end
                    
                    for ii = 1:size(sessionInfo.AnatGrps,2)
%                         if ismember(region{i}.CA1sp, sessionInfo.AnatGrps(ii).Channels)
                        if ismember(rippleChannels.Ripple_Channel, sessionInfo.AnatGrps(ii).Channels)
                            
                            p_th{i} = powerProfile_theta{i}.mean(sessionInfo.AnatGrps(ii).Channels + 1);                
                            p_sg{i} = powerProfile_sg{i}.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                            p_hg{i} = powerProfile_hg{i}.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                            channels = sessionInfo.AnatGrps(ii).Channels;
                            [~,a] = max(p_th{i}(1:find(channels == region{i}.CA1sp)));
                            region{i}.CA1so = channels(a);
                            [~,a] = max(p_th{i}(find(channels == region{i}.CA1sp):end));
                            region{i}.CA1slm = channels(a - 1+ find(channels == region{i}.CA1sp));

                            figure,
                            hold on
                            plot(1:length(channels), zscore(p_th{i}));
                            plot(1:length(channels), zscore(p_hg{i}));
                            plot(1:length(channels), zscore(p_sg{i}));
                            xlim([1 length(channels)]);
                            set(gca,'XTick',1:length(channels),'XTickLabel',channels,'XTickLabelRotation',45);
                            ax = axis;
                            xlabel(strcat('Channels (neuroscope)-Shank',num2str(ii))); ylabel('power (z)');
%                             plot(find(channels == region{i}.CA1sp)*ones(2,1),ax([3 4]),'-k');
%                             plot(find(channels == region{i}.CA1so)*ones(2,1),ax([3 4]),'--k');
%                             plot(find(channels == region{i}.CA1slm)*ones(2,1),ax([3 4]),'-.k');
                            plot(find(channels == rippleChannels.Ripple_Channel)*ones(2,1),ax([3 4]),'-k');
                            legend('4-12Hz', 'hfo', '30-60','pyr','~or','~slm');
                            saveas(gcf,['SummaryFigures\regionDef_',foldername,'.png']);
                        end
                    end
                    % power profile of pyr channel of all session
%                     lfpT = bz_GetLFP(region{i}.CA1sp,'restrict',timestamps,'noPrompts',true);
                    lfpT = bz_GetLFP(rippleChannels.Ripple_Channel,'restrict',timestamps,'noPrompts',true);
                    params.Fs = lfpT.samplingRate; params.fpass = [2 200]; params.tapers = [3 5]; params.pad = 1;
                    [S,t,f] = mtspecgramc(single(lfpT.data),[2 1], params);
                    S = log10(S); % in Db
                    S_det= bsxfun(@minus,S,polyval(polyfit(f,mean(S,1),2),f)); % detrending        
                    figure;
                    subplot(3,4,[1 2 3 5 6 7])
                    imagesc(t,f,S_det',[-1.5 1.5]);
%                     set(gca,'XTick',[]); 
                    ylabel('Freqs');
                    subplot(3,4,[4 8]);
                    plot(mean(S,1),f);
                    set(gca,'YDir','reverse','YTick',[]); xlabel('Power');
                    saveas(gcf,['SummaryFigures\spectrogram_',foldername,'.png']);
                end
            end
            % save Power Spectrum Profile theta
            save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_theta{1}.processinginfo.params.frange(1)),'_',num2str(powerProfile_theta{1}.processinginfo.params.frange(2)),'.SubSession.channelinfo.mat'],'powerProfile_theta');                
            % save Power Spectrum Profile sg
            save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_sg{1}.processinginfo.params.frange(1)),'_',num2str(powerProfile_sg{1}.processinginfo.params.frange(2)),'.SubSession.channelinfo.mat'],'powerProfile_sg');
            % save Power Spectrum Profile hg
            save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_hg{1}.processinginfo.params.frange(1)),'_',num2str(powerProfile_hg{1}.processinginfo.params.frange(2)),'.SubSession.channelinfo.mat'],'powerProfile_hg');
            save([sessionInfo.FileName,'.region.mat'],'region');
        catch
            warning('It has not been possible to run theta and power Spectrum Profile for SubSessions...')
        end
    else
        try
           [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts',true);
            % Theta Profile
            powerProfile_theta = bz_PowerSpectrumProfile([thetaFreq(1) thetaFreq(2)], 'channels',sessionInfo.channels,'showfig',true,'saveMat',false);
            % Low Gamma Profile
            powerProfile_sg = bz_PowerSpectrumProfile([sgFreq(1) sgFreq(2)], 'channels',sessionInfo.channels,'showfig',true,'saveMat',false);
            % HFO Profile
            powerProfile_hg = bz_PowerSpectrumProfile([hgFreq(1) hgFreq(2)], 'channels', sessionInfo.channels,'showfig',true,'saveMat',false);
            % get channels of interest % max theta power above pyr layer
            [~, a] = max(powerProfile_theta.mean);
            region.CA1sp = powerProfile_theta.channels(a); 
            for ii = 1:size(sessionInfo.AnatGrps,2)
                if ismember(region.CA1sp, sessionInfo.AnatGrps(ii).Channels)
                    p_th = powerProfile_theta.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                    p_hg = powerProfile_hg.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                    p_sg = powerProfile_sg.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                    channels = sessionInfo.AnatGrps(ii).Channels;
                    [~,a] = max(p_th(1:find(channels == region.CA1sp)));
                    region.CA1so = channels(a);
                    [~,a] = max(p_th(find(channels == region.CA1sp):end));
                    region.CA1slm = channels(a - 1+ find(channels == region.CA1sp));

                    figure,
                    hold on
                    plot(1:length(channels), zscore(p_th));
                    plot(1:length(channels), zscore(p_hg));
                    plot(1:length(channels), zscore(p_sg));
                    xlim([1 length(channels)]);
                    set(gca,'XTick',1:length(channels),'XTickLabel',channels,'XTickLabelRotation',45);
                    ax = axis;
                    xlabel(strcat('Channels (neuroscope)-Shank',num2str(ii))); ylabel('power (z)');
                    plot(find(channels == region.CA1sp)*ones(2,1),ax([3 4]),'-k');
                    plot(find(channels == region.CA1so)*ones(2,1),ax([3 4]),'--k');
                    plot(find(channels == region.CA1slm)*ones(2,1),ax([3 4]),'-.k');
                    legend('4-12Hz', 'hfo', '30-60','pyr','~or','~slm');
                    saveas(gcf,'SummaryFigures\regionDef.png');
                end
            end
            save([sessionInfo.FileName,'.region.mat'],'region');
            
            % power profile of pyr channel of all session
            lfpT = bz_GetLFP(region.CA1sp,'noPrompts',true);
            params.Fs = lfpT.samplingRate; params.fpass = [2 200]; params.tapers = [3 5]; params.pad = 1;
            [S,t,f] = mtspecgramc(single(lfpT.data),[2 1], params);
            S = log10(S); % in Db
            S_det= bsxfun(@minus,S,polyval(polyfit(f,mean(S,1),2),f)); % detrending        
            figure;
            subplot(3,4,[1 2 3 5 6 7])
            imagesc(t,f,S_det',[-1.5 1.5]);
            set(gca,'XTick',[]); ylabel('Freqs');
            subplot(3,4,[4 8]);
            plot(mean(S,1),f);
            set(gca,'YDir','reverse','YTick',[]); xlabel('Power');
            saveas(gcf,'SummaryFigures\spectrogramAllSession.png');
            
            % save Power Spectrum Profile theta
            save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_theta.processinginfo.params.frange(1)),'_',num2str(powerProfile_theta.processinginfo.params.frange(2)),'.channelinfo.mat'],'powerProfile_theta');                
            % save Power Spectrum Profile sg
            save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_sg.processinginfo.params.frange(1)),'_',num2str(powerProfile_sg.processinginfo.params.frange(2)),'.channelinfo.mat'],'powerProfile_sg');
            % save Power Spectrum Profile hg
            save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_hg.processinginfo.params.frange(1)),'_',num2str(powerProfile_hg.processinginfo.params.frange(2)),'.channelinfo.mat'],'powerProfile_hg');
            save([sessionInfo.FileName,'.region.mat'],'region');            
        catch
            warning('It has not been possible to run theta and gamma mod code...')
        end
    end
end

%% 3 - THETA AND GAMMA PHASE MODULATION BY SPIKES and POWER SPECTRUM PROFILE

if any(ismember(listOfAnalysis,'thetaModulation'))    
    if analyzeSubSessions
        try
            disp('Theta, gamma and HFO modulation for SubSessions ...')
            [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts',true);
            spikes = loadSpikes('getWaveformsFromDat',false);
            if ~isempty('*.channelInfo.ripples.mat')
                file = dir('*.channelInfo.ripples.mat');
                load(file.name)
            end
            if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
                load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
                count = 1;
                for i = 1:size(MergePoints.foldernames,2)
                    timestamps = MergePoints.timestamps(i,:);
                    foldername = MergePoints.foldernames{i};
                    
                    % phase modulation by spikes
                    if diffLFPs
                        filename = dir('*PhaseLockingData.SubSession.diffLFPs.mat*'); 
                        if ~isempty(filename)
                            disp(['Phase Locking Data SubSession Different LFPs already detected ! Loading file: ', filename.name]);
                            load(filename.name)
                            PLD = PhaseLockingData;
                        else
                            for ii = 1:spikes.numcells
                                lfp = bz_GetLFP(spikes.maxWaveformCh1(ii)-1,'noPrompts',true);
                                PLD{i}{ii} = bz_PhaseModulationDifferentLFPs(spikes.times{ii},lfp,[thetaFreq(1) thetaFreq(2)], 'intervals',timestamps,'plotting',false,'method','wavelet');
                            end
                            PhaseLockingData{i} = PLD{i};

                        end
                        if showFig
                            disp('Theta modulation...')
                            figure,
                            set(gcf,'Position',[100 -100 2500 1200]);
                            for jj = 1:size(spikes.UID,2)
                                subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                                if ~isempty(PLD{i}{jj}.phasebins)
                                    area([PLD{i}{jj}.phasebins; PLD{i}{jj}.phasebins + pi*2],[PLD{i}{jj}.phasedistros; PLD{i}{jj}.phasedistros],'EdgeColor','none');
                                    hold on
                                    ax = axis;
                                    x = 0:.001:4*pi;
                                    y = cos(x);
                                    y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                                    h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                                    xlim([0 4*pi]);
                                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                                    if jj == 1
                                        ylabel('prob');
                                    elseif jj == size(spikes.UID,2)
                                        set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                        xlabel('phase (rad)');
                                    else
                                        set(gca,'YTick',[],'XTick',[]);
                                    end
                                else
                                    
                                end
                            end
                            saveas(gcf,['PhaseModulationFig\thetaModDiffLFPS_',foldername,'.png']);
                        end

                        % gamma modulation
                        disp('Gamma modulation...');
                        filename = dir('*PhaseLockingData_sg.diffLFPs.SubSession.cellinfo.mat*');
                        if ~isempty(filename)
                            disp(['Phase Locking Data_sg SubSession Different LFPs already detected ! Loading file: ', filename.name]);
                            load(filename.name)
                            PLD_sg = PhaseLockingData_sg;
                        else
                            for ii=1:spikes.numcells
                                lfp = bz_GetLFP(spikes.maxWaveformCh1(ii)-1,'noPrompts',true);
                                PLD_sg{i}{ii} = bz_PhaseModulationDifferentLFPs(spikes.times{ii},lfp,[sgFreq(1) sgFreq(2)],'intervals',timestamps,'plotting',false,'method','wavelet');                              
                            end
                            PhaseLockingData_sg{i} = PLD_sg{i};
                        end
                        
                        if showFig
                            figure
                            set(gcf,'Position',[100 -100 2500 1200]);
                            for jj = 1:size(spikes.UID,2)
                                if ~isempty(PLD_sg{i}{jj}.phasebins)
                                    subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                                    area([PLD_sg{i}{jj}.phasebins; PLD_sg{i}{jj}.phasebins + pi*2],[PLD_sg{i}{jj}.phasedistros; PLD_sg{i}{jj}.phasedistros],'EdgeColor','none');
                                    hold on
                                    ax = axis;
                                    x = 0:.001:4*pi;
                                    y = cos(x);
                                    y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                                    h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                                    xlim([0 4*pi]);
                                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                                    if jj == 1
                                        ylabel('prob');
                                    elseif jj == size(spikes.UID,2)
                                        set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                        xlabel('phase (rad)');
                                    else
                                        set(gca,'YTick',[],'XTick',[]);
                                    end
                                end
                            end
                            saveas(gcf,['PhaseModulationFig\sgModDiffLFPS_',foldername,'.png']);
                        end
                    else
%                         lfpT = bz_GetLFP(region{i}.CA1sp,'noPrompts',true);
%                         PLD{i} = bz_PhaseModulation(spikes,lfpT,[thetaFreq(1) thetaFreq(2)],'plotting',false,'method','wavelet','saveMat',false);
                        lfpT = bz_GetLFP(rippleChannels.Ripple_Channel,'noPrompts',true);
                        PLD{i} = bz_PhaseModulationAbad(spikes,lfpT,[thetaFreq(1) thetaFreq(2)],'intervals',timestamps,'plotting',false,'method','wavelet','saveMat',false);
%                         lfpT = bz_GetLFP(rippleChannels.Ripple_Channel,'restrict',timestamps,'noPrompts',true);
%                         PLD{i} = bz_PhaseModulationAbadv2(spikes,lfpT,[thetaFreq(1) thetaFreq(2)],'plotting',false,'method','wavelet','saveMat',false);
%                         PLD{i} = bz_SpikeLFP(spikes,lfpT,[thetaFreq(1) thetaFreq(2)], 'plotting',false,'method',wavelet,'saveMat',false);
                        if showFig
                            figure,
                            disp('Theta modulation...')
                            set(gcf,'Position',[100 -100 2500 1200]);
                            for jj = 1:size(spikes.UID,2)                                
                                subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                                if ~isempty(PLD{i}.phasebins)
                                    area([PLD{i}.phasebins; PLD{i}.phasebins + pi*2],[PLD{i}.phasedistros(:,jj); PLD{i}.phasedistros(:,jj)],'EdgeColor','none');
                                    hold on
                                    ax = axis;
                                    x = 0:.001:4*pi;
                                    y = cos(x);
                                    y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                                    h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                                    xlim([0 4*pi]);
                                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                                    if jj == 1
                                        ylabel('prob');
                                    elseif jj == size(spikes.UID,2)
                                        set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                        xlabel('phase (rad)');
                                    else
                                        set(gca,'YTick',[],'XTick',[]);
                                    end
                                end
                            end
                            saveas(gcf,['PhaseModulationFig\thetaMod_',foldername,'.png']);
                        end

                        % gamma modulation
%                         PLD_sg{i} = bz_PhaseModulation(spikes,lfpT,[sgFreq(1) sgFreq(2)],'intervals',timestamps,'plotting',false,'method','wavelet','saveMat',false);     
                        PLD_sg{i} = bz_PhaseModulationAbad(spikes,lfpT,[sgFreq(1) sgFreq(2)],'intervals',timestamps,'plotting',false,'method','wavelet','saveMat',false);
%                         PLD{i} = bz_SpikeLFP(spikes,lfpT,[thetaFreq(1) thetaFreq(2)], 'plotting',false,'method',wavelet,'saveMat',false);

                        disp('Gamma modulation...');
                        if showFig
                            figure
                            set(gcf,'Position',[100 -100 2500 1200]);
                            for jj = 1:size(spikes.UID,2)
                                subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                                if ~isempty(PLD{i}.phasebins)
                                    area([PLD_sg{i}.phasebins; PLD_sg{i}.phasebins + pi*2],[PLD_sg{i}.phasedistros(:,jj); PLD_sg{i}.phasedistros(:,jj)],'EdgeColor','none');
                                    hold on
                                    ax = axis;
                                    x = 0:.001:4*pi;
                                    y = cos(x);
                                    y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                                    h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                                    xlim([0 4*pi]);
                                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                                    if jj == 1
                                        ylabel('prob');
                                    elseif jj == size(spikes.UID,2)
                                        set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                        xlabel('phase (rad)');
                                    else
                                        set(gca,'YTick',[],'XTick',[]);
                                    end
                                end
                            end
                            saveas(gcf,['PhaseModulationFig\gammaMod_',foldername,'.png']);
                        end
                    end    
                end
                if ~diffLFPs
                    PhaseLockingData = PLD;
                    PhaseLockingData_sg = PLD_sg;
                end
                
                % save Power Spectrum Profile theta
%                 save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_theta{1}.processinginfo.params.frange(1)),'_',num2str(powerProfile_theta{1}.processinginfo.params.frange(2)),'.SubSession.channelinfo.mat'],'powerProfile_theta');                
                % save Power Spectrum Profile sg
%                 save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_sg{1}.processinginfo.params.frange(1)),'_',num2str(powerProfile_sg{1}.processinginfo.params.frange(2)),'.SubSession.channelinfo.mat'],'powerProfile_sg');
                % save Power Spectrum Profile hg
%                 save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_hg{1}.processinginfo.params.frange(1)),'_',num2str(powerProfile_hg{1}.processinginfo.params.frange(2)),'.SubSession.channelinfo.mat'],'powerProfile_hg');
%                 save([sessionInfo.FileName,'.region.mat'],'region');
                save([sessionInfo.session.name '.PhaseLockingData.SubSession.cellinfo.mat'],'PhaseLockingData');
                save([sessionInfo.session.name '.PhaseLockingData_sg.SubSession.cellinfo.mat'],'PhaseLockingData_sg');
                    
            end
        catch
            warning('It has not been possible to run theta and gamma mod code for SubSessions...')
        end
    else
        try
            disp('Theta, gamma and HFO modulation...')
             % phase modulation by spikes
            [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts',true);
            spikes = loadSpikes('getWaveformsFromDat',false);
            if ~isempty('*.channelInfo.ripples.mat')
                file = dir('*.channelInfo.ripples.mat');
                load(file.name)
            end
            if diffLFPs
                filename = dir('*PhaseLockingData.diffLFPs.cellinfo.mat*'); 
                if ~isempty(filename)
                    disp(['Phase Locking Data Different LFPs already detected ! Loading file: ', filename.name]);
                    load(filename.name)
                    PLD = PhaseLockingData;
                else
                    for ii = 1:spikes.numcells
                        lfp = bz_GetLFP(spikes.maxWaveformCh1(ii)-1,'noPrompts',true);
                        PLD{ii} = bz_PhaseModulationDifferentLFPs(spikes.times{ii},lfp,[thetaFreq(1) thetaFreq(2)], 'plotting',false,'method','wavelet');
                    end
                    PhaseLockingData = PLD;
                    save([sessionInfo.session.name '.PhaseLockingData.diffLFPs.cellinfo.mat'],'PhaseLockingData');
                end
                disp('Theta modulation...')
                figure,
                set(gcf,'Position',[100 -100 2500 1200]);
                for jj = 1:size(spikes.UID,2)
                    subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                    area([PLD{jj}.phasebins; PLD{jj}.phasebins + pi*2],[PLD{jj}.phasedistros; PLD{jj}.phasedistros],'EdgeColor','none');
                    hold on
                    ax = axis;
                    x = 0:.001:4*pi;
                    y = cos(x);
                    y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                    h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                    xlim([0 4*pi]);
                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                    if jj == 1
                        ylabel('prob');
                    elseif jj == size(spikes.UID,2)
                        set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                        xlabel('phase (rad)');
                    else
                        set(gca,'YTick',[],'XTick',[]);
                    end
                end
                saveas(gcf,'PhaseModulationFig\thetaModDiffLFPS.png');

                % gamma modulation
                disp('Gamma modulation...');
                filename = dir('*PhaseLockingData_sg.diffLFPs.cellinfo.mat*');
                if ~isempty(filename)
                    disp(['Phase Locking Data_sg Different LFPs already detected ! Loading file: ', filename.name]);
                    load(filename.name)
                    PLD_sg = PhaseLockingData_sg;
                else
                    for ii=1:spikes.numcells
                        lfp = bz_GetLFP(spikes.maxWaveformCh1(ii)-1,'noPrompts',true);
                        PLD_sg{ii} = bz_PhaseModulationDifferentLFPs(spikes.times{ii},lfp,[sgFreq(1) sgFreq(2)],'plotting',false,'method','wavelet');     
                        PhaseLockingData_sg = PLD_sg;
                        save([sessionInfo.session.name '.PhaseLockingData_sg.diffLFPs.cellinfo.mat'],'PhaseLockingData_sg');
                    end
                end

                figure
                set(gcf,'Position',[100 -100 2500 1200]);
                for jj = 1:size(spikes.UID,2)
                    subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                    area([PLD_sg{jj}.phasebins; PLD_sg{jj}.phasebins + pi*2],[PLD_sg{jj}.phasedistros; PLD_sg{jj}.phasedistros],'EdgeColor','none');
                    hold on
                    ax = axis;
                    x = 0:.001:4*pi;
                    y = cos(x);
                    y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                    h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                    xlim([0 4*pi]);
                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                    if jj == 1
                        ylabel('prob');
                    elseif jj == size(spikes.UID,2)
                        set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                        xlabel('phase (rad)');
                    else
                        set(gca,'YTick',[],'XTick',[]);
                    end
                end
                saveas(gcf,'PhaseModulationFig\sgModDiffLFPS.png');

            else
%                 lfpT = bz_GetLFP(region.CA1sp,'noPrompts',true);
                lfpT = bz_GetLFP(rippleChannels.Ripple_Channel,'noPrompts',true);
                PLD = bz_PhaseModulation(spikes,lfpT,[thetaFreq(1) thetaFreq(2)], 'plotting',false,'method','wavelet','saveMat',false);
                disp('Theta modlation...')
                set(gcf,'Position',[100 -100 2500 1200]);
                for jj = 1:size(spikes.UID,2)
                    subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                    area([PLD.phasebins; PLD.phasebins + pi*2],[PLD.phasedistros(:,jj); PLD.phasedistros(:,jj)],'EdgeColor','none');
                    hold on
                    ax = axis;
                    x = 0:.001:4*pi;
                    y = cos(x);
                    y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                    h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                    xlim([0 4*pi]);
                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                    if jj == 1
                        ylabel('prob');
                    elseif jj == size(spikes.UID,2)
                        set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                        xlabel('phase (rad)');
                    else
                        set(gca,'YTick',[],'XTick',[]);
                    end
                end
                saveas(gcf,'SummaryFigures\thetaMod.png');

                % gamma modulation
                PLD_sg = bz_PhaseModulation(spikes,lfpT,[sgFreq(1) sgFreq(2)],'plotting',false,'method','wavelet','saveMat',false);     
                disp('Gamma modulation...');
                figure
                set(gcf,'Position',[100 -100 2500 1200]);
                for jj = 1:size(spikes.UID,2)
                    subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                    area([PLD_sg.phasebins; PLD_sg.phasebins + pi*2],[PLD_sg.phasedistros(:,jj); PLD_sg.phasedistros(:,jj)],'EdgeColor','none');
                    hold on
                    ax = axis;
                    x = 0:.001:4*pi;
                    y = cos(x);
                    y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                    h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                    xlim([0 4*pi]);
                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                    if jj == 1
                        ylabel('prob');
                    elseif jj == size(spikes.UID,2)
                        set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                        xlabel('phase (rad)');
                    else
                        set(gca,'YTick',[],'XTick',[]);
                    end
                end
                saveas(gcf,'SummaryFigures\gammaMod.png');
            end
            % save Power Spectrum Profile theta
            save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_theta.processinginfo.params.frange(1)),'_',num2str(powerProfile_theta.processinginfo.params.frange(2)),'.channelinfo.mat'],'powerProfile_theta');                
            % save Power Spectrum Profile sg
            save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_sg.processinginfo.params.frange(1)),'_',num2str(powerProfile_sg.processinginfo.params.frange(2)),'.channelinfo.mat'],'powerProfile_sg');
            % save Power Spectrum Profile hg
            save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_hg.processinginfo.params.frange(1)),'_',num2str(powerProfile_hg.processinginfo.params.frange(2)),'.channelinfo.mat'],'powerProfile_hg');
            save([sessionInfo.FileName,'.region.mat'],'region');
            if ~diffLFPs
                PhaseLockingData = PLD;
                PhaseLockingData_sg = PLD_sg;
            end
            save([sessionInfo.session.name '.PhaseLockingData.cellinfo.mat'],'PhaseLockingData');
            save([sessionInfo.session.name '.PhaseLockingData_sg.cellinfo.mat'],'PhaseLockingData_sg');
            
        catch
            warning('It has not been possible to run theta and gamma mod code...')
        end
    end
end

%% 4- BEHAVIOUR

if any(ismember(listOfAnalysis,'behaviour'))   
    try
        disp('Computing behaviour variables...')
        tracking = getSessionTracking();
        behaviour = getSessionBehaviour_v2();
        
        % PLACE CELLS SUMMARY
        if any(ismember(listOfAnalysis,'placeCells'))
            try
                spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat,'forceReload',false,'showWaveforms',showWaveforms);
                firingMaps = bz_firingMapAvg(behaviour,spikes,'saveMat',true,'speedFilter',true,'periodicAnalysis',false,'spikeShuffling',true,'numRand',1000);     
            catch
                warning('It has not been possible to run Place Cells Analysis...')
            end
        end
        if any(ismember(listOfAnalysis,'plotLinearTrack'))
            placeFieldStats = bz_findPlaceFields1D('firingMaps',firingMaps);
        end
    catch
        warning('It has not been possible to run Behaviour code...')
    end
end

%% 5  BEHAVIOUR PERFORMANCE

if any(ismember(listOfAnalysis,'performance'))
    try
        disp('Computing Performance for Experimental Paradigm (YMaze/TMaze) ...');
        tracking = getSessionTracking();
        performance = getSessionPerformance('tracking',tracking);
        
    catch
        warning('It has not been possible to run Performance code...')
    
    end
end

%% DISTANCE BY EPOCHS

if any(ismember(listOfAnalysis,'distanceByEpochs'))
    tracking = getSessionTracking();
    distanceByEpochs = computeDistanceByEpochs(tracking);
end
%% 6 - SPIKE TRAIN ANALYSIS

if any(ismember(listOfAnalysis,'spikeTrain'))
    try
        disp('Computing Spike Train Analysis...')
        spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat,'forceReload',false,'showWaveforms',showWaveforms);
        spikeTrain = bz_SpikeTrain(spikes,'analyzeSubSessions',analyzeSubSessions,'showFigure',true);
    catch
        warning('It has not been possible to run Spike Train Analysis...')
    end
end
%% 7 - CELL EXPLORER

if any(ismember(listOfAnalysis,'CellExplorer'))
    try
        disp('Running CellExplorer code...');
        session = sessionTemplate(basepath,'showGUI',false);        
        validateSessionStruct(session);        
        cell_metrics = ProcessCellMetrics('session', session,'showGUI',true);
        cell_metrics = CellExplorer('metrics',cell_metrics); 
    catch
        disp('It is not possible to run CellExplorer code...')
    end
end
%% 8 - LFP ANALYSIS
cd(basepath)
if any(ismember(listOfAnalysis,'lfp_analysis'))
    mkdir('lfpAnalysisFigures')
%     analyzeSubSessions = true;
  if analyzeSubSessions   
      try
          disp('Computing lfp_analysis for SubSessions...')
          session = loadSession;
          sessionInfo = bz_getSessionInfo(pwd,'noPrompts',true);
          shanks = sessionInfo.AnatGrps;
          shanks(excludeShanks) = [];
          
          if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
              load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
              count = 1;
              
              for ii=1:size(MergePoints.foldernames,2)
                  timestamps = MergePoints.timestamps(ii,:);
                  foldername = MergePoints.foldernames{ii};
                  % Pick Channel of Ripples and SW
%                   if ~isempty(dir([basepath filesep '*.channelinfo.SubSession.ripples.mat']))
%                      disp('Channelinfo.ripples.mat file detected ! Loading file.')
%                      file = dir([basepath filesep '*channelinfo.SubSession.ripples.mat'])
%                      load(file.name);
%                   end
                  
                  if ~isempty(dir([basepath filesep '*.channelinfo.ripples.mat']))
                     disp('Channelinfo.ripples.mat file detected ! Loading file.')
                     file = dir([basepath filesep '*channelinfo.ripples.mat'])
                     load(file.name);
                  end
                  
%                   lfp1 = bz_GetLFP(rippleChannels{ii}.Ripple_Channel,'restrict',timestamps);
%                   lfp2 = bz_GetLFP(rippleChannels{ii}.Sharpwave_Channel,'restrict',timestamps);
                  
                  lfp1 = bz_GetLFP(rippleChannels.Ripple_Channel,'restrict',timestamps);
                  lfp2 = bz_GetLFP(rippleChannels.Sharpwave_Channel,'restrict',timestamps);
                  % Compute Coherence between one channel of each shank
                  %coherencePerShank.(foldername) = bz_coherencePerShank('timestampsSubSession',timestamps,'foldername',foldername,'saveMat',false);
                  % Compute Coherence between two channels
%                   coherogram.(foldername) = bz_MTCoherogram(lfp1,lfp2,'foldername',foldername,'range',[0 200], 'window',1,'saveMat',false);
                  coherogram.(foldername) = bz_MTCoherogram_v2(lfp1,lfp2,'foldername',foldername,'range',[0 200], 'window',1,'saveMat',false);
                  
                  % Compute Cross-Frequency Coupling
                  lfp = bz_GetLFP('all','restrict',timestamps);
%                   CFCPhaseAmp_lg.(foldername) = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[sgFreq(1):1:sgFreq(2)],'phaseCh',rippleChannels{ii}.Ripple_Channel,'ampCh',rippleChannels{ii}.Ripple_Channel,'foldername',foldername,'saveMat',false);
%                   CFCPhaseAmp_hg.(foldername) = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[hgFreq(1):1:hgFreq(2)],'phaseCh',rippleChannels{ii}.Ripple_Channel,'ampCh',rippleChannels{ii}.Ripple_Channel,'foldername',foldername,'saveMat',false);
%                   CFCPhaseAmp_hfo.(foldername) = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[hfoFreq(1):1:hfoFreq(2)],'phaseCh',rippleChannels{ii}.Ripple_Channel,'ampCh',rippleChannels{ii}.Ripple_Channel,'foldername',foldername,'saveMat',false);
%                   CFCPhaseAmp.(foldername) = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[1:1:200],'phaseCh',rippleChannels{ii}.Ripple_Channel,'ampCh',rippleChannels{ii}.Ripple_Channel,'foldername',foldername,'saveMat',false);

                  CFCPhaseAmp_lg.(foldername) = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[sgFreq(1):1:sgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                  CFCPhaseAmp_hg.(foldername) = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[hgFreq(1):1:hgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                  CFCPhaseAmp_hfo.(foldername) = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[hfoFreq(1):1:hfoFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
%                   CFCPhaseAmp.(foldername) = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[1:1:200],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                  
                  
%                   PhaseAmpCoupling_lg.(foldername) = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)],'phaseCh',rippleChannels{ii}.Ripple_Channel,'ampCh',rippleChannels{ii}.Ripple_Channel,'foldername',foldername,'saveMat',false);
%                   PhaseAmpCoupling_hg.(foldername) = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [hgFreq(1) hgFreq(2)],'phaseCh',rippleChannels{ii}.Ripple_Channel,'ampCh',rippleChannels{ii}.Ripple_Channel,'foldername',foldername,'saveMat',false);
%                   PhaseAmpCoupling_hfo.(foldername) = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [hfoFreq(1) hfoFreq(2)],'phaseCh',rippleChannels{ii}.Ripple_Channel,'ampCh',rippleChannels{ii}.Ripple_Channel,'foldername',foldername,'saveMat',false);                  
%                   PhaseAmpCoupling.(foldername) = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [1 200],'phaseCh',rippleChannels{ii}.Ripple_Channel,'ampCh',rippleChannels{ii}.Ripple_Channel','foldername',foldername,'saveMat',false);

                  PhaseAmpCoupling_lg.(foldername) = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                  PhaseAmpCoupling_hg.(foldername) = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [hgFreq(1) hgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                  PhaseAmpCoupling_hfo.(foldername) = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [hfoFreq(1) hfoFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);                  
%                   PhaseAmpCoupling.(foldername) = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [1 200],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel','foldername',foldername,'saveMat',false);  
                  % Power-Power Correlation
                  %lfp = bz_GetLFP([rippleChannels{ii}.Ripple_Channel rippleChannels{ii}.Sharpwave_Channel],'restrict',timestamps);
                  %comodulogram.(foldername) = bz_Comodulogram(lfp,[],[],'foldername',foldername,'ch1',rippleChannels{ii}.Ripple_Channel, 'ch2',rippleChannels{ii}.Sharpwave_Channel,'saveMat',false);
                  % Cross-Spectral Matrix �How to plot?
                  %bz_MTCrossSpec(lfp);
                  
                  % GMI Index
                  % Maybe we can compute this index over a wide range of
                  % channels (as an input)
                  % All Channels
%                   GMI_lg.(foldername) = bz_GMI_v4([thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)],'timestamps',timestamps,'foldername',foldername,'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels,'saveMat',false);
                  % Two Channels ( rippleChannel and SharpwaveChannel)
                  GMI_lg.(foldername) = bz_GMI_v4([thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)],'timestamps',timestamps,'foldername',foldername,'phaseCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'ampCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'saveMat',false);
                  % One Channel ( Ripple Channel)
%                   GMI_lg.(foldername) = bz_GMI_v4([thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)],'timestamps',timestamps,'foldername',foldername,'phaseCh',[rippleChannels.Ripple_Channel],'ampCh',[rippleChannels.Ripple_Channel],'saveMat',false);
                  % One Channel for phase ( RippleChannel) and two channels
                  % for amp ( rippleChannel and Sharpwave Channel)
%                   GMI_lg.(foldername) = bz_GMI_v4([thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)],'timestamps',timestamps,'foldername',foldername,'phaseCh',[rippleChannels.Ripple_Channel rippleChannels.Ripple_Channel],'ampCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'saveMat',false);
                  GMI_hg.(foldername) = bz_GMI_v4([thetaFreq(1) thetaFreq(2)], [hgFreq(1) hgFreq(2)],'timestamps',timestamps,'foldername',foldername,'phaseCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'ampCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'saveMat',false);
%                   GMI_hg.(foldername) = bz_GMI_v4([thetaFreq(1) thetaFreq(2)], [hgFreq(1) hgFreq(2)],'timestamps',timestamps,'foldername',foldername,'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels,'saveMat',false);
                  GMI_hfo.(foldername) = bz_GMI_v4([thetaFreq(1) thetaFreq(2)], [hfoFreq(1) hfoFreq(2)],'timestamps',timestamps,'foldername',foldername,'phaseCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'ampCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'saveMat',false);
%                   GMI_hfo.(foldername) = bz_GMI_v4([thetaFreq(1) thetaFreq(2)], [hfoFreq(1) hfoFreq(2)],'timestamps',timestamps,'foldername',foldername,'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels,'saveMat',false);

                  % GMI Index Tort
                  % Maybe we can compute this index over a wide range of
                  % channels (as an input)
                  lfp = bz_GetLFP('all','restrict',timestamps);
                  MI_lg.(foldername) = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)], 'phaseCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'ampCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'foldername',foldername,'saveMat',false);
%                   MI_lg.(foldername) = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)], 'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels,'foldername',foldername,'saveMat',false);

                  MI_hg.(foldername) = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)],[hgFreq(1) hgFreq(2)], 'phaseCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'ampCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'foldername',foldername,'saveMat',false);
%                   MI_hg.(foldername) = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)],[hgFreq(1) hgFreq(2)], 'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels,'foldername',foldername,'saveMat',false);
                  
                  MI_hfo.(foldername) = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)], [hfoFreq(1) hfoFreq(2)], 'phaseCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'ampCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'foldername',foldername,'saveMat',false);
%                   MI_hfo.(foldername) = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)], [hfoFreq(1) hfoFreq(2)], 'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels,'foldername',foldername,'saveMat',false);
              end
              
              % Saving .mat files
              % Coherence per Shank
              if exist('coherencePerShank','var')
                save([basepath filesep sessionInfo.FileName,'.Coherence_Shanks.SubSession.lfp.mat'],'coherencePerShank');
              end
              % Coherogram
              save([basepath filesep sessionInfo.FileName,'.Coherogram.SubSession.lfp.mat'],'coherogram','-v7.3');
              % CFCPhaseAmp lg
              save([basepath filesep sessionInfo.FileName,'.CFCPhaseAmp_',num2str(sgFreq(1)),'_',num2str(sgFreq(2)),'.SubSession.lfp.mat'],'CFCPhaseAmp_lg');
              % CFCPhaseAmp hg
              save([basepath filesep sessionInfo.FileName,'.CFCPhaseAmp_',num2str(hgFreq(1)),'_',num2str(hgFreq(2)),'.SubSession.lfp.mat'],'CFCPhaseAmp_hg');
              % CFCPhaseAmp hfo
              save([basepath filesep sessionInfo.FileName,'.CFCPhaseAmp_',num2str(hfoFreq(1)),'_',num2str(hfoFreq(2)),'.SubSession.lfp.mat'],'CFCPhaseAmp_hfo');
              % CFCPhaseAmp 
              if exist('CFCPhaseAmp','var')
                save([basepath filesep sessionInfo.FileName,'.CFCPhaseAmp.SubSession.lfp.mat'],'CFCPhaseAmp');
              end
              % PhaseAmpCouplingByAmp lg
              save([basepath filesep sessionInfo.FileName,'.PhaseAmpCouplingByAmp_',num2str(sgFreq(1)),'_',num2str(sgFreq(2)),'.SubSession.lfp.mat'],'PhaseAmpCoupling_lg');
              % PhaseAmpCouplingByAmp hg
              save([basepath filesep sessionInfo.FileName,'.PhaseAmpCouplingByAmp_',num2str(hgFreq(1)),'_',num2str(hgFreq(2)),'.SubSession.lfp.mat'],'PhaseAmpCoupling_hg');
              % PhaseAmpCouplingByAmp hfo
              save([basepath filesep sessionInfo.FileName,'.PhaseAmpCouplingByAmp_',num2str(hfoFreq(1)),'_',num2str(hfoFreq(2)),'.SubSession.lfp.mat'],'PhaseAmpCoupling_hfo');
              % PhaseAmpCouplingByAmp
              if exist('PhaseAmpCoupling','var')
                save([basepath filesep sessionInfo.FileName,'.PhaseAmpCouplingByAmp.SubSession.lfp.mat'],'PhaseAmpCoupling');
              end
              % Comodulogram
              %save([basepath filesep sessionInfo.FileName,'Comodulogram'],'comodulogram');
              % GMI lg
              save([basepath filesep sessionInfo.FileName,'.GMI_',num2str(sgFreq(1)),'_',num2str(sgFreq(2)),'.SubSession.lfp.mat'],'GMI_lg');
              % GMI hg
              save([basepath filesep sessionInfo.FileName,'.GMI_',num2str(hgFreq(1)),'_',num2str(hgFreq(2)),'.SubSession.lfp.mat'],'GMI_hg');
              % GMI hfo
              save([basepath filesep sessionInfo.FileName,'.GMI_',num2str(hfoFreq(1)),'_',num2str(hfoFreq(2)),'.SubSession.lfp.mat'],'GMI_hfo');
              % GMI
              if exist('GMI','var')
                save([basepath filesep sessionInfo.FileName,'GMI.SubSession.lfp.mat'],'GMI');
              end
              % GMI Tort lg
              save([basepath filesep sessionInfo.FileName,'.GMI_Tort_',num2str(sgFreq(1)),'_',num2str(sgFreq(2)),'.SubSession.lfp.mat'],'MI_lg');
              % GMI Tort hg
              save([basepath filesep sessionInfo.FileName,'.GMI_Tort_',num2str(hgFreq(1)),'_',num2str(hgFreq(2)),'.SubSession.lfp.mat'],'MI_hg');
              % GMI Tort hfo
              save([basepath filesep sessionInfo.FileName,'.GMI_Tort_',num2str(hfoFreq(1)),'_',num2str(hfoFreq(2)),'.SubSession.lfp.mat'],'MI_hfo');
              % GMI Tort
              if exist('MI','var')
                save([basepath filesep sessionInfo.FileName,'GMI_Tort.SubSession.lfp.mat'],'MI');
              end
          end
          
      catch
          disp('It is not possible to perform lfp_analysis SubSession...')         
      end      
  else
      try
          disp('Computing lfp_analysis for Whole Session...')
          session = loadSession;
          sessionInfo = bz_getSessionInfo(pwd,'noPrompts',true);
          shanks = sessionInfo.AnatGrps;
          shanks(excludeShanks) = [];
          
         % Pick Channel of Ripples and SW
          if ~isempty(dir([basepath filesep '*.channelinfo.ripples.mat']))
             disp('Channelinfo.ripples.mat file detected ! Loading file.')
             file = dir([basepath filesep '*channelinfo.ripples.mat'])
             load(file.name);
          end
          lfp1 = bz_GetLFP(rippleChannels.Ripple_Channel);
          lfp2 = bz_GetLFP(rippleChannels.Sharpwave_Channel); 
          % Compute Coherence between one channel of each shank
          coherencePerShank = bz_coherencePerShank();
          % Compute Coherence between two channels
          bz_MTCoherogram(lfp1,lfp2,'range',[0 200], 'window',1);
          % Compute Cross-Frequency Coupling
          lfp = bz_GetLFP('all');
          CFCPhaseAmp_lg = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[sgFreq(1):1:sgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);
          CFCPhaseAmp_hg = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[hgFreq(1):1:hgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);
          CFCPhaseAmp_hfo = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[hfoFreq(1):1:hfoFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);
          CFCPhaseAmp = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[1:1:200],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);

          PhaseAmpCoupling_lg = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);
          PhaseAmpCoupling_hg = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [hgFreq(1) hgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);
          PhaseAmpCoupling_hfo = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [hfoFreq(1) hfoFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);                  
          PhaseAmpCoupling = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [1 200],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);                  
          % Power-Power Correlation
          lfp = bz_GetLFP([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel]);
          comodulogram = bz_Comodulogram(lfp,[],[],'ch1',rippleChannels.Ripple_Channel, 'ch2',rippleChannels.Sharpwave_Channel);
                  

          % GMI Index
          % Maybe we can compute this index over a wide range of
          % channels (as an input)
          lfp = bz_GetLFP('all');
          GMI_lg = bz_GMI(lfp,[thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)],'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);
%                   bz_GMI(lfp,[4 12], [30 80],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername);
          GMI_hg = bz_GMI(lfp,[thetaFreq(1) thetaFreq(2)], [hgFreq(1) hgFreq(2)],'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);
%                   bz_GMI(lfp,[4 12], [80 150],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername);
          GMI_hfo = bz_GMI(lfp,[thetaFreq(1) thetaFreq(2)],[hfoFreq(1) hfoFreq(2)],'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);
          GMI = bz_GMI(lfp,[thetaFreq(1) thetaFreq(2)], [1 200],'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);
%                   bz_GMI(lfp,[4 12], [1 200],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername);

          % GMI Index Tort
          % Maybe we can compute this index over a wide range of
          % channels (as an input)
          lfp = bz_GetLFP('all');
          MI_lg = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)], 'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);
          MI_hg = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)],[hgFreq(1) hgFreq(2)], 'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);
          MI_hfo = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)],[hfoFreq(1) hfoFreq(2)], 'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);
          MI = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)], [1 200], 'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);    
      catch
          disp('It is not possible to run lfp_analysis for Whole Session..')
      end  
  end    
end

%% NEW PROTOCOL
cd(basepath)
% pass = 15*60;
if any(ismember(listOfAnalysis,'newProtocol'))
    mkdir('newProtocolFigures')
    try
        disp('Computing Analysis for NewProtocol...');
        session = loadSession;
        sessionInfo = bz_getSessionInfo(pwd,'noPrompts',true);
        shanks = sessionInfo.AnatGrps;
        shanks(excludeShanks) = [];
        if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
            load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
            count = 1;
             if ~isempty(dir([basepath filesep '*.channelinfo.ripples.mat']))
                 disp('Channelinfo.ripples.mat file detected ! Loading file.')
                 file = dir([basepath filesep '*channelinfo.ripples.mat'])
                 load(file.name);
             end
             count1 = 1;
             for ii = 1:size(MergePoints.foldernames,2)
                 timestamps = MergePoints.timestamps(ii,:);
                 foldername = MergePoints.foldernames{ii};
                 times = zeros(round((timestamps(end)-timestamps(1))/pass),2);
                 for i = 1:size(times,1)
                     if i == 1
                         times(1,:) = [timestamps(1) timestamps(1)+pass];
                     else
                         times(i,:) = [timestamps(1)+pass*(i-1) timestamps(1)+pass*i];
                     end
                 end
                 if times(end,end) > timestamps(end)
                     times(end,end) = timestamps(end);
                 end
                 
                 for jj = 1:size(times,1)
                     lfp1 = bz_GetLFP(rippleChannels.Ripple_Channel,'restrict',times(jj,:));
                     lfp2 = bz_GetLFP(rippleChannels.Sharpwave_Channel,'restrict',times(jj,:));
                     % Compute Coherence between onw channel of each shank
%                      coherencePerShank = bz_coherenceShank();
                     % Compute Coherence between two channels
                     coherogram.(foldername){jj} = bz_MTCoherogram_NP(lfp1,lfp2,'foldername',foldername,'range',[0 200],'window',1,'saveMat',false);
                     % Compute Cross-frequency Coupling
                     lfp = bz_GetLFP('all','restrict',times(jj,:));
                     CFCPhaseAmp_lg.(foldername){jj} = bz_CFCPhaseAmp_NP(lfp,[thetaFreq(1):1:thetaFreq(2)],[sgFreq(1):1:sgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                     CFCPhaseAmp_hg.(foldername){jj} = bz_CFCPhaseAmp_NP(lfp,[thetaFreq(1):1:thetaFreq(2)],[hgFreq(1):1:hgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                     CFCPhaseAmp_hfo.(foldername){jj} = bz_CFCPhaseAmp_NP(lfp,[thetaFreq(1):1:thetaFreq(2)],[hfoFreq(1):1:hfoFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                     
                     PhaseAmpCoupling_lg.(foldername){jj} = bz_PhaseAmpCouplingByAmp_NP(lfp,[thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                     PhaseAmpCoupling_hg.(foldername){jj} = bz_PhaseAmpCouplingByAmp_NP(lfp,[thetaFreq(1) thetaFreq(2)], [hgFreq(1) hgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                     PhaseAmpCoupling_hfo.(foldername){jj} = bz_PhaseAmpCouplingByAmp_NP(lfp,[thetaFreq(1) thetaFreq(2)], [hfoFreq(1) hfoFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);   

                     GMI_lg.(foldername){jj} = bz_GMI_NP([thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)],'timestamps',times(jj,:),'foldername',foldername,'phaseCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'ampCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'saveMat',false);
                     GMI_hg.(foldername){jj} = bz_GMI_NP([thetaFreq(1) thetaFreq(2)], [hgFreq(1) hgFreq(2)],'timestamps',times(jj,:),'foldername',foldername,'phaseCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'ampCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'saveMat',false);
                     GMI_hfo.(foldername){jj} = bz_GMI_NP([thetaFreq(1) thetaFreq(2)], [hfoFreq(1) hfoFreq(2)],'timestamps',times(jj,:),'foldername',foldername,'phaseCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'ampCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'saveMat',false);

                     lfp = bz_GetLFP('all','restrict',times(jj,:));
                     MI_lg.(foldername){jj} = bz_GMITort_NP(lfp,[thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)], 'phaseCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'ampCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'foldername',foldername,'saveMat',false);
                     MI_hg.(foldername){jj} = bz_GMITort_NP(lfp,[thetaFreq(1) thetaFreq(2)],[hgFreq(1) hgFreq(2)], 'phaseCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'ampCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'foldername',foldername,'saveMat',false);
                     MI_hfo.(foldername){jj} = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)], [hfoFreq(1) hfoFreq(2)], 'phaseCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'ampCh',[rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'foldername',foldername,'saveMat',false);
                     
                 end
             end
             
             % Saving.mat files
             if exist('coherencePerShank','var')
                 save([basepath filesep sessionInfo.FileName,'.Coherence_Shanks.SubSession.lfp.mat'],'coherencePerShank');
             end
             % Coherogram
             save([basepath filesep sessionInfo.FileName,'.Coherogram.SubSession.lfp.mat'],'coherogram','-v7.3');
             % CFCPhaseAmp lg
             save([basepath filesep sessionInfo.FileName,'.CFCPhaseAmp_',num2str(sgFreq(1)),'_',num2str(sgFreq(2)),'.SubSession.lfp.mat'],'CFCPhaseAmp_lg');
             % CFCPhaseAmp hg
             save([basepath filesep sessionInfo.FileName,'.CFCPhaseAmp_',num2str(hgFreq(1)),'_',num2str(hgFreq(2)),'.SubSession.lfp.mat'],'CFCPhaseAmp_hg');
              % CFCPhaseAmp hfo
             save([basepath filesep sessionInfo.FileName,'.CFCPhaseAmp_',num2str(hfoFreq(1)),'_',num2str(hfoFreq(2)),'.SubSession.lfp.mat'],'CFCPhaseAmp_hfo');
             % CFCPhaseAmp 
             if exist('CFCPhaseAmp','var')
               save([basepath filesep sessionInfo.FileName,'.CFCPhaseAmp.SubSession.lfp.mat'],'CFCPhaseAmp');
             end
             % PhaseAmpCouplingByAmp lg
             save([basepath filesep sessionInfo.FileName,'.PhaseAmpCouplingByAmp_',num2str(sgFreq(1)),'_',num2str(sgFreq(2)),'.SubSession.lfp.mat'],'PhaseAmpCoupling_lg');
             % PhaseAmpCouplingByAmp hg
             save([basepath filesep sessionInfo.FileName,'.PhaseAmpCouplingByAmp_',num2str(hgFreq(1)),'_',num2str(hgFreq(2)),'.SubSession.lfp.mat'],'PhaseAmpCoupling_hg');
             % PhaseAmpCouplingByAmp hfo
             save([basepath filesep sessionInfo.FileName,'.PhaseAmpCouplingByAmp_',num2str(hfoFreq(1)),'_',num2str(hfoFreq(2)),'.SubSession.lfp.mat'],'PhaseAmpCoupling_hfo');
             % PhaseAmpCouplingByAmp
             if exist('PhaseAmpCoupling','var')
               save([basepath filesep sessionInfo.FileName,'.PhaseAmpCouplingByAmp.SubSession.lfp.mat'],'PhaseAmpCoupling');
             end
             % Comodulogram
             %save([basepath filesep sessionInfo.FileName,'Comodulogram'],'comodulogram');
             % GMI lg
             save([basepath filesep sessionInfo.FileName,'.GMI_',num2str(sgFreq(1)),'_',num2str(sgFreq(2)),'.SubSession.lfp.mat'],'GMI_lg');
             % GMI hg
             save([basepath filesep sessionInfo.FileName,'.GMI_',num2str(hgFreq(1)),'_',num2str(hgFreq(2)),'.SubSession.lfp.mat'],'GMI_hg');
             % GMI hfo
             save([basepath filesep sessionInfo.FileName,'.GMI_',num2str(hfoFreq(1)),'_',num2str(hfoFreq(2)),'.SubSession.lfp.mat'],'GMI_hfo');
             % GMI
             if exist('GMI','var')
               save([basepath filesep sessionInfo.FileName,'GMI.SubSession.lfp.mat'],'GMI');
             end
             % GMI Tort lg
             save([basepath filesep sessionInfo.FileName,'.GMI_Tort_',num2str(sgFreq(1)),'_',num2str(sgFreq(2)),'.SubSession.lfp.mat'],'MI_lg');
             % GMI Tort hg
             save([basepath filesep sessionInfo.FileName,'.GMI_Tort_',num2str(hgFreq(1)),'_',num2str(hgFreq(2)),'.SubSession.lfp.mat'],'MI_hg');
             % GMI Tort hfo
             save([basepath filesep sessionInfo.FileName,'.GMI_Tort_',num2str(hfoFreq(1)),'_',num2str(hfoFreq(2)),'.SubSession.lfp.mat'],'MI_hfo');
             % GMI Tort
             if exist('MI','var')
               save([basepath filesep sessionInfo.FileName,'GMI_Tort.SubSession.lfp.mat'],'MI');
             end
        end    
    catch
        disp('It is not possible to run Analysis for new Protocol...');
    end
end
%% 10 - PLOT PLACE FIELDS
if any(ismember(listOfAnalysis,'plotPlaceFields'))
    try
        disp('Plotting Place Fields')
        sessionInfo = bz_getSessionInfo(basepath);
        spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat,'forceReload',false,'showWaveforms',showWaveforms);
        tracking = getSessionTracking();
        behaviour = getSessionBehaviour_v2();
        spikeTrain = bz_SpikeTrain(spikes,'analyzeSubSessions',analyzeSubSessions,'showFigure',true);
        meanFr = bz_meanFr(spikes,'isLinearTrack',any(ismember(listOfAnalysis,'plotLinearTrack')));
        firingMaps = bz_firingMapAvg(behaviour,spikes,'saveMat',true,'speedFilter',true,'periodicAnalysis',false,'spikeShuffling',false); 
        if ~isempty([basepath filesep sessionInfo.FileName, '.cell_metrics.cellinfo.mat'])
            disp('Loading Cell Metrics...')
            file = dir([basepath filesep sessionInfo.FileName,'.cell_metrics.cellinfo.mat'])
            load(file.name)
        end
%         plot_placeFields('firingMaps',firingMaps,'spikes',spikes,'tracking',tracking,'cell_metrics',cell_metrics,'spikeTrain',spikeTrain);
        bz_plotPlaceFields('firingMaps',firingMaps,'spikes',spikes','tracking',tracking,'cell_metrics',cell_metrics,'spikeTrain',spikeTrain);        
        if any(ismember(listOfAnalysis,'plotLinearTrack'))
            bz_plotPlaceFields1D('firingMaps',firingMaps,'spikes',spikes,'tracking',tracking,'cell_metrics',cell_metrics,'spikeTrain',spikeTrain);
        end
    catch
        disp('It is not possible to run plot Place Fields')
    end
end

%% 9 - CREATION EXCEL FILE 

if any(ismember(listOfAnalysis,'excel'))
    if analyzeSubSessions
        try
            disp('Creating Excel file for SubSessions ...')
            bz_createExcelSubSessions(listOfAnalysis,'basepath',basepath,'diffLFPs',diffLFPs,'thetaFreq',thetaFreq,'sgFreq',sgFreq,'hgFreq',hgFreq,'hfoFreq',hfoFreq,'pathExcel',pathExcel,'nameExcel',nameExcel);
        catch
            warning('It has not been possible to create Excel file for SubSession...')
        end
    else
        try
            disp('Creating Excel file...')
            excel = bz_createExcel(listOfAnalysis,'basepath',basepath,'diffLFPs',diffLFPs,'thetaFreq',thetaFreq,'sgFreq',sgFreq,'hgFreq',hgFreq,'pathExcel',nameExcel);
        catch
            warning('It has not been possible to create Excel file...')
        end
    end
end
