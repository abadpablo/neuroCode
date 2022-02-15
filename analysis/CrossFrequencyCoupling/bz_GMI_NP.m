function [GMI] = bz_GMI_NP(phaserange,amprange,varargin)
% [GMI] = bz_GMI(lfp,varargin)
%
%
% This function computes the Gamma Modulation Index between a phase
% filtered signal and an amplitude filtered signal
%
%   INPUT
%
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                                                   lfp.channels
%                   lfp will be get inside this function
%
%   <options>       optional list of property-value pairs (see table below)
%
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     phaseCh      channel to compute phase. If empty take first channel
%     ampChans     channels to compute amplitude. If empty take first channel
%     method       ['wavelet'(default)|'hilbert']. Method to extract power
%        of specific bands, using wavelet (default) or hilbert
%
%     makePlot      default true
%     filterType    default 'fir1'. Method of filtering for the phase
%        bands, in case of method = 'hilbert' it also defines the filter of
%        bands in the amplitude range.
%     filterOrder   default 4. Order of the filter used for the phase
%        bands, in case of method = 'hilbert' it also defines the filter order of
%        bands in the amplitude range.
%     numBins       default 50. Number of bins to calculate the
%        Amplitude/Phase distribution.
%     perm_test     default false. To whether calculate or not a surrogate
%     test for CFC values.
%     alpha         default 0.05. Alpha for the surrogate test
%     num_inter     default 200.  Number of permutations
%     units         default MI. alternative: zscore. Whether the units are 
%     in MI or zscore after the permutation test
%    =========================================================================
%
%
%   OUTPUT
%
%
%
%
% Dependencies
%   bz_Filter
%   bz_WaveSpec
%
%   Pablo Abad PÃ©rez 2021
% Corrected for different epochs when tracking is present


%% Default params
p = inputParser();
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'sr',1250,@isnumeric);
addParameter(p,'phaseCh',[],@isnumeric);
addParameter(p,'ampCh',[],@isnumeric);
addParameter(p,'filterType','fir1',@isstring);
addParameter(p,'filterOrder',4,@isnumeric);
addParameter(p,'numBins',50,@isnumeric);
addParameter(p,'method','wavelet',@ischar);

addParameter(p,'foldername',[],@isstr);
addParameter(p,'showFig',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveFig',true,@islogical);
addParameter(p,'exist_file',false,@islogical);
addParameter(p,'speed_filter',false,@islogical); % to filter the signal for speed (HomeCage recordings do not have tracking information)Â¿?Â¿?Â¿?Â¿
addParameter(p,'noise_filter',false,@islogical);
addParameter(p,'epoch_analysis',false,@islogical); % Equal to previous parameter. If including in the analysis the filtering of the signals to remove TD ratio noise
addParameter(p,'debugMode',false,@islogical);
addParameter(p,'thres',3,@isnumeric);
addParameter(p,'timestamps',[],@isnumeric);

% params for Chronux functions

parse(p,varargin{:})
basepath = p.Results.basepath;
sr = p.Results.sr;
phaseCh = p.Results.phaseCh;
ampCh = p.Results.ampCh;
filterType = p.Results.filterType;
filterOrder = p.Results.filterOrder;
numBins = p.Results.numBins;
method = p.Results.method;

foldername = p.Results.foldername;
showFig = p.Results.showFig;
saveMat = p.Results.saveMat;
saveFig = p.Results.saveFig;
exist_file = p.Results.exist_file;
speed_filter = p.Results.speed_filter;
noise_filter = p.Results.noise_filter;
epoch_analysis = p.Results.epoch_analysis;
debugMode = p.Results.debugMode;
thres = p.Results.thres;
timestamps = p.Results.timestamps;

%% Load sessionInfo
sessionInfo = bz_getSessionInfo(basepath);
shanks = sessionInfo.AnatGrps;
%% Load tracking
tracking = getSessionTracking();
%% Load MergePoints
if ~isempty(dir([basepath filesep sessionInfo.FileName '.MergePoints.events.mat']))
    file = dir([basepath filesep sessionInfo.FileName '.MergePoints.events.mat']);
    load(file.name)
end
%% Input Channels
lfp_phase = bz_GetLFP(phaseCh,'restrict',timestamps);
lfp_amp = bz_GetLFP(ampCh,'restrict',timestamps);

if ismember(foldername,tracking.folders)
    for i=1:length(tracking.folders)
        if isequal(foldername,tracking.folders{i})
            trackingEvent = i;
        end
    end
    for ii=1:length(lfp_phase.channels)
        % Speed Limits
        ts = InIntervals(tracking.timestamps,timestamps);
        v = tracking.position.speed(ts);
        timestamp = tracking.timestamps(ts);
%         v = tracking.position.speed(tracking.events.subSessionsMask == trackingEvent);
%         timestamp = tracking.timestamps(tracking.events.subSessionsMask == trackingEvent);
%         ts = InIntervals(tracking.timestamps,timestamps);
%         v = v(ts);
        Vi = 0.05;
        Vs = 0.5;
        indS1 = find(v<Vi);
        indS2 = find(v>Vs);
        mat = v;
        mat = mat';
        mat2 = mat;
        mat2(indS1) = 0;
        mat2(indS2) = 0;
        mat2(isnan(mat2)) = 0;
    
        % Finding epochs
        for i=1:length(MergePoints.foldernames)
            if strcmpi(MergePoints.foldernames{i},foldername)
                file_number = i;
            end
        end
        vvv = bwlabel(squeeze(mat2));
        MINPFSIZE = tracking.samplingRate(trackingEvent)*2;
        % We need to make timestamps
        tr_timestamps = tracking.timestamps(ts);
%         tr_timestamps = tracking.timestamps(tracking.events.subSessionsMask == trackingEvent);
%         tr_timestamps = tr_timestamps(ts);
        t = lfp_phase.timestamps;
        [n2,bin2] = histc(tr_timestamps,t);
    
        % Epoch obtention
        for jj = 1:max(vvv(:))
            ttt = sum(vvv(:) == jj);
            [rb,cb] = find(bwlabel(vvv) == jj);
            if ttt > MINPFSIZE
                if bin2(cb(1)) > 0 && bin2(cb(end)) > 0 
                    time_tracking = tr_timestamps(cb(1):cb(end));
                    lfp_phase_epoch = lfp_phase.data(bin2(cb(1)):bin2(cb(end)),ii);
                    lfp_amp_epoch = lfp_amp.data(bin2(cb(1)):bin2(cb(end)),ii);
                    time_epoch = lfp_phase.timestamps(bin2(cb(1)):bin2(cb(end)));
                    [Eband_phase{ii},RatTD_phase{ii},pphase{ii},fphase{ii}] = bz_MakPower2Chronux(lfp_phase_epoch);
                    [Eband_amp{ii},RatTD_amp{ii},pamp{ii},famp{ii}] = bz_MakPower2Chronux(lfp_amp_epoch);

                    if RatTD_phase{ii} > thres && RatTD_amp{ii} > thres
                        [DegMod,DegreeS,FrecG,GamPower,PowLoc3,ThetFr,GammaFr,AmpT,AmpG] = bz_HilbRelated(lfp_phase_epoch',lfp_amp_epoch',time_epoch','amprange',amprange);
                        DegMod_epoch{ii}{jj} = DegMod;
                        DegreeS_epoch{ii}{jj} = DegreeS;
                        FrecG_epoch{ii}{jj} = FrecG;
                        GamPower_epoch{ii}{jj} = GamPower;
                        PowLoc3_epoch{ii}{jj} = PowLoc3;
                        ThetFr_epoch{ii}{jj} = ThetFr;
                        GammaFr_epoch{ii}{jj} = GammaFr;
                        AmpT_epoch{ii}{jj} = AmpT;
                        AmpG_epoch{ii}{jj} = AmpG;

                    end
                end
            end
        end
    
        if exist('DegMod_epoch','var')
            % Organizing data
            DegMod_M{ii} = cell2mat(DegMod_epoch{ii});
            DegreeS_M{ii} = cell2mat(DegreeS_epoch{ii});
            FrecG_M{ii} = cell2mat(FrecG_epoch{ii}); %Gamma Frequencies
            GamPower_M{ii} = cell2mat(GamPower_epoch{ii});
            PowLoc3_M{ii} = cell2mat(PowLoc3_epoch{ii});
            ThetFr_M{ii} =mean(cell2mat(ThetFr_epoch{ii}));
            GammaFr_M{ii} =mean(cell2mat(GammaFr_epoch{ii}));
            AmpT_M{ii} = mean(cell2mat(AmpT_epoch{ii}));
            AmpG_M{ii} = mean(cell2mat(AmpG_epoch{ii}));
                
    
            % Hilbert transform plotting and modulation index calculation
            intSize=18;
            [GMI{ii}, GFMI{ii}, MeanFGa{ii}, MeanPowGa{ii}, MeanFStd{ii}, MeanPowStd{ii}, MaxF{ii}, MaxFLoc{ii}, MaxP{ii}, MaxPLoc{ii}, VarRag2{ii}, bindeRad{ii}, MaxdegG{ii}] = bz_HistHilbert(GamPower_M{ii},FrecG_M{ii},GammaFr_epoch{ii},DegreeS_M{ii},PowLoc3_M{ii},intSize);
            NoA=(find(isnan(VarRag2{ii})));
            VarRag2{ii}(NoA)=0;% eliminate non numeric values result of bining
            % Circular Statistics
            % Preparing data for Circular statistic

            binde{ii}=rad2deg(bindeRad{ii} )+180;
            edges=(-pi:0.3307:pi);
            data2{ii}=deg2rad(cell2mat(PowLoc3_epoch{ii})-180);
            hist{ii}=histc(data2{ii},edges); 
            hist{ii}(end)=[];
            [pval{ii},z{ii}]=circ_rtest((edges(1:end-1)+edges(2:end))/2,hist{ii});  %% 
            [mu{ii} ul{ii} ll{ii}] = circ_mean((edges(1:end-1)+edges(2:end))/2,hist{ii},2);
            mvl{ii} = circ_r((edges(1:end-1)+edges(2:end))/2,hist{ii},[],2);% mean vector length
            mu_out{ii} = circ_rad2ang(mu{ii})+180;
            ul_out{ii} = circ_rad2ang(ul{ii})+180;
            ll_out{ii} = circ_rad2ang(ll{ii})+180;
        end
    end
    
    if exist('DegMod_epoch','var')
        if showFig
            binde = binde{1};
            meanpowga = cell2mat(MeanPowGa);
            stdpowga = cell2mat(MeanPowStd);
            minpowga = min(meanpowga);
            maxpowga = max(meanpowga);
            minpowstd = min(stdpowga);
            maxpowstd = max(stdpowga);
            if length(lfp_phase.channels) == length(sessionInfo.channels) && length(lfp_amp.channels) == length(sessionInfo.channels)
                figure,
                set(gcf,'Position',get(0,'ScreenSize'))
                for jj=1:size(shanks,2)
                    count = jj;
                    for ii=1:length(shanks(jj).Channels)
                        subplot(sessionInfo.nChannels/sessionInfo.nElecGps,sessionInfo.nChannels/sessionInfo.nElecGps,count)
                        errorbar([binde binde+binde(end)],[MeanPowGa{shanks(jj).Channels(ii)+1} MeanPowGa{shanks(jj).Channels(ii)+1}],[MeanPowStd{shanks(jj).Channels(ii)+1} MeanPowStd{shanks(jj).Channels(ii)+1}]);
                        axis([0 720 minpowga-maxpowstd*2 maxpowga+maxpowstd*2]) 
                        if jj == 1 && ii ~= length(shanks(jj).Channels)
                            ylabel('Gamma Power');
                            set(gca,'XTick',[]);
                        elseif jj == 1 && ii == length(shanks(jj).Channels)
                            xlabel('Theta Degrees');
                            ylabel('GammaPower');
                        elseif jj~= 1 && ii == length(shanks(jj).Channels)
                            xlabel('Theta Degrees')
                            set(gca,'YTick',[])
                        else
                            set(gca,'XTick',[])
                            set(gca,'YTIck',[])
                        end
                        title(['GMI:', num2str(GMI{shanks(jj).Channels(ii)+1}),' Ch:', num2str(shanks(jj).Channels(ii)+1)])
                        count = count+size(shanks,2);
                    end
                end
            elseif length(lfp_phase.channels) == 1 && length(lfp_amp.channels) == 1
                figure,
                set(gcf,'Position',get(0,'ScreenSize'))
                errorbar([binde binde+binde(end)],[MeanPowGa{1} MeanPowGa{1}],[MeanPowStd{1} MeanPowStd{1}]);
                axis([0 720 minpowga-maxpowstd*2 maxpowga+maxpowstd*2]) 
                title(['GMI:', num2str(GMI{1}),' Ch Phase: ', num2str(phaseCh), ' Ch Amp :', num2str(ampCh)])

            elseif length(lfp_phase.channels) > 1 && length(lfp_amp.channels) > 1 && length(lfp_phase.channels) < length(sessionInfo.channels) && length(lfp_amp.channels) < length(sessionInfo.channels)
                figure,
                set(gcf,'Position',get(0,'ScreenSize'))
                for i=1:length(lfp_phase.channels)
                    subplot(length(lfp_phase.channels), 1, i)
                    errorbar([binde binde+binde(end)],[MeanPowGa{i} MeanPowGa{i}],[MeanPowStd{i} MeanPowStd{i}]);
                    axis([0 720 minpowga-maxpowstd*2 maxpowga+maxpowstd*2]) 
                    title(['GMI:', num2str(GMI{i}),' Ch Phase: ', num2str(phaseCh(i)), ' Ch Amp :', num2str(ampCh(i))])
                end
            end
            if saveFig
                if ~isempty(foldername)
                    saveas(gcf,['lfpAnalysisFigures\GMI.',foldername,'_',num2str(amprange(1)),num2str(amprange(end)),'_',num2str(timestamps(1)),'-',num2str(timestamps(end)),'.png']);
                else
                    saveas(gcf,['lfpAnalysisFigures\GMI_',num2str(amprange(1)),num2str(amprange(end)),'.png'])
                end
            end
        end
    end
    
else
    for i=1:length(lfp_phase.channels)
        [DegMod{i},DegreeS{i},FrecG{i},GamPower{i},PowLoc3{i},ThetFr{i},GammaFr{i},AmpT{i},AmpG{i}] = bz_HilbRelated(lfp_phase.data(:,i)',lfp_amp.data(:,i)',lfp_phase.timestamps','amprange',amprange);
        % Hilbert transform plotting and modulation index calculation
        intSize = 18;
        [GMI{i},GFMI{i},MeanFGa{i},MeanPowGa{i},MeanFStd{i},MeanPowStd{i},MaxF{i},MaxFLoc{i},MaxP{i},MaxPLoc{i},VarRag2{i},bindeRad{i},maxdegG{i}] = bz_HistHilbert(GamPower{i},FrecG{i},GammaFr{i},DegreeS{i},PowLoc3{i},intSize); 
        NoA=(find(isnan(VarRag2{i})));
        VarRag2{i}(NoA)=0;% eliminate non numeric values result of bining
        % Circular Statistics
        % Preparing data for Circular statistic
        binde{i}=rad2deg(bindeRad{i} )+180;
        edges=(-pi:0.3307:pi);
        data2{i}=deg2rad(PowLoc3{i}-180);
        hist{i}=histc(data2{i},edges); 
        hist{i}(end)=[];
        [pval{i},z{i}]=circ_rtest((edges(1:end-1)+edges(2:end))/2,hist{i});  %% 
        [mu{i} ul{i} ll{i}] = circ_mean((edges(1:end-1)+edges(2:end))/2,hist{i},2);
        mvl{i} = circ_r((edges(1:end-1)+edges(2:end))/2,hist{i},[],2);% mean vector length
        mu_out{i} = circ_rad2ang(mu{i})+180;
        ul_out{i} = circ_rad2ang(ul{i})+180;
        ll_out{i} = circ_rad2ang(ll{i})+180;
    end
    if showFig
        binde = binde{1};
        meanpowga = cell2mat(MeanPowGa);
        stdpowga = cell2mat(MeanPowStd);
        minpowga = min(meanpowga);
        maxpowga = max(meanpowga);
        minpowstd = min(stdpowga);
        maxpowstd = max(stdpowga);
        if length(lfp_phase.channels) == length(sessionInfo.channels) && length(lfp_amp.channels) == length(sessionInfo.channels)
            figure,
            set(gcf,'Position',get(0,'ScreenSize'))
            for jj=1:size(shanks,2)
                count = jj;
                for ii=1:length(shanks(jj).Channels)
                    subplot(sessionInfo.nChannels/sessionInfo.nElecGps,sessionInfo.nChannels/sessionInfo.nElecGps,count)
                    errorbar([binde binde+binde(end)],[MeanPowGa{shanks(jj).Channels(ii)+1} MeanPowGa{shanks(jj).Channels(ii)+1}],[MeanPowStd{shanks(jj).Channels(ii)+1} MeanPowStd{shanks(jj).Channels(ii)+1}]);
                    axis([0 720 minpowga-maxpowstd*2 maxpowga+maxpowstd*2]) 
                    if jj == 1 && ii ~= length(shanks(jj).Channels)
                        ylabel('Gamma Power');
                        set(gca,'XTick',[]);
                    elseif jj == 1 && ii == length(shanks(jj).Channels)
                        xlabel('Theta Degrees');
                        ylabel('GammaPower');
                    elseif jj~= 1 && ii == length(shanks(jj).Channels)
                        xlabel('Theta Degrees')
                        set(gca,'YTick',[])
                    else
                        set(gca,'XTick',[])
                        set(gca,'YTIck',[])
                    end
                    title(['GMI:', num2str(GMI{shanks(jj).Channels(ii)+1}),' Ch:', num2str(shanks(jj).Channels(ii)+1)])
                    count = count+size(shanks,2);
                end
            end
        elseif length(lfp_phase.channels) == 1 && length(lfp_amp.channels) == 1
            figure,
            set(gcf,'Position',get(0,'ScreenSize'))
            errorbar([binde binde+binde(end)],[MeanPowGa{1} MeanPowGa{1}],[MeanPowStd{1} MeanPowStd{1}]);
            axis([0 720 minpowga-maxpowstd*2 maxpowga+maxpowstd*2]) 
            title(['GMI:', num2str(GMI{1}),' Ch Phase: ', num2str(phaseCh), ' Ch Amp :', num2str(ampCh)])
            
        elseif length(lfp_phase.channels) > 1 && length(lfp_amp.channels) > 1 && length(lfp_phase.channels) < length(sessionInfo.channels) && length(lfp_amp.channels) < length(sessionInfo.channels)
            figure,
            set(gcf,'Position',get(0,'ScreenSize'))
            for i=1:length(lfp_phase.channels)
                subplot(length(lfp_phase.channels), 1, i)
                errorbar([binde binde+binde(end)],[MeanPowGa{i} MeanPowGa{i}],[MeanPowStd{i} MeanPowStd{i}]);
                axis([0 720 minpowga-maxpowstd*2 maxpowga+maxpowstd*2]) 
                title(['GMI:', num2str(GMI{i}),' Ch Phase: ', num2str(phaseCh(i)), ' Ch Amp :', num2str(ampCh(i))])
            end
        end
        if saveFig
            if ~isempty(foldername)
                saveas(gcf,['lfpAnalysisFigures\GMI.',foldername,'_',num2str(amprange(1)),num2str(amprange(end)),'_',num2str(timestamps(1)),'-',num2str(timestamps(end)),'.png']);
            else
                saveas(gcf,['lfpAnalysisFigures\GMI_',num2str(amprange(1)),num2str(amprange(end)),'.png'])
            end
        end
    end
end
 
%% save matlab structure
if exist('DegMod_epoch','var')
    gmi = GMI;
    GMI = [];
    GMI.GMI = gmi;
    GMI.GFMI = GFMI;
    if ~isempty(foldername)
        GMI.folder = foldername;
    end
    GMI.Channels.phaseCh = phaseCh;
    GMI.Channels.ampCh = ampCh;
    GMI.binde = binde;
    GMI.MeanPowGa = MeanPowGa;
    GMI.MeanPowStd = MeanPowStd;
    GMI.MeanFGa = MeanFGa;
    GMI.MeanFStd = MeanFStd;
else
    GMI.GMI = [];
    GMI.GFMI = [];
    if ~isempty(foldername)
        GMI.folder = foldername;
    end
    GMI.Channels.phaseCh = [];
    GMI.Channels.ampCh = [];
    GMI.binde = [];
    GMI.MeanPowGa = [];
    GMI.MeanPowStd = [];
    GMI.MeanFGa = [];
    GMI.MeanFStd = [];
end

if saveMat
    if ~isempty(foldername)
        try
                save([basepath filesep sessionInfo.FileName, '.', foldername, '.GMI_', num2str(amprange(1)), num2str(amprange(end)),'.lfp.mat'],'GMI');
        catch
            save([basepath filesep sessionInfo.FileName, '.', foldername, '.GMI_', num2str(amprange(1)), num2str(amprange(end)),'.lfp.mat'],'GMI','-v7.3')
        end
    else
        try
           save([basepath filesep sessionInfo.FileName,'.GMI_', num2str(amprange(1)), num2str(amprange(end)),'.lfp.mat'],'GMI') 
        catch
            save([basepath filesep sessionInfo.FileName,'.GMI_', num2str(amprange(1)), num2str(amprange(end)),'.lfp.mat'],'GMI','-v7.3')
        end
    end
end

end


