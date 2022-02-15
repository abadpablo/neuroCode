function [GMI] = bz_GMI_v2(lfp,phaserange,amprange,varargin)
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
%   Pablo Abad Pérez 2021

%% Default params
p = inputParser();
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'sr',1250,@isnumeric);
addParameter(p,'phaseCh',lfp.channels(1),@isnumeric);
addParameter(p,'ampCh',lfp.channels(1),@isnumeric);
addParameter(p,'filterType','fir1',@isstring);
addParameter(p,'filterOrder',4,@isnumeric);
addParameter(p,'numBins',50,@isnumeric);
addParameter(p,'method','wavelet',@ischar);

addParameter(p,'foldername',[],@isstr);
addParameter(p,'showFig',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveFig',true,@islogical);
addParameter(p,'exist_file',false,@islogical);
addParameter(p,'speed_filter',false,@islogical); % to filter the signal for speed (HomeCage recordings do not have tracking information)¿?¿?¿?¿
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
if isequal(length(phaseCh),length(ampCh)) && length(phaseCh) ~= 1
    disp('Same number of phaseChannels and amplitudeChannels. More than 1 channel')
elseif length(phaseCh) == 1 && ~isequal(length(phaseCh),length(ampCh))
    disp('One phaseChannel (ref) and more than 1 channel for amplitude')
elseif length(phaseCh) == 1 && length(ampCh) == 1
    disp('Computing GMI for only 1 channel')
end

if ismember(foldername,tracking.folders)
    
    if phaseCh == ampCh
        lfp = bz_GetLFP(ampCh,'restrict',timestamps);
    end
    
    cd(foldername)
    if ~isempty('*Tracking.Behavior.mat')
        file = dir('*Tracking.Behavior.mat');
        load(file.name)
    end
    % Speed Limits
    v = tracking.position.speed;
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
    MINPFSIZE = tracking.samplingRate*2;
    % We need to make timestamps
    tracking.timestamps = tracking.timestamps + MergePoints.timestamps(file_number,1);
    t = lfp.timestamps;
    [n2,bin2] = histc(tracking.timestamps,t);
    
    % Epoch obtention
    for jj = 1:max(vvv(:))
        ttt = sum(vvv(:) == jj);
        [rb,cb] = find(bwlabel(vvv) == jj);
        if ttt > MINPFSIZE
            time_tracking = tracking.timestamps(cb(1):cb(end));
            lfp_epoch = lfp.data(bin2(cb(1)):bin2(cb(end)));
            time_epoch = lfp.timestamps(bin2(cb(1)):bin2(cb(end)));
            [Eband,RatTD,pCA1,fCA1] = bz_MakPower2Chronux(lfp_epoch);
            
            if RatTD > thres
                [DegMod,DegreeS,FrecG,GamPower,PowLoc3,ThetFr,GammaFr,AmpT,AmpG] = bz_HilbRelated(lfp_epoch',lfp_epoch',time_epoch','amprange',amprange);
                DegMod_epoch{jj} = DegMod;
                DegreeS_epoch{jj} = DegreeS;
                FrecG_epoch{jj} = FrecG;
                GamPower_epoch{jj} = GamPower;
                PowLoc3_epoch{jj} = PowLoc3;
                ThetFr_epoch{jj} = ThetFr;
                GammaFr_epoch{jj} = GammaFr;
                AmpT_epoch{jj} = AmpT;
                AmpG_epoch{jj} = AmpG;
                
                % ---- HilbRelated_abad variables ----
                DegMod_M = cell2mat(DegMod_epoch);
                DegreeS_M = cell2mat(DegreeS_epoch);
                FrecG_M = cell2mat(FrecG_epoch); %Gamma Frequencies
                GamPower_M = cell2mat(GamPower_epoch);
                PowLoc3_M = cell2mat(PowLoc3_epoch);
                ThetFr_M =mean(cell2mat(ThetFr_epoch));
                GammaFr_M =mean(cell2mat(GammaFr_epoch));
                AmpT_M = mean(cell2mat(AmpT_epoch));
                AmpG_M = mean(cell2mat(AmpG_epoch));
                
                % Hilbert transform plotting and modulation index calculation
                intSize=18;
                [GMI, GFMI, MeanFGa, MeanPowGa, MeanFStd, MeanPowStd, MaxF, MaxFLoc, MaxP, MaxPLoc, VarRag2, bindeRad, MaxdegG] = bz_HistHilbert(GamPower_M,FrecG_M,GammaFr_epoch,DegreeS_M,PowLoc3_M,intSize);
                NoA=(find(isnan(VarRag2)));
                VarRag2(NoA)=0;% eliminate non numeric values result of bining
                % Circular Statistics
                % Preparing data for Circular statistic
                GMI_out{jj} = GMI;
                GFMI_out{jj} = GFMI;
                MeanFGa_out{jj} = MeanFGa;
                MeanPowGa_out{jj} = MeanPowGa;
                MeanFStd_out{jj} = MeanFStd;
                MeanPowStd_out{jj} = MeanPowStd;
                MaxF_out{jj} = MaxF;
                MaxFLoc_out{jj} = MaxFLoc;
                MaxP_out{jj} = MaxP;
                MaxPLoc_out{jj} = MaxPLoc;
                VarRag2_out{jj} = VarRag2;
                bindeRad_out{jj} = bindeRad;
                MaxdegG_out{jj} = MaxdegG;
                
                binde=rad2deg(bindeRad )+180;
                edges=(-pi:0.3307:pi);
                data2=deg2rad(cell2mat(PowLoc3_epoch)-180);
                hist=histc(data2,edges); 
                hist(end)=[];
                [pval{jj},z{jj}]=circ_rtest((edges(1:end-1)+edges(2:end))/2,hist);  %% 
                [mu ul ll] = circ_mean((edges(1:end-1)+edges(2:end))/2,hist,2);
                mvl{jj} = circ_r((edges(1:end-1)+edges(2:end))/2,hist,[],2);% mean vector length
                mu_out{jj} = circ_rad2ang(mu)+180;
                ul_out{jj} = circ_rad2ang(ul)+180;
                ll_out{jj} = circ_rad2ang(ll)+180;

            end
        end
    end

    %% No se si esto está bien
    pval_ = cell2mat(pval);
    z_ = cell2mat(z);
    mvl_ = cell2mat(mvl);
    mu_ = cell2mat(mu_out);
    ul_ = cell2mat(ul_out);
    ll_ = cell2mat(ll_out);
    
    pval_mean = mean(pval_);
    z_mean = mean(z_);
    mvl_mean = mean(mvl_);
    mu_mean = mean(mu_);
    ul_mean = mean(ul_);
    ll_mean = mean(ll_);
    
    GMI_ = cell2mat(GMI_out); % 
    GFMI_ = cell2mat(GFMI_out); % 
    MeanFGa_ = cell2mat(MeanFGa_out); % 
    MeanPowGa_ = cell2mat(MeanPowGa_out);
    MeanFStd_ = cell2mat(MeanFStd_out);
    MeanPowStd_ = cell2mat(MeanPowStd_out);
    MaxF_ = cell2mat(MaxF_out);
    MaxP_ = cell2mat(MaxP_out);
    
    MaxF_mean = mean(MaxF_);
    MaxP_mean = mean(MaxP_);
    index = length(MeanFGa_out{1}); 
    
    a = 1;
    b = 21;
    clear MeanFGa_out
    clear MeanPowGa_out
    clear MeanPowStd_out
    clear MeanFStd_out
    for i = 1:(length(MeanFGa_)/21)
        MeanFGa_out(i,1:21) = MeanFGa_(a:b);
        MeanPowGa_out(i,1:21) = MeanPowGa_(a:b);
        MeanPowStd_out(i,1:21) = MeanPowStd_(a:b);
        MeanFStd_out(i,1:21) = MeanFStd_(a:b);
    %     VarRag2_out(i,1:21) = VarRag2_(a:b);
    %     bindeRad_out(i,1:21) = bindeRad_(a:b);
         a = a+21;
         b = b+21;
    end
    
    intSize = 18;
    binde = 0:intSize:360;
    
    GMI_mean = mean(GMI_);
    GFMI_mean = mean(GFMI_);
    MeanFGa_mean = mean(MeanFGa_out);
    MeanFStd_mean = mean(MeanFStd_out);
    MeanPowGa_mean = mean(MeanPowGa_out);
    MeanPowStd_mean = mean(MeanPowStd_out);
    figure,
    subplot(2,1,1); errorbar(binde,MeanPowGa_mean,MeanPowStd_mean,'LineWidth',2),
    title(['Gamma Power/Theta Modulation ',' GMI  ',num2str(GMI_mean)]);
    
end

%% Filter signal for theta 4-12 Hz
filtered_phase = bz_Filter(lfp,'passband',phaserange,'filter',filterType,'order',filterOrder,'channels',phaseCh);
sigphase = filtered_phase.phase;
gradthe = rad2deg(sigphase);
gradthe2 = gradthe+180;

filtered_amp = bz_Filter(lfp,'passband',amprange,'filter',filterType,'order',filterOrder,'channels',ampCh);
sigamp = filtered_amp.data;


%% Find peaks
for ii=1:length(filtered_phase.channels)
    [mmth{ii} ppth{ii}] = findpeaks(filtered_phase.data(:,ii)); % theta
end

for jj=1:length(filtered_amp.channels)
    [mm{jj} pp{jj}] = findpeaks(filtered_amp.data(:,jj),'minpeakheight',mean(filtered_amp.data(:,jj)));
end

if debugMode
    for i=1:length(filtered_phase.channels)
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        plot(lfp.timestamps,lfp.data(:,i),'k')
        hold on
        plot(filtered_amp.timestamps,filtered_amp.data(:,i),'r')
        hold on
        plot(filtered_amp.timestamps(pp{i}),filtered_amp.data(pp{i},i),'x')
        hold on
        plot(filtered_phase.timestamps,filtered_phase.data(:,i),'b')
        hold on
        plot(filtered_phase.timestamps(ppth{i}),filtered_phase.data(ppth{i},i),'x');
        xlim([0 1])
    end
end

%% Correction of the distance between theta peaks
for ii=1:length(ppth)
    Vft{ii} = diff(sort(ppth{ii})); % distance between theta detected peaks
    inTi1 = find(Vft{ii} >= 1/phaserange(1)*sr); % Remove those peaks that have a separation of more than 2.5 seconds (1/4 Hz)
    Vft{ii}(inTi1)=[];
    inTs1 = find(Vft{ii}<=(1/phaserange(end)*sr)); % Remove those peaks that have a separation of less than 0.08 seconds (1/12 Hz)
    Vft{ii}(inTs1)=[];

    ThetFr{ii} = mean(1./((Vft{ii}) / sr)); % Mean of the frequency at which theta peks occur
    AmpT{ii} = mean(filtered_phase.data(ppth{ii},ii)); % Mean valur of the theta amplitude at the peaks
end
%% Correction of the distance between gamma peaks
for ii=1:length(pp)
    Vfg{ii} = diff(sort(pp{ii}));
    inGi1 = find(Vfg{ii}>=(round(1/amprange(1)*sr)));
    Vfg{ii}(inGi1) = [];
    inGs1 = find(Vfg{ii}<=(1/amprange(end)*sr));
    Vfg{ii}(inGs1) = [];

    GammFr{ii} = mean(1./((Vfg{ii})/sr)); % Mean of the frequency at which gamma peaks occur
    AmpG{ii} = mean(filtered_amp.data(pp{ii},ii)); % Mean value of the gamma amplitude at the peaks
end

%% Some other variables
for ii=1:length(pp)
    DegreeS{ii} = gradthe2(sort(pp{ii}(2:end,1))); % Degrees of theta in which gamma peaks happen
    GamPower{ii} = filtered_amp.data(sort(pp{ii}),ii)/sr; % Gamma amplitude
    FrecG{ii} = 1./(diff(sort(pp{ii}))/sr);
    [MaxFLoc{ii}, MaxF{ii}] = max(FrecG{ii});
    [MaxP{ii},MaxPLoc{ii}] = max(GamPower{ii});

    binde2 = 1:10:360; % Vector from 0 to 360 in intervals of 10
    binde = linspace(1,360,50); % Vector from 0 to 360 in 50 bins
    [DegMod{ii} Val1{ii}]=histc(gradthe2(:,ii),binde); % bin of all angles obtained from theta peaks with the degrees vector
end

%% Theta cicle breakdown
for jj=1:length(filtered_phase.channels)
    [Negmmth{jj} Negppth{jj}] = findpeaks(-filtered_phase.data(:,jj));
    
    if debugMode
        figure;
        set(gcf,'Position',get(0,'ScreenSize'))
        plot(filtered_phase.timestamps,filtered_phase.data(:,jj),'k');title('Theta filtered signal');
        text(filtered_phase.timestamps(Negppth{jj}),filtered_phase.data(Negppth{jj},jj),'x');
        xlim([0 1])
    end

    Sorth{jj} = sort(Negppth{jj});
    for ii=1:length(Sorth{jj})-1

       sigth = filtered_phase.data(Sorth{jj}(ii): Sorth{jj}(ii+1),jj);
       siggam = filtered_amp.data(Sorth{jj}(ii): Sorth{jj}(ii+1),jj);
       siggradthe2 = gradthe2(Sorth{jj}(ii): Sorth{jj}(ii+1),jj);

       [MaxP  PowLoc]=max(siggam); 
       PowLoc2{ii}=PowLoc;

    end

    PowLoc3{jj} = cell2mat(PowLoc2);
end

%% HIST HILBERT
for ii=1:length(DegreeS)
    [nn2{ii},bbin2{ii}] = histc(DegreeS{ii},binde);

    for jkl=1:length(binde)
        % Colocamos cada ángulo de theta de un pico gamma detectado (juntados
        % los epochs), dentro de cada uno de los 21 bines. Y agrupamos los
        % ángulos de theta de acuerdo a los bines correspondientes
        inx = find(bbin2{ii}==jkl);
        MeanPow = GamPower{ii}(inx); % Gamma 
        MeanPowGa{ii}(jkl) = mean(MeanPow);
        MeanPowStd{ii}(jkl) = std(MeanPow)/sqrt(length(MeanPow)) ;
        MeanF = FrecG{ii}(inx);
        MeanFStd{ii}(jkl) = std( MeanF)/sqrt(length( MeanF)) ;
        MeanFGa{ii}(jkl) = mean(MeanF);
    end


    GMI{ii} = (max(MeanPowGa{ii})-min(MeanPowGa{ii}))/(max(MeanPowGa{ii})+min(MeanPowGa{ii})); % Modulation index as in Monyer
    GFMI{ii} = (max(MeanFGa{ii})-min(MeanFGa{ii}))/(max(MeanFGa{ii})+min(MeanFGa{ii})); %

    [MaxF, MaxFLoc] = max(MeanFGa{ii});
    [MaxP, MaxPLoc] = max(MeanPowGa{ii});

    % Data for Rayleight
    bindeRad = deg2rad(binde-180);
    VarRag = deg2rad(PowLoc3{ii}-180);
    [Radnn2,Radbbin2] = histc(VarRag,bindeRad); % bins the spike time of all timestams unit to the position time interval

    for jkll=1:length(bindeRad)
        inxR = find(Radbbin2 == jkll);
        VarRag2(jkll) = mean(VarRag(inxR));
    end
end


%% Circular statistics
clear p
for ii=1:length(PowLoc3)
    % Preparaing data for circular statistics
    binde = rad2deg(bindeRad)+180;
    % edges=(-pi:0.3307:pi);
    edges = linspace(-pi,pi,50);
    data2 = deg2rad(PowLoc3{ii}-180);
    hist = histc(data2,edges); 
    hist(end) = [];

    [p{ii},z{ii}]=circ_rtest((edges(1:end-1)+edges(2:end))/2,hist);
    % figure,
    % polarplot(hist);
    % figure,
    % rose(PowLoc3,18)


    % mean circular 
    [mu{ii} ul{ii} ll{ii}] = circ_mean((edges(1:end-1)+edges(2:end))/2,hist,2);
    mvl{ii} = circ_r((edges(1:end-1)+edges(2:end))/2,hist,[],2);% mean vector length
    mu_out{ii} = circ_rad2ang(mu{ii})+180;
    ul_out{ii} = circ_rad2ang(ul{ii})+180;
    ll_out{ii} = circ_rad2ang(ll{ii})+180;
end



%% PLOTTING

if showFig
    
    if length(phaseCh) == 1 && length(ampCh) == 1
    
        figure
        set(gcf,'Position',get(0,'ScreenSize'))
        subplot(2,1,1);errorbar( binde, MeanPowGa{1},MeanPowStd{1} ,'LineWidth',2);title(['Gamma Power/Theta Modulation ',' GMI  ',num2str(GMI{1}),' Ch: ', num2str(phaseCh)]);
        axis([0 380 min(MeanPowGa{1})-max(MeanPowStd{1})*2  max(MeanPowGa{1})+max(MeanPowStd{1})*2 ]);
        subplot(2,1,2);errorbar( binde,MeanFGa{1}, MeanFStd{1},'LineWidth',2 );title(['Gamma Frequency/Theta Modulation ',' GFMI  ',num2str(GFMI{1}),' Ch: ', num2str(phaseCh)]);
        axis([0 380 min(MeanFGa{1})-max(MeanFStd{1})*2  max(MeanFGa{1})+max(MeanFStd{1})*2  ]);
        if saveFig
            if ~isempty(foldername)
                saveas(gcf,['lfpAnalysisFigures\GMI.',foldername,'_',num2str(amprange(1)),num2str(amprange(end)),'.png'])
            else
                saveas(gcf,['lfpAnalysisFigures\GMI','_',num2str(amprange(1)),num2str(amprange(end)),'.png'])
            end
        end
    else
        meanpowga = cell2mat(MeanPowGa);
        stdpowga = cell2mat(MeanPowStd);
        minpowga = min(meanpowga);
        maxpowga = max(meanpowga);
        minpowstd = min(stdpowga);
        maxpowstd = max(stdpowga);
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for jj=1:size(shanks,2)
            count = jj;
            for ii=1:length(shanks(jj).Channels) 
                subplot(sessionInfo.nChannels/sessionInfo.nElecGps,sessionInfo.nChannels/sessionInfo.nElecGps,count)
                errorbar(binde,MeanPowGa{shanks(jj).Channels(ii)+1},MeanPowStd{shanks(jj).Channels(ii)+1} );
                axis([0 380 min(MeanPowGa{shanks(jj).Channels(ii)+1})-max(MeanPowStd{shanks(jj).Channels(ii)+1})*2  max(MeanPowGa{shanks(jj).Channels(ii)+1})+max(MeanPowStd{shanks(jj).Channels(ii)+1})*2 ]);
                axis([0 380 minpowga-maxpowstd*2 maxpowga+maxpowstd*2])                
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
        
        if saveFig
            if ~isempty(foldername)
                saveas(gcf,['lfpAnalysisFigures\GMI_allCh.',foldername,'_',num2str(amprange(1)),num2str(amprange(end)),'.png'])
            else
                saveas(gcf,['lfpAnalysisFigures\GMI_allCh','_',num2str(amprange(1)),num2str(amprange(end)),'.png'])
            end
        end
        
        
    end    
end




%% save matlab structure
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


