function [DegMod,DegreeS,FrecG,GamPower,PowLoc3,ThetFr,GammaFr,AmpT,AmpG] = bz_HilbRelated(phase,amp,t,varargin)
%% Default and Params
p = inputParser;

addParameter(p,'fs',1250,@isnumeric);
addParameter(p,'amprange',[],@isnumeric);
addParameter(p,'phaserange',[4 12], @isnumeric);
addParameter(p,'verbose',false,@islogical);

parse(p,varargin{:})

fs = p.Results.fs;
amprange = p.Results.amprange;
phaserange = p.Results.phaserange;
verbose = p.Results.verbose;

%% Phase Filter 
BC = fir1(200,[phaserange(1)/(fs/2) phaserange(2)/(fs/2)],'bandpass');
filteredPhase = filtfilt(BC,[1],phase);
%% Amplitude Filter
B = fir1(400,[amprange(1)/(fs/2) amprange(2)/(fs/2)],'bandpass');
filteredAmp = filtfilt(B,[1],amp);


%% Hilber transform for phase signal
filteredPhase_hilbert = rad2deg(angle(hilbert(filteredPhase)))+180;

% Phase peaks for frequency calculation
[mmth ppth] = findpeaks(filteredPhase);
% Amp peaks
[mm pp] = findpeaks(filteredAmp,'minpeakheight',mean(filteredAmp));
% Phase Correction
Vft = diff(sort(ppth));
inTi1=find(Vft>=(1/phaserange(1)*fs));
inTs1=find(Vft<=(1/phaserange(2)*fs));
inT = [inTi1 inTs1];
% Vft(inTi1) = [];
% Vft(inTs1)=[];
Vft(inT) = [];
ThetFr=mean(1./((Vft) /fs ));
AmpT=mean(filteredPhase(ppth));

% Amp Correction
Vfg=diff(sort(pp));
inGi1=find(Vfg>=(round(1/amprange(1)*fs)));
inGs1=find(Vfg<=(round(1/amprange(2)*fs)));
% inG = [inGi1; inGs1];
inG = [inGi1, inGs1];
% Vfg(inGi1)=[];
% Vfg(inGs1)=[];
Vfg(sort(inG)) = [];

if verbose % for plotting
    figure;
    plot(t,phase);hold
    plot(t,filteredAmp,'r'); 
    text(t(pp),filteredAmp(pp),'x')
    plot(t,filteredPhase,'k');title('Theta(black) and gamma (red) filtered signal and peak detention');
    text(t(ppth),filteredPhase(ppth),'t')
    %close
end

GammaFr=mean(1./((Vfg)/fs));
AmpG=mean(filteredAmp(pp));

DegreeS=filteredPhase_hilbert(sort(pp(1,2:end)));
GamPower=filteredAmp(sort(pp(1,2:end)));
FrecG = 1./(diff(sort(pp))/fs);
[MaxFLoc,  MaxF]=max(FrecG);
binde = 1:10:360;
[DegMod Val1]=histc(filteredPhase_hilbert,binde);

%% Phase cycle breakdown
[Negmmth Negppth] = findpeaks(-filteredPhase);
if verbose
    figure;plot(t,filteredPhase,'k');title('Theta filtered signal');
    text(t(Negppth),filteredPhase(Negppth),'x');
end

Sorth=sort(Negppth);
for ijkl=1:length(Sorth)-1
    sigth=filteredPhase(Sorth(ijkl): Sorth(ijkl+1));
    siggam=filteredAmp(Sorth(ijkl): Sorth(ijkl+1));
    siggradthe2=filteredPhase_hilbert(Sorth(ijkl): Sorth(ijkl+1));
    [MaxP  PowLoc]=max(siggam); % text(PowLoc,sigth( PowLoc),'X?')
    PowLoc2{ijkl}=PowLoc;   
end

PowLoc3=cell2mat(PowLoc2);

end

