function [GMI, GFMI, MeanFGa, MeanPowGa,  MeanFStd,   MeanPowStd, MaxF,  MaxFLoc, MaxP, MaxPLoc,  VarRag2, bindeRad, MaxdegG] = bz_HistHilbert(GamPowerM,FrecGM,GammaFQ,DegreeS,PowLoc44,intSize,varargin)

% Defaults and Params

p = inputParser;

addParameter(p,'fs',1250,@isnumeric);

parse(p,varargin{:})

fs = p.Results.fs;

%%
binde=0:intSize:360;
[nn2,bbin2] = histc(DegreeS,binde);

for jkl=1:length(binde)
    % Colocamos cada ángulo de theta de un pico gamma detectado (juntados
    % los epochs), dentro de cada uno de los 21 bines. Y agrupamos los
    % ángulos de theta de acuerdo a los bines correspondientes
    inx = find(bbin2==jkl);
    MeanPow=GamPowerM(inx); % Gamma 
    MeanPowGa(jkl) = mean(MeanPow);
    MeanPowStd(jkl)= std(MeanPow)/sqrt(length(MeanPow)) ;
%     MeanPowStd(jkl) = std(MeanPow); %edited by pablo
    MeanF = FrecGM(inx);
    MeanFStd(jkl)=std( MeanF)/sqrt(length( MeanF)) ;
    MeanFGa(jkl)=mean(MeanF);
end

GMI=(max(MeanPowGa)-min(MeanPowGa))/(max(MeanPowGa)+min(MeanPowGa)); % Modulation index as in Monyer
GFMI=(max(MeanFGa)-min(MeanFGa))/(max(MeanFGa)+min(MeanFGa)); %
% Modulation index s in Monyer

[MaxF, MaxFLoc]=max(MeanFGa);
[MaxP, MaxPLoc]=max(MeanPowGa);
MaxdegG=binde(MaxPLoc);

%% Data for Rayleight
bindeRad=deg2rad(binde-180);
VarRag=deg2rad(PowLoc44)-180;
[Radnn2,Radbbin2] = histc(VarRag,bindeRad); % bins the spike time of all timestams unit to the position time interval

for jkll=1:length(bindeRad);
    inxR=find(Radbbin2==jkll);
      VarRag2(jkll)=mean(VarRag(inxR));
   

end

end

