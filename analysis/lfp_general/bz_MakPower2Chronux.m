function [Eband, RatTD, p,f] = bz_MakPower2Chronux(data,varargin)

%% Default Params
p = inputParser;

addParameter(p,'fs',1250,@isnumeric);
addParameter(p,'Band',[0 200], @isnumeric);
addParameter(p,'showFig',false,@islogical);

parse(p,varargin{:})

fs = p.Results.fs;
Band = p.Results.Band;
showFig = p.Results.showFig;

params = struct('Fs',fs,'pad',1,'fpass',Band,'err',[1 0.05]);

[p,f] = mtspectrumc(data,params);

FreqBands = [1 3; 4 12; 13 16; 17 29; 30 65; 66 130; 150 185];
nBands = size(FreqBands,1);
Eband = zeros(1,nBands);

for k = 1:nBands
    iDx = f >= FreqBands(k,1) & f < FreqBands(k,2);
    Eband(k) = 1*trapz(f(iDx),p(iDx)); % trapz function obtains the integral, (area) below the periodogram for each band
end

RatTD = Eband(2)/Eband(1);


end

