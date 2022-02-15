function [] = referenceChannelFromDat(basepath,varargin)
% Edit dat files with different options
%
% USAGE
%   referenceChannelFromDat(basepath,varargin)
%
% INPUT
%   basepath    If not provided, takes pwd
%   refCh       refCh to work as reference channel to substract to all of
%                   the rest
%   keepDat     Default, true.
%
%
% Pablo Abad - December 2021, based on Manuel Valero-BuzsakiLab 2021
%               removeNoiseFromDat
%
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'ch','all',@isnumeric);
addParameter(p,'refCh',[],@isnumeric);
addParameter(p,'keepDat',false,@islogical);

warning('Performing reference Channel substraction !! Dat file will be compromised!!');
parse(p,varargin{:})
basepath = p.Results.basepath;
ch = p.Results.ch;
refCh = p.Results.refCh;
keepDat = p.Results.keepDat;

% Get elements
prevPath = pwd;
cd(basepath);

% Get sessionInfo and session
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
session = sessionTemplate(pwd,'showGUI',false);

fileTargetAmplifier = dir('amplifier*.dat');
if isempty(fileTargetAmplifier)
    filename = split(pwd,filesep); filename = filename{end};
    fileTargetAmplifier = dir([filename '*.dat']);
end

if size(fileTargetAmplifier,1) == 0
    error('Dat file not found!!');
end

if ischar(ch) && strcmpi(ch,'all')
    ch = 1:length(sessionInfo.channels);
end

nChannels = sessionInfo.nChannels;
duration = 60;
frequency = sessionInfo.rates.wideband;
fid = fopen(fileTargetAmplifier.name,'r'); filename = fileTargetAmplifier.name;
C = strsplit(fileTargetAmplifier.name,'.dat'); filenameOut = [C{1} '_temp.dat'];
fidOutput = fopen(filenameOut,'a');

while 1
    data = fread(fid,[nChannels frequency*duration],'int16');
    if isempty(data)
        break;
    end
    
        m_data = (data(refCh,:));

    for ii = 1:length(ch)
        data(ch(ii),:) = int16(double(data(ch(ii),:)) - double(m_data));
    end
    fwrite(fidOutput,data,'int16');
end
fclose(fid);
fclose(fidOutput);

if ~keepDat
    copyfile(filename, [C{1} '_original.dat']);
end
delete(filename);
movefile(filenameOut, filename);

cd(prevPath);



end 

