function removeNoiseFromDat_subfolders(basepath,varargin)
% Edit dat files with different options.
%
% USAGE
%   editDatFile(basepath,ints,varargin)
% 
% INPUT
% basepath      If not provided, takes pwd
% threshold     Intervals to be edited.
% method        'substractMedian' or 'substratMean' (defaut)
% keepDat       Default, false.
%
% <optional>
% option        'remove' or 'zeroes' (default). 
%
% Manu Valero-BuzsakiLab 2021
%
%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'ch','all');
addParameter(p,'method','substractMean',@ischar);
addParameter(p,'keepDat',false,@islogical);

warning('Performing median/mean substraction!! Dat file will be compromised!! ');
parse(p,varargin{:});
ch = p.Results.ch;
method = p.Results.method;
basepath = p.Results.basepath;
keepDat = p.Results.keepDat;

% Get elements
prevPath = pwd;
cd(basepath);

xml = LoadParameters;
frequency = xml.rates.wideband;
session = loadSession();
sessionInfo = bz_getSessionInfo();

if exist([basepath filesep strcat(session.general.name,'.MergePoints.events.mat')],'file')
	load(strcat(session.general.name,'.MergePoints.events.mat'));
end

if ischar(ch) && strcmpi(ch, 'all')
    ch = 1:length(xml.channels);
end
nChannels = xml.nChannels;


for i=1:length(MergePoints.foldernames)
    cd(MergePoints.foldernames{i})
    fileTargetAmplifier = dir('amplifier*.dat');
    if isempty(fileTargetAmplifier)
        filename = split(pwd,filesep); 
        filename = filename{end};
        fileTargetAmplifier = dir([filename '*.dat']);
    end
    if size(fileTargetAmplifier,1) == 0
        error('Dat file not found!!');
    end
    disp('Creating .dat back up...')
    copyfile(strcat(filename,'.dat'),strcat(filename,'_original.dat'));
    
    tic
    m = memmapfile(fullfile(basepath,strcat(filename,'.dat')),'Format','int16','Writable',true);
    data = reshape(m.Data,nChannels,[]);
    toc
    
    fid = fopen(fileTargetAmplifier.name,'r'); 
    filename = fileTargetAmplifier.name;
    C = strsplit(fileTargetAmplifier.name,'.dat'); 
    filenameOut = [C{1} '_temp.dat'];
    fidOutput = fopen(filenameOut,'a');
    
    
end
    
    
    











while 1
    data = fread(fid,[nChannels frequency*duration],'int16');
    
    if isempty(data)
        break;
    end
    
    if strcmpi('substractMedian',method)
        m_data = median(data(ch,:));
    elseif strcmpi('substractMean',method)
        m_data = mean(data(ch,:));
    end

    for ii = 1:length(ch)
        data(ch(ii),:) = int16(double(data(ch(ii),:)) - double(m_data));
    end
    fwrite(fidOutput,data,'int16');
end
fclose(fid);
fclose(fidOutput);

% if ~keepDat
%     copyfile(filename, [C{1} '_original.dat']);
% end
% delete(filename);
% movefile(filenameOut, filename);

cd(prevPath);

end

