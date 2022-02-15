function [performance] = getSessionPerformance(varargin)
%
% Gets performance for each sub-session that contains tracking (only for
% TMaze and YMaze) and concatenate all of them so they are aligned with LFP
% and spikes.
%
% Performance can be calculated with the digitalIn from AnyMaze, but we
% need to create an alternative because in some experiments the digital
% Inputs were not working.
%
% USAGE
%
%   [performance] = getSessionPerformance(varargin)
%
% INPUTS
%   basePath       -(default: pwd) basePath for the recording file, in buzcode format:
%   saveMat        - default true
%   forceReload    - default false
%
% OUTPUT
%       - Performance.behaviour output structure, with the fields:
%   
%   
%   
%   
%   
%   
%   
%   
%   
%    

%   HISTORY:
%     - Pablo Abad 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Arm 1 (main): Digital Input 3
% - Arm 2 (right): Digital Input 4
% - Arm 3 (left) : Digital Input 5
% - Center (center): Digital Input 6

%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'forceReload',false,@islogical)
addParameter(p,'tracking',[],@isstruct);



parse(p,varargin{:});
basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;
tracking = p.Results.tracking;

if isempty(tracking)
    tracking = getSessionTracking();
end

%% In case performance already exists 
if ~isempty(dir([basepath filesep '*Performance.Behavior.mat'])) || forceReload
    disp('Performance already detected! Loading file.');
    file = dir([basepath filesep '*Performance.Behavior.mat']);
    load(file.name);
    return
end

[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);

%% Detect YMaze or TMaze and analyze the performance
if exist([basepath filesep strcat(sessionInfo.session.name, '.MergePoints.events.mat')],'file')
    load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
end

count = 0;
count1 = 0;
for ii=1:length(MergePoints.foldernames)
    if ismember(MergePoints.foldernames{ii},tracking.folders)
        count1 = count1+1;
        % YMaze Performance
        if strcmpi(tracking.apparatus{count1}.name,'YMaze Apparatus') || strcmpi(tracking.apparatus{count1}.name,'YMaze')
            
            count = count+1;
            cd(tracking.folders{count1})
            digitalIn = bz_getDigitalIn('all');
            tracking_subfolder = getSessionTracking();
            try
               fprintf('Computing YMaze performance in %s folder \n',tracking.folders{count1});
               tempPerformance{count} =  performanceYMaze('digitalIn',digitalIn,'tracking',tracking_subfolder);
               performanceFolder(count) = ii; 

            catch
                error('It is not possible to run Performance code...')

            end

            cd(basepath)
        
            % TMaze Performance
        elseif strcmpi(tracking.apparatus{count1}.name,'TMaze')
            count = count+1;
            cd(tracking.folders{count1})
            digitalIn = bz_getDigitalIn('all');
            tracking_subfolder = getSessionTracking();
            try
                fprintf('Computing TMaze performance in %s folder \n',tracking.folders{count1});
                tempPerformance{count} = performanceTMaze('digitalIn', digitalIn, 'tracking', tracking_subfolder);
                performanceFolder(count) = ii;
            catch
                
            end
            cd (basepath)
            
%         elseif strcmpi(tracking.apparatus{count1}.name,'Open Field')
%             count = count+1;
%             cd(tracking.folders{count1})
%             digitalIn = bz_getDigitalIn('all');
%             tracking_subfolder = getSessionTracking();
%             try
%                 fprintf('Computing Open Field Performance in %s folder \n', tracking.folders{count1});
%                 tempPerformance{count} = performanceOpenField('digitalIn', digitalIn, 'tracking', tracking_subfolder);
%                 performanceFolder(count) = ii;
%             catch
% 
%             end
%             
%             cd(basepath)
        elseif strcmpi(tracking.apparatus{count1}.name,'Object Recognition')
            count = count+1;
            cd(tracking.folders{count1})
            digitalIn = bz_getDigitalIn('all');
            tracking_subfolder = getSessionTracking();
            try
                fprintf('Computing TMaze performance in %s folder \n',tracking.folders{count1});
                tempPerformance{count} = performanceObjectRecognition('digitalIn',digitalIn,'tracking',tracking_subfolder);
                performanceFolder(count) = ii;
            catch
            end
            
            cd(basepath)
        end
    end

end

%% Concatenate and sync timestamps

ts = []; subSessions = []; maskSessions = [];
% YMaze variables
ts_transitions = [];
ts_transitions_zero = [];
ts_right = [];
ts_wrong = [];
group1_sample = []; group1_choice = [];
group2_sample = []; group2_choice = [];
group3_sample = []; group3_choice = [];
group4_sample = []; group4_choice = [];
times_stemArm = [];
times_rightArm = [];
times_leftArm = [];

% Object Recognition variables

transition = [];
ts_Object1 = [];
ts_Object2 = [];



for ii=1:length(performanceFolder)
    if strcmpi(MergePoints.foldernames{performanceFolder(ii)},tempPerformance{ii}.folder)
        if strcmpi(tempPerformance{ii}.paradigm,'YMaze') 
            sumTs{ii} = tempPerformance{ii}.timestamps + MergePoints.timestamps(performanceFolder(ii),1);
            sumTs_transitions{ii} = tempPerformance{ii}.transitions_times_out_nonzero + MergePoints.timestamps(performanceFolder(ii),1);
            sumTs_transitions_zero{ii} = tempPerformance{ii}.transitions_times_out + MergePoints.timestamps(performanceFolder(ii),1);
            subSessions = [subSessions; MergePoints.timestamps(performanceFolder(ii),1:2)];
            maskSessions = [maskSessions; ones(size(sumTs{ii}))*ii];
            ts_stemArm{ii} = tempPerformance{ii}.times_stem' + MergePoints.timestamps(performanceFolder(ii),1);
            ts_leftArm{ii} = tempPerformance{ii}.times_left' + MergePoints.timestamps(performanceFolder(ii),1);
            ts_rightArm{ii} = tempPerformance{ii}.times_right' + MergePoints.timestamps(performanceFolder(ii),1);


            ts = [ts, sumTs{ii}'];
            ts_transitions = [ts_transitions; sumTs_transitions{ii}'];
            ts_transitions_zero = [ts_transitions_zero; sumTs_transitions_zero{ii}'];
            times_stemArm = [times_stemArm; ts_stemArm{ii}'];
            times_leftArm = [times_leftArm; ts_leftArm{ii}'];
            times_rightArm = [times_rightArm; ts_rightArm{ii}'];

            for i=1:tempPerformance{ii}.sampleVSchoice.group1.trials
                group1_sample{ii}{i} = tempPerformance{ii}.sampleVSchoice.group1.sample{i}(3,:) + MergePoints.timestamps(performanceFolder(ii),1);
                group1_choice{ii}{i} = tempPerformance{ii}.sampleVSchoice.group1.choice{i}(3,:) + MergePoints.timestamps(performanceFolder(ii),1);
            end

            for i=1:tempPerformance{ii}.sampleVSchoice.group2.trials
                group2_sample{ii}{i} = tempPerformance{ii}.sampleVSchoice.group2.sample{i}(3,:) + MergePoints.timestamps(performanceFolder(ii),1);
                group2_choice{ii}{i} = tempPerformance{ii}.sampleVSchoice.group2.choice{i}(3,:) + MergePoints.timestamps(performanceFolder(ii),1);
            end

            for i=1:tempPerformance{ii}.sampleVSchoice.group3.trials
                group3_sample{ii}{i} = tempPerformance{ii}.sampleVSchoice.group3.sample{i}(3,:) + MergePoints.timestamps(performanceFolder(ii),1);
                group3_choice{ii}{i} = tempPerformance{ii}.sampleVSchoice.group3.choice{i}(3,:) + MergePoints.timestamps(performanceFolder(ii),1);
            end

            for i=1:tempPerformance{ii}.sampleVSchoice.group4.trials
                group4_sample{ii}{i} = tempPerformance{ii}.sampleVSchoice.group4.sample{i}(3,:) + MergePoints.timestamps(performanceFolder(ii),1);
                group4_choice{ii}{i} = tempPerformance{ii}.sampleVSchoice.group4.choice{i}(3,:) + MergePoints.timestamps(performanceFolder(ii),1);
            end

            for i=1:length(tempPerformance{ii}.right_times)
                ts_right{ii}{i} = tempPerformance{ii}.right_times{i} + MergePoints.timestamps(performanceFolder(ii),1);
            end
            
            if iscell(tempPerformance{ii}.wrong_times)
                for i = 1:length(tempPerformance{ii}.wrong_times)
                    ts_wrong{ii}{i} = tempPerformance{ii}.wrong_times{i} + MergePoints.timestamps(performanceFolder(ii),1);
                end
            else
                ts_wrong = NaN;
            end
        
        
        elseif strcmpi(tempPerformance{ii}.paradigm,'Object Recognition')
            sumTs{ii} = tempPerformance{ii}.timestamps + MergePoints.timestamps(performanceFolder(ii),1);
            sumTs_Object1{ii} = tempPerformance{ii}.ts.Object1 + MergePoints.timestamps(performanceFolder(ii),1);
            sumTs_Object2{ii} = tempPerformance{ii}.ts.Object2 + MergePoints.timestamps(performanceFolder(ii),1);
            
            ts = [ts; sumTs{ii}];
            ts_Object1 = [ts_Object1 sumTs_Object1{ii}];
            ts_Object2 = [ts_Object2 sumTs_Object2{ii}];
            
            subSessions = [subSessions; MergePoints.timestamps(performanceFolder(ii),1:2)];
            maskSessions = [maskSessions; ones(size(sumTs{ii}))*ii];
%             entrances_Object1{ii} = tempPerformance{ii}.entrances.Object1;
%             entrances_Object2{ii} = tempPerformance{ii}.entrances.Object2;
%             score_Object1{ii} = tempPerformance{ii}.score.Object1;
%             score_Object2{ii} = tempPerformance{ii}.score.Object2;
%             time_Object1{ii} = tempPerformance{ii}.time.Object1;
%             time_Object2{ii} = tempPerformance{ii}.time.Object2;
            
            transition{ii} = tempPerformance{ii}.transition;
            transition{ii}(3,:) = MergePoints.timestamps(performanceFolder(ii),1) + transition{ii}(3,:);
        end    
    end
end



    %% Concatenating performance fields...

% General Paradigm
folder = []; paradigm = [];

% YMaze Paradigm
score = [];
transitions_out_nonzero = []; transitions_fr_out_nonzero = []; transitions_times_out_nonzero = [];
right_transitions_epoch = []; right_transitions_fr = []; wrong_transitions_epoch = []; wrong_transitions_fr = [];
entrances = []; transitions_out = []; transitions_fr_out = []; transitions_times_out = [];
frames_stem = []; frames_right = []; frames_left = []; times_stem = []; times_right = []; times_left = [];
ballistic = [];  doubted = []; ballistic_perc = []; doubted_perc = []; times = [];
right_epoch = []; right_fr = []; right_times = []; wrong_epoch = []; wrong_fr = []; wrong_times = []; rightTriades = []; wrongTriades = []; sampleVSchoice = []; 

% Object Recognition Paradigm
entrances = [];
score = [];
time = [];
folder = [];
paradigm = [];


for ii=1:size(tempPerformance,2)
    if strcmpi(tempPerformance{ii}.paradigm,'YMaze') 
        score = [score; tempPerformance{ii}.performance];

        transitions_out_nonzero{ii} = [tempPerformance{ii}.transitions_out_nonzero];
        transitions_fr_out_nonzero{ii} = tempPerformance{ii}.transitions_fr_out_nonzero;

        right_transitions_epoch{ii} = tempPerformance{ii}.right_transitions_epoch;
        right_transitions_fr{ii} = tempPerformance{ii}.right_transitions_fr;

        wrong_transitions_epoch{ii} = tempPerformance{ii}.wrong_transitions_epoch;
        wrong_transitions_fr{ii} = tempPerformance{ii}.wrong_transitions_fr;

        entrances{ii} = tempPerformance{ii}.entrances;

        transitions_out{ii} = tempPerformance{ii}.transitions_out;
        transitions_fr_out{ii} = tempPerformance{ii}.transitions_fr_out;

        frames_stem{ii} = tempPerformance{ii}.frames_stem;
        frames_right{ii} = tempPerformance{ii}.frames_right;
        frames_left{ii} = tempPerformance{ii}.frames_left;

        ballistic{ii} = tempPerformance{ii}.ballistic;
        doubted{ii} = tempPerformance{ii}.doubted;
        ballistic_perc{ii} = tempPerformance{ii}.ballistic_perc;
        doubted_perc{ii} = tempPerformance{ii}.doubted_perc;
        times{ii} = tempPerformance{ii}.times;

        right_epoch{ii} = tempPerformance{ii}.right_epoch;
        right_fr{ii} = tempPerformance{ii}.right_fr;
        wrong_epoch{ii} = tempPerformance{ii}.wrong_epoch;
        wrong_fr{ii} = tempPerformance{ii}.wrong_fr;
        rightTriades{ii} = tempPerformance{ii}.rightTriades;
        wrongTriades{ii} = tempPerformance{ii}.wrongTriades;
        sampleVSchoice{ii} = tempPerformance{ii}.sampleVSchoice;
        folder{ii} = tempPerformance{ii}.folder;
        paradigm{ii} = tempPerformance{ii}.paradigm;
        
    elseif strcmpi(tempPerformance{ii}.paradigm,'Object Recognition')
        entrances{ii} = tempPerformance{ii}.entrances;
        score{ii} = tempPerformance{ii}.score;
        time{ii} = tempPerformance{ii}.time;
        folder{ii} = tempPerformance{ii}.folder;
        paradigm{ii} = tempPerformance{ii}.paradigm;
    end
end

if ~isempty(sampleVSchoice)
    for i=1:length(sampleVSchoice)
        if ~isempty(sampleVSchoice{i})
            for j=1:length(sampleVSchoice{i}.group1.sample)
                sampleVSchoice{i}.group1.sample{j}(3,:) = group1_sample{i}{j};
                sampleVSchoice{i}.group1.choice{j}(3,:) = group1_choice{i}{j};
            end
            for j=1:length(sampleVSchoice{i}.group2.sample)
                sampleVSchoice{i}.group2.sample{j}(3,:) = group2_sample{i}{j};
                sampleVSchoice{i}.group2.choice{j}(3,:) = group2_choice{i}{j};
            end
            if sampleVSchoice{i}.group3.trials > 0 
                for j=1:length(sampleVSchoice{i}.group3.sample)
                    sampleVSchoice{i}.group3.sample{j}(3,:) = group3_sample{i}{j};
                    sampleVSchoice{i}.group3.choice{j}(3,:) = group3_choice{i}{j};
                end
            else
                sampleVschoice{i}.group3.sample = NaN;
                sampleVschoice{i}.group3.choice = NaN;
            end
            for j=1:length(sampleVSchoice{i}.group4.sample)
                sampleVSchoice{i}.group4.sample{j}(3,:) = group4_sample{i}{j};
                sampleVSchoice{i}.group4.choice{j}(3,:) = group4_choice{i}{j};
            end
        end
    end
end

%% CREATING RIGHT AND WRONG FROM SAMPLEVSCHOICE

r = [];
w = [];
right = [];
wrong = [];

for ii=1:size(tempPerformance,2)
    if strcmpi(tempPerformance{ii}.paradigm,'YMaze')
    
        % Group1
        for i=1:tempPerformance{ii}.sampleVSchoice.group1.trials

            right_sample_timestamp = sampleVSchoice{ii}.group1.sample{i}(3,2);
            right_choice_timestamp = sampleVSchoice{ii}.group1.choice{i}(3,2);        
            r = [r; right_sample_timestamp; right_choice_timestamp]; 
        end

        % Group2
        for i=1:tempPerformance{ii}.sampleVSchoice.group2.trials
            right_choice_timestamp = sampleVSchoice{ii}.group2.choice{i}(3,2);
            r = [r; right_choice_timestamp];

            wrong_sample_timestamp = sampleVSchoice{ii}.group2.sample{i}(3,2);
            w = [w; wrong_sample_timestamp];
        end

        % Group3
        for i=1:tempPerformance{ii}.sampleVSchoice.group3.trials
            right_sample_timestamp = sampleVSchoice{ii}.group3.sample{i}(3,2);
            r = [r; right_sample_timestamp];

            wrong_choice_timestamp = sampleVSchoice{ii}.group3.choice{i}(3,2);
            w = [w; wrong_choice_timestamp];

        end

        % Group4
        for i=1:tempPerformance{ii}.sampleVSchoice.group4.trials

            wrong_sample_timestamp = sampleVSchoice{ii}.group4.sample{i}(3,2);
            wrong_choice_timestamp = sampleVSchoice{ii}.group4.choice{i}(3,2);        
            w = [w; wrong_sample_timestamp; wrong_choice_timestamp]; 
        end

        right{ii}.timestamps = unique(r);
        wrong{ii}.timestamps = unique(w);
    end
end





%% OUTPUT
performance = [];


performance.score = score;
performance.entrances = entrances;
performance.folder = folder;
performance.paradigm = paradigm;
performance.timestamps = ts;
performance.ts = sumTs;
performance.subSessions = subSessions;
performance.maskSessions = maskSessions;

% Object Recogntion
if exist('time','var')
    performance.time = time;
end
if exist('ts_Object1','var')
    performance.ts_Object1 = ts_Object1;
end
if exist('ts_Object2','var')
    performance.ts_Object2 = ts_Object2;
end
if exist('transition','var')
    performance.transition = transition;
end

% YMaze Performance
if ~isempty(transitions_out_nonzero)
    performance.transitions_out_nonzero = transitions_out_nonzero;
end
if ~isempty(transitions_fr_out_nonzero)
    performance.transitions_fr_out_nonzero = transitions_fr_out_nonzero;
end
if exist('sumTs_transitions','var') 
    performance.transitions_times_out_nonzero = sumTs_transitions;
end
if ~isempty(transitions_out)
    performance.transitions_out = transitions_out;
end
if ~isempty(transitions_fr_out)
    performance.transitions_fr_out = transitions_fr_out;
end
if exist('sumTs_transitions_zero','var')
    performance.transitions_times_out = sumTs_transitions_zero;
end
if ~isempty(right_transitions_epoch)
    performance.right_transitions_epoch = right_transitions_epoch;
end
if ~isempty(right_transitions_fr)
    performance.right_transitions_fr = right_transitions_fr;
end
if ~isempty(right_epoch)
    performance.right_epoch = right_epoch;
end
if ~isempty(right_fr)
    performance.right_fr = right_fr;
end
if ~isempty(ts_right)
    performance.right_times = ts_right;
end
if ~isempty(rightTriades)
    performance.rightTriades = rightTriades;
end
if ~isempty(wrong_transitions_epoch)
    performance.wrong_transitions_epoch = wrong_transitions_epoch;
end
if ~isempty(wrong_transitions_fr)
    performance.wrong_transitions_fr = wrong_transitions_fr;
end
if ~isempty(wrong_epoch)
    performance.wrong_epoch = wrong_epoch;
end
if ~isempty(wrong_fr)
    performance.wrong_fr = wrong_fr;
end
if ~isempty(ts_wrong)
    performance.wrong_times = ts_wrong;
end
if ~isempty(wrongTriades)
    performance.wrongTriades = wrongTriades;
end
if ~isempty(frames_stem)
    performance.frames_stem = frames_stem;
end
if ~isempty(frames_left)
    performance.frames_left = frames_left;
end
if ~isempty(frames_right)
    performance.frames_right = frames_right;
end
if exist('ts_stemArm','var')
    performance.times_stem = ts_stemArm;
end
if exist('ts_leftArm','var')
    performance.times_left = ts_leftArm;
end
if exist('ts_rightArm','var')
    performance.times_right = ts_rightArm;
end
if ~isempty(ballistic)
    performance.ballistic = ballistic;
end
if ~isempty(doubted)
    performance.doubted = doubted;
end
if ~isempty(ballistic_perc)
    performance.ballistic_perc = ballistic_perc;
end
if ~isempty(doubted_perc)
    performance.doubted_perc = doubted_perc;
end
if ~isempty(times)
    performance.times = times;
end
if ~isempty(sampleVSchoice)
    performance.sampleVSchoice = sampleVSchoice;
end
if ~isempty(right)
    performance.right = right;
end
if ~isempty(wrong)
	performance.wrong = wrong;
end

%% save performance 
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
if saveMat
    save([basepath filesep sessionInfo.FileName '.Performance.Behavior.mat'],'performance');
end

end
