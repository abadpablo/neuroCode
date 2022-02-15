%% ANALISIS GLUN3

% Update Exp Folder
updateExpFolder('\\DISCOVERY_ONE\Sub\HPS23','F:\data\HPS23');
arrangeSessionFolder('C:\Projects\GLUN3\IPO447');
bpath= 'C:\Projects\GLUN3\IPO430';
createFiles('basepath',bpath);
%% PREPROCESSING

% HPS22 100621 MK801
bpath = 'F:\data\HPS22\HPS22_100621_sess26';
% bz_PreprocessSession('basepath',bpath);
% Compute Ripples
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
% Compute general analysis
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','mk801');

% HPS22 150621 VEHICLE
bpath = 'F:\data\HPS22\HPS22_150621_sess27';
% bz_PreprocessSession('basepath',bpath);
% Compute Ripples
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});
% Compute general analysis RippleChannel 30 SharpWaveChannel 12
bz_analyzeSession('basepath',bpath,'forceReloadRipples',true,'selectedRippleChannel',30,'selectedSWChannel',12,'analyzeSubSessions',true,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','vehicle');

% HPS22 160621 KETAMINE
bpath = 'F:\data\HPS22\HPS22_160621_sess28';
% bz_PreprocessSession('basepath',bpath);
% Compute Ripples
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});
% Compute general analysis
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','ketamine');

% HPS23 090621 MK801
bpath = 'F:\data\HPS23\HPS23_090621_sess9';
% bz_PreprocessSession('basepath',bpath);
% Compute Ripples
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});
% Compute general analysis
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN32_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','mk801');

% HPS23 110621 VEHICLE
bpath = 'F:\data\HPS23\HPS23_110621_sess10';
% bz_PreprocessSession('basepath',bpath);
% Compute Ripples
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});
% Compute general analysis
bz_analyzeSession('basepath',bpath,'forceReloadRipples',true,'selectedRippleChannel',46,'selectedSWChannel',30,'analyzeSubSessions',true,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','vehicle');

% HPS23 160621 KETAMINE
bpath = 'F:\data\HPS23\HPS23_160621_sess11';
% bz_PreprocessSession('basepath',bpath);
% Compute Ripples
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});
% Compute general analysis
bz_analyzeSession('basepath',bpath,'forceReloadRipples',true,'analyzeSubSessions',true,'selectedRippleChannel',46,'selectedSWChannel',30,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','ketamine');

% HPS24 230621 KETAMINE
bpath = 'F:\data\HPS24\HPS24_230621_sess1';
% bz_PreprocessSession('basepath',bpath);
% Compute Ripples
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});
% Compute general analysis
bz_analyzeSession('basepath',bpath,'forceReloadRipples',true,'selectedRippleChannel',39,'selectedSWChannel',24,'analyzeSubSessions',true,'selectedRippleChannel',36,'selectedSWChannel',26,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','ketamine');

% HPS24 280621 VEHICLE
bpath = 'F:\data\HPS24\HPS24_280621_sess3';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'forceReloadRipples',true,'selectedRippleChannel',39,'selectedSWChannel',24,'analyzeSubSessions',true,'selectedRippleChannel',36,'selectedSWChannel',26,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','spikeTrain','distanceByEpochs','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','newProtocol','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','vehicle');

% HPS24 290621 MK801
bpath = 'F:\data\HPS24\HPS24_290621_sess4';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'forceReloadRipples',true,'selectedRippleChannel',39,'selectedSWChannel',24,'analyzeSubSessions',true,'selectedRippleChannel',36,'selectedSWChannel',26,'nameExcel','GLUN3_161221','exclude',{'thetaMOdulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','mk801');

% HPS25 240621 KETAMINE
bpath = 'F:\data\HPS25\HPS25_240621_sess2';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'forceReloadRipples',true,'selectedRippleChannel',46,'selectedSWChannel',16,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','ketamine');

% HPS25 280621 VEHICLE 
bpath = 'F:\data\HPS25\HPS25_280621_sess3';
% bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'forceReloadRipples',true,'selectedRippleChannel',60,'selectedSWChannel',3,'analyzeSubSessions',true,'nameExcel','GLUN3_030122','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs','newProtocol'});
bz_BaselineVsDrug('basepath',bpath,'drug','vehicle');

% HPS25 050721 MK801
bpath = 'F:\data\HPS25\HPS25_050721_sess4';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'forceReloadRipples',true,'selectedRippleChannel',60,'selectedSWChannel',3,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','mk801');

% IPO429 111021 VEHICLE
bpath = 'C:\Projects\GLUN3\IPO429\IPO429_111021_sess1';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[2 3 4 5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',15,'selectedSWChannel',17,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','vehicle');

% IPO429 131021 KETAMINE
bpath = 'C:\Projects\GLUN3\IPO429\IPO429_131021_sess2';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[2 3 4 5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',15,'selectedSWChannel',17,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','ketamine');

% IPO429 151021 MK801
bpath = 'C:\Projects\GLUN3\IPO429\IPO429_151021_sess3';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[2 3 4 5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',15,'selectedSWChannel',17,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','performance','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','mk801','yMaze',false);

% IPO430 111021 VEHICLE
bpath = 'C:\Projects\GLUN3\IPO430\IPO430_111021_sess1';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[2 3 4 5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',15,'selectedSWChannel',17,'nameExcel','GLUN3_030122','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs','newProtocol'});
bz_BaselineVsDrug('basepath',bpath,'drug','vehicle');

% IPO430 131021 KETAMINE
bpath = 'C:\Projects\GLUN3\IPO430\IPO430_131021_sess2';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[2 3 4 5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',15,'selectedSWChannel',17,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','ketamine','yMaze',false);

% IPO430 151021 MK801
bpath = 'C:\Projects\GLUN3\IPO430\IPO430_151021_sess3';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[2 3 4 5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',15,'selectedSWChannel',17,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','mk801','yMaze',false);

% IPO447 181121 VEHICLE CA1-PFC 
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_181121_sess1';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',11,'selectedSWChannel',3,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','vehicle','yMaze',false);

% IPO447 191121 MK801 CA1-PFC
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_191121_sess2';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',11,'selectedSWChannel',3,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','mk801','yMaze',false);

% IPO447 221121 KETAMINE CA1-PFC
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_221121_sess3';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',11,'selectedSWChannel',3,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','ketamine','yMaze',false);

% IPO447 241121 KETAMINE AGAIN CA1-PFC
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_241121_sess4';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',11,'selectedSWChannel',3,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','ketamine','yMaze',false);

% IPO447 091221 VEHICLE NEW PROTOCOL CA1-PFC
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_091221_sess5';
% bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',11,'selectedSWChannel',3,'nameExcel','GLUN3_NewProtocol','exclude',{'thetaModulation','placeCells','performance','plotLinearTrack','plotPlaceFields','digitalPulses','lfp_analysis','CellExplorer'});

% IPO447 141221 KETAMINE NEW PROTOCOL CA1-PFC
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_141221_sess8';
% bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',11,'selectedSWChannel',3,'nameExcel','GLUN3_NewProtocol','exclude',{'thetaModulation','placeCells','performance','plotLinearTrack','plotPlaceFields','digitalPulses','lfp_analysis','CellExplorer'});

% IPO447 161221 MK801 NEW PROTOCOL ( good experiment)
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_161221_sess9';
% bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',11,'selectedSWChannel',3,'nameExcel','GLUN3_NewProtocol','exclude',{'thetaModulation','placeCells','performance','plotLinearTrack','plotPlaceFields','digitalPulses','lfp_analysis','CellExplorer'});

% IPO149 231021 MK801 
bpath = 'C:\Projects\GLUN3\IPO149\IPO149_231021_sess1';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',1,'selectedSWChannel',0,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','mk801','yMaze',false);

% IPO149 251021 KETAMINE 
bpath = 'C:\Projects\GLUN3\IPO149\IPO149_251021_sess2';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',1,'selectedSWChannel',0,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','ketamine','yMaze',false);

% IPO149 281021 VEHICLE
bpath = 'C:\Projects\GLUN3\IPO149\IPO149_281021_sess3';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',1,'selectedSWChannel',0,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','vehicle','yMaze',false);

% IPO150 261021 MK801
bpath = 'C:\Projects\GLUN3\IPO150\IPO150_261021_sess1';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',1,'selectedSWChannel',0,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','mk801','yMaze',false);

% IPO150 291021 KETAMINE
bpath = 'C:\Projects\GLUN3\IPO150\IPO150_291021_sess3';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',1,'selectedSWChannel',0,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','ketamine','yMaze',false);

% IPO150 091121 VEHICLE
bpath = 'C:\Projects\GLUN3\IPO150\IPO150_091121_sess2';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'forceReloadRipples',true,'selectedRippleChannel',1,'selectedSWChannel',0,'nameExcel','GLUN3_161221','exclude',{'thetaModulation','placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});
bz_BaselineVsDrug('basepath',bpath,'drug','vehicle','yMaze',false);

% IPO150 301121 KETAMINE NEW PROTOCOL
bpath = 'C:\Projects\GLUN3\IPO150\IPO150_301121_sess4';
% bz_PreprocessSession('basepath',bpath);
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','plotPlaceFields'});
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'placeCells','plotLinearTrack','plotPlaceFields','digitalPulses','CellExplorer','distanceByEpochs'});

% IPO150 021221 MK801 NEW PROTOCOL




% IPO150 121221 VEHICLE NEW PROTOCOL



