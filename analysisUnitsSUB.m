%% ANALYSIS UNITS SUB PROJECT
% HPS22 040621 Linear Track + OF
bpath = 'C:\Projects\SUB\HPS22\HPS22_040621_sess25';
% bz_PreprocessSession('basepath',bpath);
% Ripple Channel 3, SharpWave Channel 13
% bz_analyzeSession('basepath',bpath,'forceReloadRipples',false,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','selectedRippleChannel',3,'selectedSWChannel',13);
cd(bpath)
sessionInfo = bz_getSessionInfo(bpath);
sessionInfo.AnatGrps(1).region = 'SUB';
sessionInfo.AnatGrps(2).region = 'SUB';
sessionInfo.AnatGrps(3).region = 'SUB';
sessionInfo.AnatGrps(4).region = 'SUB';
sessionInfo.AnatGrps(5).region = 'CA1';
sessionInfo.AnatGrps(6).region = 'CA1';
sessionInfo.AnatGrps(7).region = 'CA1';
sessionInfo.AnatGrps(8).region = 'CA1';
save([sessionInfo.FileName,'.sessionInfo.mat'],'sessionInfo');
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','selectedRippleChannel',3,'selectedSWChannel',13,'exclude',{'digitalPulses','performance','distanceByEpochs','lfp_analysis'});

% HPS22 210521 Linear Track + OF
bpath = 'C:\Projects\SUB\HPS22\HPS22_210521_sess17';
bz_PreprocessSession('basepath',bpath);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% HPS23 040621 Linear Track + OF
bpath = 'C:\Projects\SUB\HPS23\HPS23_040621_sess8';
% bz_PreprocessSession('basepath',bpath);
% Compute Ripples
% bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});
% RippleChannel: 47 SHarWaveChannel:31
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','selectedRippleChannel',47,'selectedSWChannel',31,'exclude',{'digitalPulses','performance','distanceByEpochs','lfp_analysis'});

% VB20R3
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_131221_sess1';
bz_PreprocessSession('basepath',bpath);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

%% VB20R3 151221
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_151221_sess3';
bz_PreprocessSession('basepath',bpath);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% medianSubstraction
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_151221_sess4';
bz_PreprocessSession('basepath',bpath,'medianSubstr',true);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});


% VB20R3 181221
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_181221_sess4';
% bz_PreprocessSession('basepath',bpath);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% VB20R3 181221 referenced
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_181221_sess5';
% bz_PreprocessSession('basepath',bpath,'refCh',11);
bz_PreprocessSession('basepath',bpath);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% medianSubstr
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_181221_sess5';
bz_PreprocessSession('basepath',bpath,'medianSubstr',true);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% VB20R3 191221 Linear Track + OF 
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_191221_sess5';
bz_PreprocessSession('basepath',bpath);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% VB20R3 191221 referenced Linear Track + OF referenced
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_191221_sess6';
% bz_PreprocessSession('basepath',bpath,'refCh',11);
bz_PreprocessSession('basepath',bpath);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% medianSubstr
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_191221_sess6';
bz_PreprocessSession('basepath',bpath,'medianSubstr',true);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% VB20R3 211221 Linear Track + OF 
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_211221_sess6';
bz_PreprocessSession('basepath',bpath);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% VB20R3 211221 Linear Track + OF referenced
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_211221_sess7';
% bz_PreprocessSession('basepath',bpath,'refCh',11);
bz_PreprocessSession('basepath',bpath);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% VB20R3 211221 Linear Track + OF medianSubstraction
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_211221_sess8';
% bz_PreprocessSession('basepath',bpath,'refCh',11);
bz_PreprocessSession('basepath',bpath,'medianSubstr',true);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});


% VB20R3 221221 Linear Track + OF
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_221221_sess7';
% bz_PreprocessSession('basepath',bpath,'refCh',11);
bz_PreprocessSession('basepath',bpath);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% VB20R3 221221 Linear Track + OF medianSubstr
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_221221_sess8';
% bz_PreprocessSession('basepath',bpath,'refCh',11);
bz_PreprocessSession('basepath',bpath,'medianSubstr',true);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% VB20R3 231221 OF
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_231221_sess8';
% bz_PreprocessSession('basepath',bpath,'refCh',11);
bz_PreprocessSession('basepath',bpath);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% VB20R3 231221 OF medianSubstr
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_231221_sess9';
% bz_PreprocessSession('basepath',bpath,'refCh',11);
bz_PreprocessSession('basepath',bpath,'medianSubstr',true);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});


%%
% prueba refCh
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_151221_sess5';
% bz_PreprocessSession('basepath',bpath,'refCh',16);
bz_PreprocessSession('basepath',bpath);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

bpath = 'C:\Projects\SUB\VB20R3\VB20R3_181221_sess6';
bz_PreprocessSession('basepath',bpath,'refCh',16);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

bpath = 'C:\Projects\SUB\VB20R3\VB20R3_191221_sess7';
bz_PreprocessSession('basepath',bpath,'refCh',16);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

bpath = 'C:\Projects\SUB\VB20R3\VB20R3_211221_sess9';
% bz_PreprocessSession('basepath',bpath,'refCh',16);
bz_PreprocessSession('basepath',bpath);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

bpath = 'C:\Projects\SUB\VB20R3\VB20R3_221221_sess9';
bz_PreprocessSession('basepath',bpath,'refCh',16);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

bpath = 'C:\Projects\SUB\VB20R3\VB20R3_231221_sess10';
bz_PreprocessSession('basepath',bpath,'refCh',16);
% Compute Ripples
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});


%% VB20R3_151221_sess4 prueba withouth PHY
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_151221_sess4';
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','performance','distancebyEpochs','lfp_analysis','newProtocol'});


bpath = 'C:\Projects\SUB\VB20R3\VB20R3_151221_sess5';
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','unitsSUB','pathExcel','C:\Projects\SUB','exclude',{'digitalPulses','powerSpectrumProfile','performance','distancebyEpochs','lfp_analysis','newProtocol'});


% VB20R3_291221 OR
bpath = 'C:\Projects\SUB\VB20R3\VB20R3_291221_sess9';
bz_PreprocessSession('basepath',bpath,'refCh',16);
