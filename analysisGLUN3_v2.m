% IPO429 111021 VEHICLE
bpath = 'C:\Projects\GLUN3\IPO429\IPO429_111021_sess1';
bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[2 3 4 5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO429 131021 KETAMINE
bpath = 'C:\Projects\GLUN3\IPO429\IPO429_131021_sess2';
bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[2 3 4 5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO429 151021 MK801
bpath = 'C:\Projects\GLUN3\IPO429\IPO429_151021_sess3';
bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[2 3 4 5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO430 111021 VEHICLE
bpath = 'C:\Projects\GLUN3\IPO430\IPO430_111021_sess1';
bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[2 3 4 5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','plotLinearTrack','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','newProtocol','excel','plotPlaceFields'});

% IPO430 131021 KETAMINE
bpath = 'C:\Projects\GLUN3\IPO430\IPO430_131021_sess2';
bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[2 3 4 5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO430 151021 MK801
bpath = 'C:\Projects\GLUN3\IPO430\IPO430_151021_sess3';
bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[2 3 4 5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO447 181121 VEHICLE
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_181121_sess1';
% bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO447 191121 MK801
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_191121_sess2';
% bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO447 221121 KETAMINE
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_221121_sess3';
% bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO447 241121 KETAMINE AGAIN
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_241121_sess4';
% bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO447 091221 VEHICLE NEW PROTOCOL
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_091221_sess5';
bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO447 101221 MK801 NEW PROTOCOL ( CREO QUE NO HA FUNCIONADO)
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_101221_sess6';
bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO447 101221 MK801 NEW PROTOCOL ( VARIOS REGISTROS Y DOS INYECCIONES)
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_121221_sess7';
bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO447 141221 KETAMINE NEW PROTOCOL
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_141221_sess8';
bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO447 181221 MK801 NEW PROTOCOL ( good experiment)
bpath = 'C:\Projects\GLUN3\IPO447\IPO447_161221_sess9';
bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO149 231021 MK801
bpath = 'C:\Projects\GLUN3\IPO149\IPO149_231021_sess1';
% bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO149 251021 KETAMINE
bpath = 'C:\Projects\GLUN3\IPO149\IPO149_251021_sess2';
% bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO149 281021 VEHICLE
bpath = 'C:\Projects\GLUN3\IPO149\IPO149_281021_sess3';
% bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO150 261021 MK801
bpath = 'C:\Projects\GLUN3\IPO150\IPO150_261021_sess1';
% bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO150 291021 KETAMINE
bpath = 'C:\Projects\GLUN3\IPO150\IPO150_291021_sess3';
% bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO150 091121 VEHICLE
bpath = 'C:\Projects\GLUN3\IPO150\IPO150_091121_sess2';
% bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO150 301121 KETAMINE NEW PROTOCOL
bpath = 'C:\Projects\GLUN3\IPO150\IPO150_301121_sess4';
bz_PreprocessSession('basepath',bpath,'getDigital',false);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO150 021221 MK801 NEW PROTOCOL
bpath = 'C:\Projects\GLUN3\IPO150\IPO150_021221_sess5';
bz_PreprocessSession('basepath',bpath,'getDigital',false);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});

% IPO150 121221 VEHICLE NEW PROTOCOL
bpath = 'C:\Projects\GLUN3\IPO150\IPO150_121221_sess6';
bz_PreprocessSession('basepath',bpath,'getDigital',false);
bz_analyzeSession('basepath',bpath,'excludeShanks',[5 6 7 8],'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','distanceByEpochs','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});
