% Prueba DistanceByEpochs
bpath = 'C:\Projects\GLUN3\IPO150\IPO150_301121_sess4';
bz_analyzeSession('basepath',bpath,'analyzeSubSessions',true,'nameExcel','GLUN3','exclude',{'spikes','meanCrossCorr','digitalPulses','ripples','powerSpectrumProfile','thetaModulation','behaviour','placeCells','performance','spikeTrain','CellExplorer','lfp_analysis','excel','plotPlaceFields'});
