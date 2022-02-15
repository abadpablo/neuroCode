function bz_plotPlaceFields1D(varargin)

% Modification of plot_placeFields that takes into account if there is
% linear Track ( 2 conditions )



p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'firingMaps',[],@isstruct);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'tracking',[],@isstruct);
addParameter(p,'cell_metrics',[],@isstruct);
addParameter(p,'spikeTrain',[],@isstruct);

parse(p,varargin{:})

basepath = p.Results.basepath;
firingMaps = p.Results.firingMaps;
spikes = p.Results.spikes;
tracking = p.Results.tracking;
cell_metrics = p.Results.cell_metrics;
spikeTrain = p.Results.spikeTrain;

% Load sessionInfo and session
sessionInfo = bz_getSessionInfo();
session = loadSession;

mkdir('placeFields1D')

for i=1:length(tracking.apparatus)
    if strcmpi(tracking.apparatus{i}.name,'Linear Track  N-S')
        linearPos = i;
    end
end
conditions = [linearPos linearPos+1];


%% Need to include only spikeTrain analysis for folders with tracking
for i=1:length(tracking.folders)
    folder = tracking.folders{i};   
    if strcmpi(tracking.apparatus{i}.name,'Linear Track  N-S')
        spkTrain{i} = spikeTrain.(folder);
        spkTrain{i+1} = spikeTrain.(folder);
    end
end

% for i=1:length(tracking.folders)
%     timestamps{i} = tracking.events.subSessions;
% end

timestamps(linearPos,:) = tracking.events.subSessions(linearPos,:);

for i=1:length(firingMaps.rateMaps)
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        subplot(2,4,1)
        plot(spikes.filtWaveform{i})
        title(['Cell:' , num2str(i)])

        subplot(2,4,2)
        %ACG
        area(cell_metrics.acg.narrow(:,i),'LineStyle','none')
        title(['Cell type: ', num2str(cell_metrics.putativeCellType{i})])


        subplot(2,4,3)
        % ISI
        area(cell_metrics.isi.log10(:,i),'LineStyle','none')
        plot(cell_metrics.isi.log10(:,i))
        title(['Burst Index: ', num2str(spkTrain{linearPos}{i}.BurstIndex)])

        subplot(2,4,4)
        % Trilateralization
        plot(cell_metrics.general.chanCoords.x,cell_metrics.general.chanCoords.y,'.k'), hold on
        plot(cell_metrics.trilat_x(i),cell_metrics.trilat_y(i),'ob'), xlabel('x position (µm)'), ylabel('y position (µm)')
        title([ ' Shank:' ,num2str(spikes.shankID(i)),])

        subplot(2,4,5)
        [n,bin] = histc(spikes.times{i}(spikes.times{i} > timestamps(1,1) & spikes.times{i} < timestamps(1,2)),tracking.timestamps(tracking.events.subSessionsMask == linearPos));
        plot(tracking.position.x(tracking.events.subSessionsMask == linearPos),tracking.position.y(tracking.events.subSessionsMask == linearPos),'Color',[0.5 0.5 0.5])
        hold on
        view(0,-90)
        xToPlot = tracking.position.x(tracking.events.subSessionsMask == linearPos);
        yToPlot = tracking.position.y(tracking.events.subSessionsMask == linearPos);
        plot(xToPlot(bin(bin>0)),yToPlot(bin(bin>0)),'.','MarkerEdgeColor','r','MarkerSize',15);
        hold off
        set(gca,'DataAspectRatio',[1 1 1]);
        title(['Mean Firing Rate: ' , num2str(firingMaps.stats{i}{linearPos}.meanFr)])

        subplot(2,4,6)
        % Occupancy
        for j = 1:length(conditions)
            occupancy(j,:) = firingMaps.occupancy{i}{j};
        end
        if strcmpi(firingMaps.params.analysis,'tint')
            imagesc(occupancy);
            colormap(jet(15)),colorbar, shading flat
            title('occupancy')
                set(gca,'DataAspectRatio',[1 1 1]);
        elseif strcmpi(firingMaps.params.analysis,'buzcode')
            imagesc(occupancy);
            colormap(jet(15)),colorbar, shading flat
            view(0,-90)
            title('occupancy')
                set(gca,'DataAspectRatio',[1 1 1]);
        end


        subplot(2,4,7)
        %count
        for j = 1:length(conditions)
            countMaps(j,:) = firingMaps.countMaps{i}{j};
        end
        if strcmpi(firingMaps.params.analysis,'tint')
            imagesc(countMaps);
            colormap(jet(15)), colorbar, shading flat
            title('count')
            set(gca,'DataAspectRatio',[1 1 1]); 
        elseif strcpi(firingMaps.params.analysis,'buzcode')
            imagesc(countMaps);
            colormap(jet(15)), colorbar, shading flat
            title('count')
            view(0,-90)
            set(gca,'DataAspectRatio',[1 1 1]);
        end

        subplot(2,4,8)
        % rateMap
        for j = 1:length(conditions)
            rateMaps(j,:) = firingMaps.rateMaps{i}{j};
        end
        if strcmpi(firingMaps.params.analysis,'tint')
            imagesc(rateMaps);
            colormap(jet(15)), colorbar, shading flat
            title('ratemap')
            set(gca,'DataAspectRatio',[1 1 1]); 
        elseif strcmpi(firingMaps.params.analysis,'buzcode')
            imagesc(rateMaps);
            colormap(jet(15)), colorbar, shading flat
            title('ratemap')
            view(0,-90)
            set(gca,'DataAspectRatio',[1 1 1]); 
        end

        saveas(gcf,['placeFields1D\placeFields_Cell',num2str(i),'folder',tracking.folders{linearPos},'.png']);
end
    
% Now let's plot all the units together for both directions

figure,
set(gcf,'Position',get(0,'ScreenSize'))
for i = 1:length(conditions)
    subplot(length(conditions),1,i)
    for j = 1:length(firingMaps.rateMaps)
        rateMapsAll{i}(j,:) = firingMaps.rateMaps{j}{i};
    end
    
    imagesc(rateMapsAll{i})
    colormap(jet(15)), colorbar, shading flat 
end
saveas(gcf,['placeFields1D\placeFields_allCells.png'])

% Let's try to order firing Fields in 1D
figure,
set(gcf,'Position',get(0,'ScreenSize'))
    
for i = 1:length(rateMapsAll)
    subplot(length(conditions),1,i)
    rmapall = rateMapsAll{i};
    rmapall_order{i} = zeros(size(rmapall));
    for j = 1:size(rmapall,1)
        [a(j),b(j)] = find(rmapall == max(rmapall(j,:)))
    end
    [b_order,order] = sort(b);
    for k = 1:length(order)
        rmapall_order{i}(k,:) = rmapall(order(k),:);
    end
    imagesc(rmapall_order{i})
    colormap(jet(15)),colorbar,shading flat
end
saveas(gcf,['placeFields1D\placeFields_allCellsOrdered.png'])

% Now let's plot it normalized by the maximum firing rate per neuron

figure,
set(gcf,'Position',get(0,'ScreenSize'))
for i = 1:length(rmapall_order)
    subplot(length(rmapall_order),1,i)
    for j = 1:size(rmapall_order{i},1)
        rmapall_order{i}(j,:) = rmapall_order{i}(j,:) / max(rmapall_order{i}(j,:));
    end
    imagesc(rmapall_order{i})
    colormap(jet(15)), colormap(jet(15)), colorbar
end
saveas(gcf,['placeFields1D\placeFields_allCellsOrderedNormalized.png'])

close all

end






