function  [] = plotRasterAbad(spkdata,varargin)

figure,
set(gcf,'Position',get(0,'ScreenSize'))
hold on;

for iCell = 1:length(spkdata)
    spks = spkdata{iCell}';
    xspikes = repmat(spks,3,1);
    yspikes = nan(size(xspikes));
    
    if ~isempty(yspikes)
        yspikes(1,:) = iCell - 1;
        yspikes(2,:) = iCell;
    end
    
    plot(xspikes,yspikes,'Color','k')

end

