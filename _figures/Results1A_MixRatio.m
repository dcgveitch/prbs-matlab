%% Plot boxplots of general grouped results
%% Flow type, Seq Period, Seq Length

% Make boxplots of the non-subsmapled theoretical results
% Interested only in no-noise with differnet combinations of zones

clear
[d_upperPath, d_folder, ~] = fileparts(pwd);
if d_folder(2)=='_', d_folderTS=d_folder(5:15);
else d_folderTS=d_folder(1:11); end

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'), '-regexp', '^(?!r_flowSim)...');
cd P6;
load(strcat(d_folderTS(1:11), '__outP6Single.mat'));

nX=4;
nXpad=2;
nY=1;

figPosition=[5, 5, 16, 7];

sizeXMOff=0.05;
sizeYMOff=0.05;
sizeXTotal=0.94;
sizeYTotal=0.83;
sizeXPad=0.01;
sizeYPad=0.01;

sizeXM=(sizeXTotal-((nX-1)*sizeXPad)-(nX/nXpad-1)*sizeXPad)/nX; 
sizeYM=(sizeYTotal-(nY-1)*sizeYPad)/nY;

fileSaveName='MixRatio.pdf';
figOut=figure;
set(figOut,'Units','centimeters');
set(figOut,'Position', figPosition);
set(figOut,'Units','pixels');

colours=pmkmp(5,'IsoL');
sizeMult=[10 20 10 20];

for d_i=1:4        
    handAxesM(d_i)=axes;
    out_resultsCombined{d_i}=flipud(out_resultsCombined{d_i});
    h=scatter(out_resultsCombined{d_i}(:,4),out_resultsCombined{d_i}(:,2)*100,out_resultsCombined{d_i}(:,1)*sizeMult(d_i),'k','filled');
%     h=scatter(out_resultsCombined{d_i}(:,4),out_resultsCombined{d_i}(:,2)*100,out_resultsCombined{d_i}(:,1)*10,colours(out_resultsCombined{d_i}(:,3)+1,:),'filled');
    ylim([-5 5]);
    xlim([0 49]);
    set(handAxesM(d_i),'box','on','Position',[sizeXMOff+rem(d_i-1,nX)*(sizeXM+sizeXPad)+floor(rem(d_i-1,nX)/nXpad)*sizeXPad 1-(sizeYMOff+(ceil(d_i/nX))*(sizeYM+sizeYPad)) sizeXM sizeYM]);
    legend OFF;

    legendText={'2','3','5','8'};  
    hXLabel = xlabel('Int/Ext Flow Ratio');

    switch d_i
        case 1
            title('\rm{Total}');
            hYLabel = ylabel('Flow Error (%)');
            set(hYLabel, 'Units', 'Normalized', 'Position', [-0.12, 0.5, 0]);
            text((max(xlim)-min(xlim))*(1+((sizeXPad/2)/(sizeXTotal/nX)))+min(xlim), (max(ylim)-min(ylim))*1.0+min(ylim), '\bf{-External-}', 'FontSize', 9, 'VerticalAlignment' , 'bottom', 'HorizontalAlignment' , 'center');
        case 2
            title('\rm{Individual}');
            set(gca,'YTickLabel',[]);
        case 3
            title('\rm{Total}');
            set(gca,'YTickLabel',[]);
            text((max(xlim)-min(xlim))*(1+((sizeXPad/2)/(sizeXTotal/nX)))+min(xlim), (max(ylim)-min(ylim))*1.0+min(ylim), '\bf{-Internal-}', 'FontSize', 9, 'VerticalAlignment' , 'bottom', 'HorizontalAlignment' , 'center');
        case 4
            title('\rm{Individual}');
            set(gca,'YTickLabel',[]);
%             legendflex(legendText,'xscale', 0.5, 'title', {'nZones'}, 'padding', [10 10 10], 'box', 'on', 'ncol', 2);
    end

    set(gca, ...
      'TickDir'     , 'in'         , ...
      'TickLength'  , [0.01 0]         , ...
      'XMinorTick'  , 'off'         , ...
      'YMinorTick'  , 'off'         , ...
      'YGrid'       , 'on'          , ...
      'XGrid'       , 'on'          , ...
      'GridColor'   , [0.5 0.5 0.5] , ...
      'GridAlpha'   , 1             , ...
      'YTick'       , -5:2.5:5);      
end

export_fig('filename', fileSaveName, '-nocrop');
close(gcf);

