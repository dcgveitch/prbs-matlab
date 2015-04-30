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

nX=3;
nY=1;

sizeXMOff=0.07;
sizeYMOff=0.07;
sizeXTotal=0.9;
sizeYTotal=0.78;
sizeXPad=0.01;
sizeYPad=0.01;

sizeXM=(sizeXTotal-(nX-1)*sizeXPad)/nX; 
sizeYM=(sizeYTotal-(nY-1)*sizeYPad)/nY;

fileSaveName='MixRatio.pdf';
figOut=figure;
set(figOut,'Units','centimeters');
set(figOut,'Position', [5, 5, 17, 7]);
set(figOut,'Units','pixels');

for d_i=1:3        
    handAxesM(d_i)=axes;
    h=gscatter(out_resultsCombined(:,5),out_resultsCombined(:,d_i)*100,out_resultsCombined(:,4))
    ylim([-5 5]);
    xlim([0 49]);
    set(handAxesM(d_i),'color','none','Position', [sizeXMOff+rem(d_i-1,3)*(sizeXM+sizeXPad) 1-(sizeYMOff+(ceil(d_i/3))*(sizeYM+sizeYPad)) sizeXM sizeYM]);
    legend OFF;

    legendText={'2','3','5','8'};  
    hXLabel = xlabel('Internal/External Flow Ratio');

    switch d_i
        case 1
            title('Total Dwelling Ventilation');
            hYLabel = ylabel('Flow Weighted Error (%)');
            set(hYLabel, 'Units', 'Normalized', 'Position', [-0.12, 0.5, 0]);
        case 2
            title('Zonal Infiltration');
            set(gca,'YTickLabel',[]);
        case 3
            title('Individual Airflows');
            set(gca,'YTickLabel',[]);
            legendflex(legendText,'xscale', 0.5, 'title', {'nZones'}, 'padding', [10 10 10], 'box', 'on', 'ncol', 2);
    end

    set(gca, ...
      'TickDir'     , 'in'         , ...
      'TickLength'  , [0.01 0]         , ...
      'XMinorTick'  , 'off'         , ...
      'YMinorTick'  , 'off'         , ...
      'YGrid'       , 'on'          , ...
      'GridColor'   , [0.5 0.5 0.5] , ...
      'GridAlpha'   , 1             , ...
      'YTick'       , -5:2.5:5);      
end

export_fig('filename', fileSaveName, '-nocrop');
close(gcf);

