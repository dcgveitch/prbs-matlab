%% Overlapping nature of cross correlation impulses
%% Use with data from '3_1_140922T1601_Impulses'

clear
[d_upperPath, d_folder, ~] = fileparts(pwd);
if d_folder(2)=='_', d_folderTS=d_folder(5:15);
else d_folderTS=d_folder(1:11); end

cd Results;
load(strcat(d_folderTS(1:11), '__outP1.mat'));
load(strcat(d_folderTS(1:11), '__outP2.mat'));
load(strcat(d_folderTS(1:11), '__outP3.mat'));
load(strcat(d_folderTS(1:11), '_setup.mat'), '-regexp', '^(?!r_flowSim)...');

seqdT=30/3600;
seqLength=12;
seqRepeat=2;
nZones=2;

input=downsample(out_impulseSim{1}(:,[1 3]),30);

for i=1:ceil(seqLength/seqRepeat)
    for j=1:nZones
        count=(i-1)*nZones+j;
        impulseT(((seqRepeat*((i-1)+(j-1)/nZones))/seqdT)+1:((seqRepeat*((i-1)+(j-1)/nZones))/seqdT)+length(input),count)=input(:,j);
    end
end

impulseT=impulseT(1:seqLength/seqdT,:);

nX=1;
nY=2;

sizeXMOff=0.07;
sizeYMOff=0.07;
sizeXTotal=0.9;
sizeYTotal=0.78;
sizeXPad=0.01*17/14;
sizeYPad=0.02;

sizeXM=(sizeXTotal-(nX-1)*sizeXPad)/nX;
sizeYM=(sizeYTotal-(nY-1)*sizeYPad)/nY;

colours=pmkmp(6,'CubicL');

figOut=figure;
set(figOut,'Units','centimeters');
set(figOut,'Position', [5, 5, 9, 7]);
set(figOut,'Units','pixels');

set(gcf,'Renderer','Painters');    

for d_i=1:2
    handAxesM(d_i)=axes('color','none');
    if (d_i==2)
        areaP=area(impulseT,'linestyle','none');
        coloursB=brighten(colours,0.75);
        for d_j=1:size(impulseT,2)/nZones
            for d_k=1:nZones
                count=(d_j-1)*nZones+d_k;
                set(areaP(count),'facecolor',coloursB(d_k+4,:));
            end
        end
        hold on
        x=3:1080;
        y=sum(impulseT,2);
        lineP=plot(x,y(3:1080),'--k');
        xlim([0 960]);
        ylim([0 2000]);
        ylab=ylabel({'Cross-Correlation' 'Output'});
        set(ylab,'Position',get(ylab,'Position')+[105 0 0]);
        xlabel('Monitoring Time (hours)');
        text(max(xlim)*0.01, max(ylim)*0.75,{'Equal to superposition'; 'of impulse responses'},'FontSize', 7, 'VerticalAlignment', 'Bottom');
        set(gca,...
        'TickLength'    , [0.01 0] , ...
        'XTickLabel'    , {'0','2','4','6','8','10'} , ...    
        'XTick'         , [0:240:1200] , ...
        'YTickLabel'    , [] , ...
        'YTick'         , [0:500:2000] , ...
        'GridColor'     , [0.5 0.5 0.5] , ...
        'GridAlpha'     , 0.85);
    elseif (d_i==1)
        outlineP=plot(impulseT);
        xlim([0 960]);
        ylim([0 2000]);
        for d_j=1:size(impulseT,2)/nZones
            for d_k=1:nZones
                count=(d_j-1)*nZones+d_k;
                set(outlineP(count),'color',colours(d_k+4,:));
            end
        end
        line([240 240],[1050 2000],'Color','k','LineWidth', 0.2);
        line([480 480],[1050 2000],'Color','k','LineWidth', 0.2);
        line([720 720],[1050 2000],'Color','k','LineWidth', 0.2);
        line([960 960],[1050 2000],'Color','k','LineWidth', 0.2);
        ylab=ylabel({'Z1 Impulse' 'Responses'});
        set(ylab,'Position',get(ylab,'Position')+[105 0 0]);
        set(gca,...
        'TickLength'    , [0.01 0] , ...
        'XTickLabel'    , [] , ...
        'XTick'         , [0:120:1200] , ...
        'YTickLabel'    , [] , ...
        'YTick'         , [0:500:2000] , ...
        'GridColor'     , [0.5 0.5 0.5] , ...
        'GridAlpha'     , 0.85);
    end
    grid on;
    set(handAxesM(d_i),'Position', [sizeXMOff+rem(d_i-1,nX)*(sizeXM+sizeXPad) 1-(sizeYMOff+ceil(d_i/nX)*(sizeYM+sizeYPad)) sizeXM sizeYM], 'color', 'none');
    set(gca, 'layer', 'top');
end

d_legend{1}=['Z1'];
d_legend{2}=['Z2'];

legendflex(handAxesM(1),d_legend,'xscale', 0.5, 'title', {'Input Location'}, 'padding', [10 10 10], 'box', 'on', 'ncol', 2);

[Xf2a, Yf2a] = ds2nfu(handAxesM(1), 240, 1100);
[Xf2b, Yf2b] = ds2nfu(handAxesM(1), 340, 1100);
[Xf3a, Yf3a] = ds2nfu(handAxesM(1), 480, 1100);
[Xf3b, Yf3b] = ds2nfu(handAxesM(1), 580, 1100);
[Xf4a, Yf4a] = ds2nfu(handAxesM(1), 720, 1100);
[Xf4b, Yf4b] = ds2nfu(handAxesM(1), 820, 1100);
annotation('arrow', [Xf2a,Xf2b], [Yf2a,Yf2b],'Color','k', 'LineWidth', 0.2);
annotation('arrow', [Xf3a,Xf3b], [Yf3a,Yf3b],'Color','k', 'LineWidth', 0.2);
annotation('arrow', [Xf4a,Xf4b], [Yf4a,Yf4b],'Color','k', 'LineWidth', 0.2);

POSf = ds2nfu(handAxesM(1), [240-10 1650 200 100]');
str={'PRBS', 'Repeats'};
annotation('textbox', POSf',...
           'String', str,...
           'LineSTyle', 'none',...
           'FontSize', 7,...
           'HorizontalAlignment', 'left',...
           'VerticalAlignment', 'top');


fileSaveName=['2_Superposition.pdf'];
export_fig('filename', fileSaveName, '-nocrop');
close(gcf);

