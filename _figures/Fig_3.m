%% Plot boxplots of general grouped results
%% Flow type, Seq Period, Seq Length

[d_upperPath, d_folder, ~] = fileparts(pwd);
if d_folder(2)=='_', d_folderTS=d_folder(5:15);
else d_folderTS=d_folder(1:11); end

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'), '-regexp', '^(?!r_flowSim)...');

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

colours=pmkmp(6,'CubicL');

% Select grouping dimensions
groupDims=[9 2 1];

% What's included in the summary
d_reqSeqLength=1:length(unique(r_seqLength));
d_reqSeqPeriod=1:length(unique(r_seqPeriod));
d_reqNZones=1:length(unique(r_nZones));
d_reqSolve=[1];
d_reqImp=[1];
d_reqNSeqA=[1];
d_reqTSeqA=[1];
d_reqConc=[4];
d_reqFlowType=[1 2 3];

d_dir{1}='P6';
d_figAve=1;

d_out=1

figOut=figure;
set(figOut,'Units','centimeters');
set(figOut,'Position', [5, 5, 17, 7]);
set(figOut,'Units','pixels');
set(gcf,'Renderer','Painters');

P7_Grouping; 
    
for d_i=1:3        
    d_summary=out_summary(find(out_summary(:,1)==d_i),:);
    flowStats=d_summary(:,11:15)'*100;
    nTotal=sum(d_summary(:,6));
    grouping1=d_summary(:,3)';
    grouping2=d_summary(:,2)';
    nGrouping1=length(unique(d_summary(:,3)));
    nGrouping2=length(unique(d_summary(:,2)));

    handAxesM(d_i)=axes;
    data=rand(20,size(flowStats,2));
    h=boxplotDV(data,{grouping2,grouping1},'factorgap',[3 0.1],'colors',repmat(flipud(colours(2:2+(nGrouping1-1),:)),nGrouping2,1),'factorseparator',[1], 'widths', 1, 'boxstyle','outline');
    set(h(1,:),{'Ydata'},num2cell(flowStats(end-1:end,:),1)')
    set(h(2,:),{'Ydata'},num2cell(flowStats(2:-1:1,:),1)')
    set(h(3,:),{'Ydata'},num2cell(flowStats([end end],:),1)')
    set(h(4,:),{'Ydata'},num2cell(flowStats([1 1],:),1)')
    set(h(5,:),{'Ydata'},num2cell(flowStats([2 4 4 2 2],:),1)')
    set(h(6,:),{'Ydata'},num2cell(flowStats([3 3],:),1)')
    set(h(7,:),{'Visible'},{'off'})
    ylim([-100 100]);
    set(handAxesM(d_i),'color','none','Position', [sizeXMOff+rem(d_i-1,3)*(sizeXM+sizeXPad) 1-(sizeYMOff+(ceil(d_i/3))*(sizeYM+sizeYPad)) sizeXM sizeYM]);

    h = findobj(gca,'Tag','Box');
    for x=1:length(h)
        boxShade=patch(get(h(x),'XData'),get(h(x),'YData'), colours(rem(x-1,nGrouping1)+2,:), 'linestyle', 'none');
        uistack(boxShade,'bottom');
    end
    delete(h);
    set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-','LineWidth',1)
    set(findobj(gcf,'Tag','Upper Adjacent Value'),'LineStyle','-','LineWidth',1)
    set(findobj(gcf,'Tag','Lower Adjacent Value'),'LineStyle','-','LineWidth',1)
    set(findobj(gcf,'Tag','Median'),'Color','Red','LineStyle','-','LineWidth',0.5)

    graphX=max(xlim);
    graphXint=graphX/6;

    set(gca,'xTick', graphXint/2:graphXint:graphX, 'XTickLabel',{'2' '3' '4' '6' '12' '24'})
    legendText={'15','31','63','127','255'};
    legendText=legendText(1:nGrouping1);

    switch d_i
        case 1
            title('Total Dwelling Flow', 'FontSize', 10);
            hYLabel = ylabel('Flow Weighted Error (%)');
            set(hYLabel, 'Units', 'Normalized', 'Position', [-0.12, 0.5, 0]);
            hXLabel = xlabel('PRBS Period (hours)');
        case 2
            title('Total Zonal Flows','FontSize', 10);
            set(gca,'YTickLabel',[]);
            hXLabel = xlabel('PRBS Period (hours)');
        case 3
            title('Individual Flows','FontSize', 10);
            set(gca,'YTickLabel',[]);
            hXLabel = xlabel('PRBS Period (hours)');
            legendflex(legendText,'xscale', 0.5, 'title', {'PRBS Length (bits)'}, 'padding', [10 10 10], 'box', 'on', 'ncol', 2);
            % handAxesM(3),
    end
    text((max(xlim)-min(xlim))*0.8+min(xlim), (max(ylim)-min(ylim))*0.05+min(ylim), ['n = ' separatethousands(nTotal,',',0)], 'FontSize', 7, 'VerticalAlignment' , 'middle', 'HorizontalAlignment' , 'center', 'BackgroundColor', 'White');

    set(gca, ...
      'TickDir'     , 'out'         , ...
      'TickLength'  , [0 0]         , ...
      'XMinorTick'  , 'off'         , ...
      'YMinorTick'  , 'off'         , ...
      'YGrid'       , 'on'          , ...
      'YTick'       , -100:20:100   , ...
      'Layer'       , 'top');
end
    
fileSaveName=['Plot_' num2str(d_out) '.pdf'];
export_fig('filename', fileSaveName, '-nocrop');
close(gcf);

