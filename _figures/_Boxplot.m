%% Base boxplot script for general grouped results

% Generate grouped 
P7_Grouping;

if (~isempty(out_summary))
    figOut=figure;
    set(figOut,'Units','centimeters');
    set(figOut,'Position', figPosition);
    set(figOut,'Units','pixels');
    
    % Make list of max combinations in figure
    groupingFull2=req_i{groupDims(2)};
    groupingFull3=req_i{groupDims(3)};
    groupingFull=[];
    for d_j=1:length(groupingFull2)
        for d_k=1:length(groupingFull3)
            groupingFull(end+1,:)=[groupingFull2(d_j) groupingFull3(d_k)];
        end
    end

    for d_i=req_i{groupDims(1)}
        d_summary=out_summary(find(out_summary(:,1)==d_i),:);
        groupingMiss=setdiff(groupingFull,d_summary(:,2:3),'rows');

        % Pad missing data
        if (~isempty(groupingMiss))
            for d_j=1:size(groupingMiss,1)
                d_summary(end+1,1)=d_i;
                d_summary(end,2:3)=groupingMiss(d_j,1:2);
                d_summary(end,4:5)=1;
                d_summary(end,6)=0;
                d_summary(end,7:end)=999;
            end
            d_summary=sortrows(d_summary,[2 3]);
        end

        flowStats=d_summary(:,11:15)'*100;
        nTotal=sum(d_summary(:,6));
        grouping1=d_summary(:,3)';
        grouping2=d_summary(:,2)';
        nGrouping1=length(unique(d_summary(:,3)));
        nGrouping2=length(unique(d_summary(:,2)));

        handAxesM(d_i)=axes;
        data=rand(20,size(flowStats,2));
        h=boxplot(data,{grouping2,grouping1},'factorgap',[3 0.1],'colors',repmat(flipud(colours(2:2+(nGrouping1-1),:)),nGrouping2,1), 'factorseparator', [1], 'widths', 1, 'boxstyle','outline');
        set(h(1,:),{'Ydata'},num2cell(flowStats(end-1:end,:),1)')
        set(h(2,:),{'Ydata'},num2cell(flowStats(2:-1:1,:),1)')
        set(h(3,:),{'Ydata'},num2cell(flowStats([end end],:),1)')
        set(h(4,:),{'Ydata'},num2cell(flowStats([1 1],:),1)')
        set(h(5,:),{'Ydata'},num2cell(flowStats([2 4 4 2 2],:),1)')
        set(h(6,:),{'Ydata'},num2cell(flowStats([3 3],:),1)')
        set(h(7,:),{'Visible'},{'off'})
        ylim([-100 100]);
        set(handAxesM(d_i),'color','none','Position', [sizeXMOff+rem(d_i-1,nX)*(sizeXM+sizeXPad)+floor(rem(d_i-1,nX)/nXpad)*sizeXPad 1-(sizeYMOff+(ceil(d_i/nX))*(sizeYM+sizeYPad)) sizeXM sizeYM]);

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
                title('\rm{Total}');
                hYLabel = ylabel('Flow Weighted Error (%)');
                set(hYLabel, 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0]);
                hXLabel = xlabel('PRBS Period (hours)');
                text((max(xlim)-min(xlim))*(1+((sizeXPad/2)/(sizeXTotal/nX)))+min(xlim), (max(ylim)-min(ylim))*1.0+min(ylim), '\bf{-External-}', 'FontSize', 9, 'VerticalAlignment' , 'bottom', 'HorizontalAlignment' , 'center');
            case 2
                title('\rm{Individual}');
                set(gca,'YTickLabel',[]);
                hXLabel = xlabel('PRBS Period (hours)');
            case 3
                title('\rm{Total}');
                set(gca,'YTickLabel',[]);
                hXLabel = xlabel('PRBS Period (hours)');
                text((max(xlim)-min(xlim))*(1+((sizeXPad/2)/(sizeXTotal/nX)))+min(xlim), (max(ylim)-min(ylim))*1.0+min(ylim), '\bf{-Internal-}', 'FontSize', 9, 'VerticalAlignment' , 'bottom', 'HorizontalAlignment' , 'center');
            case 4
                title('\rm{Individual}');
                set(gca,'YTickLabel',[]);
                hXLabel = xlabel('PRBS Period (hours)');
                legendflex(legendText,'xscale', 0.5, 'title', {'PRBS Length (bits)'}, 'padding', [10 10 10], 'box', 'on', 'ncol', 2);
        end
        text((max(xlim)-min(xlim))*0.96+min(xlim), (max(ylim)-min(ylim))*0.055+min(ylim), ['n = ' separatethousands(nTotal,',',0)], 'FontSize', 7, 'VerticalAlignment' , 'middle', 'HorizontalAlignment' , 'right', 'BackgroundColor', 'White');

        set(gca, ...
          'TickDir'     , 'in'         , ...
          'TickLength'  , [0.01 0]         , ...
          'XMinorTick'  , 'off'         , ...
          'YMinorTick'  , 'off'         , ...
          'YGrid'       , 'on'          , ...
          'GridColor'   , [0.5 0.5 0.5] , ...
          'GridAlpha'   , 1             , ...
          'YTick'       , -100:20:100);    
    end

    export_fig('filename', fileSaveName1, '-nocrop');
    close(gcf);
end

