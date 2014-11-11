for d_zone=[1 2 3 4 5]
    for d_imp=[2]
        for d_conc=[3 4]
            g_data=out_groupResults{1,d_zone}{d_imp,d_conc};
            periods=unique(g_data(:,2));
            
            if d_zone==5
                g_data=[15 2 zeros(1,33); g_data];
                for i=5:-1:1
                    g_data=[g_data(1:i*4+1,:); 15 periods(i+1) zeros(1,33); g_data(i*4+2:end,:)];
                end
            end
            
            grouping1=g_data(:,1)';
            grouping2=g_data(:,2)';
            
            data=rand(20,30);

            flowStats=cell(1,3);
            nTotal(1)=sum(g_data(:,3));
            flowStats{1}=g_data(:,8:12)'*100;
            nTotal(2)=sum(g_data(:,14));
            flowStats{2}=g_data(:,19:23)'*100;
            nTotal(3)=sum(g_data(:,25));
            flowStats{3}=g_data(:,30:34)'*100;

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

            figOut=figure;
            set(figOut,'Units','centimeters');
            set(figOut,'Position', [5, 5, 17, 7]);
            set(figOut,'Units','pixels');

            set(gcf,'Renderer','Painters');    

            for d_i=1:3
                handAxesM(d_i)=axes;
    %             h=boxplot(data,{grouping2,grouping1},'factorgap',[3 0.1],'colors',repmat(flipud(colours(2:6,:)),6,1),'factorseparator',[1], 'widths', 1, 'boxstyle','outline');
                h=boxplot(data,{grouping2,grouping1},'factorgap',[3 0.1],'colors',repmat(flipud(colours(2:6,:)),6,1),'factorseparator',[1], 'widths', 1, 'boxstyle','outline');
                set(h(1,:),{'Ydata'},num2cell(flowStats{d_i}(end-1:end,:),1)')
                set(h(2,:),{'Ydata'},num2cell(flowStats{d_i}(2:-1:1,:),1)')
                set(h(3,:),{'Ydata'},num2cell(flowStats{d_i}([end end],:),1)')
                set(h(4,:),{'Ydata'},num2cell(flowStats{d_i}([1 1],:),1)')
                set(h(5,:),{'Ydata'},num2cell(flowStats{d_i}([2 4 4 2 2],:),1)')
                set(h(6,:),{'Ydata'},num2cell(flowStats{d_i}([3 3],:),1)')
                set(h(7,:),{'Visible'},{'off'})
                ylim([-100 100]);
                set(handAxesM(d_i),'color','none','Position', [sizeXMOff+rem(d_i-1,3)*(sizeXM+sizeXPad) 1-(sizeYMOff+(ceil(d_i/3))*(sizeYM+sizeYPad)) sizeXM sizeYM]);
                % set(handAxesM(j),'XTickLabel',[],'YTickLabel',[], 'Position', [sizeXMOff+rem(j-1,3)*sizeXM 1-(sizeYMOff+(ceil(j/3))*sizeYM) sizeXM sizeYM], 'color', 'none');

                h = findobj(gca,'Tag','Box');
                for x=1:length(h)
                    boxShade=patch(get(h(x),'XData'),get(h(x),'YData'), colours(rem(x-1,5)+2,:), 'linestyle', 'none');
%                     boxShade=patch(get(h(x),'XData'),get(h(x),'YData'), colours(rem(x-1,2)+2,:), 'linestyle', 'none');
                    uistack(boxShade,'bottom');
                end
                delete(h);
                set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-','LineWidth',1)
                set(findobj(gcf,'Tag','Upper Adjacent Value'),'LineStyle','-','LineWidth',1)
                set(findobj(gcf,'Tag','Lower Adjacent Value'),'LineStyle','-','LineWidth',1)
                set(findobj(gcf,'Tag','Median'),'Color','Red','LineWidth',0.5)

                graphX=max(xlim);
                graphXint=graphX/6;

                set(gca,'xTick', graphXint/2:graphXint:graphX, 'XTickLabel',{'2' '3' '4' '6' '12' '24'})
    %             set(gca,'FontName','Helvetica','FontSize',16);

                switch d_i
                    case 1
                        title('Total Dwelling Flow');
                        hYLabel = ylabel('Flow Weighted Error (%)');
                        hXLabel = xlabel('PRBS Sequence Period (hours)');
                    case 2
                        title('Total Zonal Flows');
                        set(gca,'YTickLabel',[]);
                        hXLabel = xlabel('PRBS Sequence Period (hours)');
                    case 3
                        title('Individual Flows');
                        set(gca,'YTickLabel',[]);
                        hXLabel = xlabel('PRBS Sequence Period (hours)');
                        hleg = legend('15 bit','31 bit','63 bit','127 bit','255 bit');
                        set(hleg,'FontSize',7);
                end
                text((max(xlim)-min(xlim))*0.68+min(xlim), (max(ylim)-min(ylim))*0.07+min(ylim), ['n = ' separatethousands(nTotal(d_i),',',0)], 'FontSize', 7,'VerticalAlignment' , 'top');

                set(gca, ...
                  'TickDir'     , 'out'     , ...
                  'TickLength'  , [0 0] , ...
                  'XMinorTick'  , 'off'      , ...
                  'YMinorTick'  , 'off'      , ...
                  'YGrid'       , 'on'      , ...
                  'YTick'       , -100:20:100, ...
                  'LineWidth'   , 0.5         );
            end
            fileSaveName=[num2str(d_conc) '_' num2str(d_imp) '_' num2str(d_zone) '.pdf'];
            export_fig('filename', fileSaveName, '-nocrop');
            close(gcf);  
        end
    end
end


