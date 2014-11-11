%% Graphing Code
data=rand(20,24);

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

groupDims=[7 9 1 2 6];

d_fig=4;
d_reqSeqLength=1:length(unique(r_seqLength));
d_reqSeqPeriod=1:length(unique(r_seqPeriod));
d_reqNZones=1:length(unique(r_nZones));
d_reqSolve=[1];
d_reqImp=[1];
d_reqNSeqA=[1 2 3 4];
d_reqTSeqA=[1 2];
d_reqConc=[3 4];
d_reqFlowType=[1 2 3];

for d_zone=1
    if d_zone~=1, d_reqNZones=d_zone-1; end
    for d_tSeqA=2
        d_reqTSeqA=d_tSeqA;
        for d_conc=4
            d_reqConc=d_conc;
            P7_GroupingN;
            
            d_summary=sortrows(out_summary,[17 18]);

            %% Add in stats for empty boxplots

            figOut=figure;
            set(figOut,'Units','centimeters');
            set(figOut,'Position', [5, 5, 17, 7]);
            set(figOut,'Units','pixels');

            set(gcf,'Renderer','Painters');    

            for d_i=1:3
                d_summaryF=d_summary(find(d_summary(:,2)==d_i),:);
                d_summaryF=[d_summaryF(1,:); d_summaryF(1,:); d_summaryF(1,:);d_summaryF(1,:); d_summaryF(1,:); d_summaryF(1,:); d_summaryF(2:end,:)];
                d_summaryF=[d_summaryF(1:8,:); d_summaryF(8,:); d_summaryF(9:end,:)];
                d_summaryF=[d_summaryF(1:10,:); d_summaryF(10,:); d_summaryF(10,:); d_summaryF(11:end,:)];
                d_summaryF=[d_summaryF(1:17,:); d_summaryF(17,:); d_summaryF(18:end,:)];
                d_summaryF([2 3 4 5 6 9 11 12 18],:)=0;
                d_summaryF=d_summaryF(1:24,:);
                d_summaryF(:,3)=repmat([1; 2; 3],8 ,1);
                d_summaryF(:,4)=repmat([1; 1; 1; 2; 2; 2],4 ,1);
                d_summaryF(1:6,17)=2;
                d_summaryF(7:12,17)=4;
                d_summaryF(13:18,17)=8;

                flowStats=d_summaryF(:,11:15)'*100;
                nTotal=sum(d_summaryF(:,6));
                grouping1=d_summaryF(:,3)';
                grouping2=d_summaryF(:,17)';
                grouping3=d_summaryF(:,4)';
                
                handAxesM(d_i)=axes;
                h=boxplotDV(data,{grouping2, grouping3, grouping1},'factorgap',[5.74 5.74 0.1],'colors',repmat(flipud(colours(2:4,:)),8,1),'factorseparator',[1], 'widths', 1, 'boxstyle','outline');
                set(h(1,:),{'Ydata'},num2cell(flowStats(end-1:end,:),1)')
                set(h(2,:),{'Ydata'},num2cell(flowStats(2:-1:1,:),1)')
                set(h(3,:),{'Ydata'},num2cell(flowStats([end end],:),1)')
                set(h(4,:),{'Ydata'},num2cell(flowStats([1 1],:),1)')
                set(h(5,:),{'Ydata'},num2cell(flowStats([2 4 4 2 2],:),1)')
                set(h(6,:),{'Ydata'},num2cell(flowStats([3 3],:),1)')
                set(h(7,:),{'Visible'},{'off'})
                ylim([-100 100]);
                set(handAxesM(d_i),'color','none','Position', [sizeXMOff+rem(d_i-1,3)*(sizeXM+sizeXPad) 1-(sizeYMOff+(ceil(d_i/3))*(sizeYM+sizeYPad)) sizeXM sizeYM]);
                set(gca, 'Layer','top')

                h = findobj(gca,'Tag','Box');
                for x=1:length(h)
                    boxShade(rem(x-1,3)+1)=patch(get(h(x),'XData'),get(h(x),'YData'), colours(rem(x-1,3)+2,:), 'linestyle', 'none');
                    uistack(boxShade(rem(x-1,3)+1),'bottom');
                end

                delete(h);
                set(findobj(gca,'-regexp','Tag','\w*Whisker'),{'LineStyle'},{'-','-','-','-','-','-','none','-','-','-','-','-','none','none','-','none','-','-','none','none','none','none','none','-','-','-','-','-','-','-','none','-','-','-','-','-','none','none','-','none','-','-','none','none','none','none','none','-'}','LineWidth',1);
                set(findobj(gca,'Tag','Upper Adjacent Value'),{'LineStyle'},{'-','-','-','-','-','-','none','-','-','-','-','-','none','none','-','none','-','-','none','none','none','none','none','-'}','LineWidth',1)
                set(findobj(gca,'Tag','Lower Adjacent Value'),{'LineStyle'},{'-','-','-','-','-','-','none','-','-','-','-','-','none','none','-','none','-','-','none','none','none','none','none','-'}','LineWidth',1)
                set(findobj(gca,'Tag','Median'),'LineStyle','-',{'Color'},{'Red','Red','Red','Red','Red','Red','Black','Red','Red','Red','Red','Red','Black','Black','Red','Black','Red','Red','Black','Black','Black','Black','Black','Red'}','LineWidth',0.5)

                graphX=max(xlim);
                graphXint=graphX/4;
                
                for x=1:4
                    dtShade=patch([graphXint*(x-0.5) graphXint*(x-0.5) graphXint*x graphXint*x], [-100 100 100 -100], [0.9 0.9 0.9], 'linestyle', 'none');
                    uistack(dtShade,'bottom');
                end

                set(gca,'xTick', graphXint/2:graphXint:graphX, 'XTickLabel',{'2' '4' '8' '16'})
                
                text(graphXint*1/4,90,{'dt=','4m'},'HorizontalAlignment', 'Center', 'FontSize', 7, 'FontAngle', 'italic');
                text(graphXint*3/4,90,{'dt=','8m'},'HorizontalAlignment', 'Center', 'FontSize', 7, 'FontAngle', 'italic');

                switch d_i
                    case 1
                        title('Total Dwelling Flow');
                        hYLabel = ylabel('Flow Weighted Error (%)');
                    case 2
                        title('Total Zonal Flows');
                        set(gca,'YTickLabel',[]);
                        hXLabel = xlabel('Data Analysis Period (hours)');
                        xlabh = get(gca,'XLabel');
                        get(xlabh,'Position');
                    case 3
                        title('Individual Flows');
                        set(gca,'YTickLabel',[]);
                        hleg = legend(boxShade, '15 bit','31 bit','63 bit');
                        set(hleg,'FontSize',7);
                end
                text((max(xlim)-min(xlim))*0.68+min(xlim), (max(ylim)-min(ylim))*0.07+min(ylim), ['n = ' separatethousands(nTotal,',',0)], 'FontSize', 7,'VerticalAlignment' , 'top');

                set(gca, ...
                  'TickDir'     , 'out'     , ...
                  'TickLength'  , [0 0] , ...
                  'XMinorTick'  , 'off'      , ...
                  'YMinorTick'  , 'off'      , ...
                  'YGrid'       , 'on'      , ...
                  'YTick'       , -100:20:100, ...
                  'LineWidth'   , 0.5         );
              
            fileSaveName=['Data_' num2str(d_i) '.xlsx'];
            xlswrite(fileSaveName,d_summary);
            
            end
            fileSaveName=[num2str(d_zone) '_' num2str(d_tSeqA) '_' num2str(d_conc) '_Test.pdf'];
            export_fig('filename', fileSaveName, '-nocrop');
            close(gcf);
        end
    end
end


