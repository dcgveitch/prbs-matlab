%% Graphing Code
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

groupDims=[9 2 1];

d_fig=7;
d_reqSeqLength=1:length(unique(r_seqLength(1:38)));
d_reqSeqPeriod=1:length(unique(r_seqPeriod(1:38)));
d_reqNZones=[1 2 3];
d_reqSolve=[];
d_reqImp=[];
d_reqNSeqA=[1];
d_reqTSeqA=[1];
d_reqConc=[];
d_reqFlowType=[1 2 3];

for d_zone=1
    if d_zone~=1, d_reqNZones=d_zone-1; end
    for d_solve=[1 4 5]
        d_reqSolve=d_solve;
        for d_conc=[3 4]
            d_reqConc=d_conc;
            for d_imp=2
                d_reqImp=d_imp;

                P7_Grouping;
                d_summary=out_summary;

                %% Add in stats for empty boxplots

                figOut=figure;
                set(figOut,'Units','centimeters');
                set(figOut,'Position', [5, 5, 17, 7]);
                set(figOut,'Units','pixels');

                set(gcf,'Renderer','Painters');    

                for d_i=1:3
                    d_summaryF=d_summary(find(d_summary(:,1)==d_i),:);

                    flowStats=d_summaryF(:,11:15)'*100;
                    nTotal=sum(d_summaryF(:,6));
                    grouping1=d_summaryF(:,2)';
                    grouping2=d_summaryF(:,3)';
                    
                    data=rand(20,size(flowStats,2));

                    handAxesM(d_i)=axes;
                    h=boxplotDV(data,{grouping1, grouping2},'factorgap',[5.74 0.01],'colors',repmat(flipud(colours(5:6,:)),2,1),'factorseparator',[1], 'widths', 1, 'boxstyle','outline');
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
                        boxShade(rem(x-1,3)+1)=patch(get(h(x),'XData'),get(h(x),'YData'), colours(rem(x-1,2)+5,:), 'linestyle', 'none');
                        uistack(boxShade(rem(x-1,3)+1),'bottom');
                    end

                    delete(h);
                    set(findobj(gca,'-regexp','Tag','\w*Whisker'),'LineStyle','-','LineWidth',1);
                    set(findobj(gca,'Tag','Upper Adjacent Value'),'LineStyle','-','LineWidth',1)
                    set(findobj(gca,'Tag','Lower Adjacent Value'),'LineStyle','-','LineWidth',1)
                    set(findobj(gca,'Tag','Median'),'LineStyle','-','Color','Red','LineWidth',0.5)

                    graphX=max(xlim);
                    graphXint=graphX/6;

                    set(gca,'xTick', graphXint/2:graphXint:graphX, 'XTickLabel',{'2' '3' '4' '6' '12' '24'})


                    switch d_i
                        case 1
                            title('Total Dwelling Flow');
                            hYLabel = ylabel('Flow Weighted Error (%)');
                        case 2
                            title('Total Zonal Flows');
                            set(gca,'YTickLabel',[]);
                            hXLabel = xlabel('PRBS Period (hours)');
                            xlabh = get(gca,'XLabel');
                            get(xlabh,'Position');
                        case 3
                            title('Individual Flows');
                            set(gca,'YTickLabel',[]);
                            hleg = legend(boxShade, '15 bit','31 bit');
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
                fileSaveName=[num2str(d_zone) '_' num2str(d_solve) '_' num2str(d_conc) '_' num2str(d_imp)' '_Test.pdf'];
                export_fig('filename', fileSaveName, '-nocrop');
                close(gcf);
            end
        end
    end
end


