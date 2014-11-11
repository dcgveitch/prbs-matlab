for perm=4
    clear d_*;
    d_nZones=2;
    nX=1;
    nY=4.5;

    sizeXMOff=0.07;
    sizeYMOff=0.07;
    sizeXTotal=0.9;
    sizeYTotal=0.78;
    sizeXPad=0;
    sizeYPad=0;

%     for i=1:length(nX)
%         sizeXM(i)=(sizeXTotal-(nX(i)-1)*sizeXPad)/nX(i); 
%     end
    sizeXM=(sizeXTotal-(nX-1)*sizeXPad)/nX; 
    sizeYM(1:4)=(sizeYTotal-(nY-1)*sizeYPad)/nY;
    sizeYM(1)=sizeYM(1)*1.5;

    colours=pmkmp(6,'CubicL');  

    figOut=figure;
    set(figOut,'Units','centimeters');
    set(figOut,'Position', [5, 5, 17, 7]);
    set(figOut,'Units','pixels');

    d_flowCount=1;

    for d_flowType=1:3
        d_flow=out_varFlow;
        if d_flowType==1
            d_flowRef=1;
            handAxesM(d_flowType,d_flowRef)=axes;
            d_plot=d_flow{1,1}{1,1}(:,2);
            plot(d_plot);
            set(handAxesM(d_flowType,d_flowRef),'color','none', 'XTickLabel', [], 'YTickLabel', [],'Position', [sizeXMOff 1-(sizeYMOff+sum(sizeYM(1:d_flowCount))) sizeXM sizeYM(d_flowCount)]);
            ylim([0 180*1.5]);
            xlim([0 length(d_plot)]);
            set(gca,'xTick', 0:2880:40320);
            set(gca,'yTick', [0 180*1.5]);
            grid on;
            d_flowCount=d_flowCount+1;
        end
        if d_flowType==2
            for d_flowRef=1:size(d_flow{2,1},2)
                handAxesM(d_flowType,d_flowRef)=axes;
                d_plot=d_flow{2,1}{1,d_flowRef}(:,2);
                plot(d_plot,'r');
                hold on;
                zoneFlow=plot(d_flow{2,1}{2,d_flowRef}(:,2));
                set(handAxesM(d_flowType,d_flowRef),'color','none', 'XTickLabel', [], 'YTickLabel', [], 'Position', [sizeXMOff 1-(sizeYMOff+sum(sizeYM(1:d_flowCount))) sizeXM sizeYM(d_flowCount)]);
                ylim([0 180]);
                xlim([0 length(d_plot)]);
                set(gca,'xTick', 0:2880:40320);
                set(gca,'yTick', [0 180]);
                grid on;
                d_flowCount=d_flowCount+1;
            end
        end    
        if d_flowType==3
            d_indices=reshape(1:d_nZones^2,d_nZones,d_nZones)';
            d_indicePairs(:,1)=nonzeros(tril(d_indices,-1));
            d_indicePairs(:,2)=nonzeros(triu(d_indices,1)');
            for d_flowRef=1:size(d_indicePairs,1)
                handAxesM(d_flowType,d_flowRef)=axes;
                d_plot=d_flow{3,1}{floor((d_indicePairs(d_flowRef,1)-1)/d_nZones)+1,rem(d_indicePairs(d_flowRef,1)-1,d_nZones)+1}(:,2);
                plot(d_plot,'r');
                hold on;
                plot(d_flow{3,1}{floor((d_indicePairs(d_flowRef,2)-1)/d_nZones)+1,rem(d_indicePairs(d_flowRef,2)-1,d_nZones)+1}(:,2));
                if size(d_indicePairs,1)==1, d_flowPos=2; else d_flowPos=d_flowRef; end
                set(handAxesM(d_flowType,d_flowRef),'color','none', 'XTickLabel', [], 'YTickLabel', [], 'Position', [sizeXMOff 1-(sizeYMOff+sum(sizeYM(1:d_flowCount))) sizeXM sizeYM(d_flowCount)]);
                ylim([0 180]);
                xlim([0 length(d_plot)]);
                set(gca,'xTick', 0:2880:40320);
                set(gca,'yTick', [0 180]);
                grid on;
                d_flowCount=d_flowCount+1;
            end
            set(gca,'xTick', 0:2880:40320, 'XTickLabel',{'0' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' })
            hXLabel = xlabel('Time (days)');
        end
    end
    legendflex(handAxesM(2,1),{'Exfiltration' 'Infiltration'}, 'ref', handAxesM(1,1));

    fileSaveName=['Perm_' num2str(perm) '.pdf'];
    export_fig('filename', fileSaveName, '-nocrop');
%     close(gcf);
end
        
% 
%             
%     
% 
% for d_i=1:3
%     P7_Grouping;    
%     d_summary=out_summary(find(out_summary(:,1)==d_i),:);
%     flowStats=d_summary(:,11:15)'*100;
%     nTotal=sum(d_summary(:,6));
%     grouping1=d_summary(:,3)';
%     grouping2=d_summary(:,2)';
%     
%     handAxesM(d_i)=axes;
%     h=boxplotDV(data,{grouping2,grouping1},'factorgap',[3 0.1],'colors',repmat(flipud(colours(2:6,:)),6,1),'factorseparator',[1], 'widths', 1, 'boxstyle','outline');
%     set(h(1,:),{'Ydata'},num2cell(flowStats(end-1:end,:),1)')
%     set(h(2,:),{'Ydata'},num2cell(flowStats(2:-1:1,:),1)')
%     set(h(3,:),{'Ydata'},num2cell(flowStats([end end],:),1)')
%     set(h(4,:),{'Ydata'},num2cell(flowStats([1 1],:),1)')
%     set(h(5,:),{'Ydata'},num2cell(flowStats([2 4 4 2 2],:),1)')
%     set(h(6,:),{'Ydata'},num2cell(flowStats([3 3],:),1)')
%     set(h(7,:),{'Visible'},{'off'})
%     ylim([-100 100]);
%     set(handAxesM(d_i),'color','none','Position', [sizeXMOff+rem(d_i-1,3)*(sizeXM+sizeXPad) 1-(sizeYMOff+(ceil(d_i/3))*(sizeYM+sizeYPad)) sizeXM sizeYM]);
%     % set(handAxesM(j),'XTickLabel',[],'YTickLabel',[], 'Position', [sizeXMOff+rem(j-1,3)*sizeXM 1-(sizeYMOff+(ceil(j/3))*sizeYM) sizeXM sizeYM], 'color', 'none');
% 
%     h = findobj(gca,'Tag','Box');
%     for x=1:length(h)
%         boxShade=patch(get(h(x),'XData'),get(h(x),'YData'), colours(rem(x-1,5)+2,:), 'linestyle', 'none');
% %                     boxShade=patch(get(h(x),'XData'),get(h(x),'YData'), colours(rem(x-1,2)+2,:), 'linestyle', 'none');
%         uistack(boxShade,'bottom');
%     end
%     delete(h);
%     set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-','LineWidth',1)
%     set(findobj(gcf,'Tag','Upper Adjacent Value'),'LineStyle','-','LineWidth',1)
%     set(findobj(gcf,'Tag','Lower Adjacent Value'),'LineStyle','-','LineWidth',1)
%     set(findobj(gcf,'Tag','Median'),'Color','Red','LineWidth',0.5)
% 
%     graphX=max(xlim);
%     graphXint=graphX/6;
% 
%     set(gca,'xTick', graphXint/2:graphXint:graphX, 'XTickLabel',{'2' '3' '4' '6' '12' '24'})
% %             set(gca,'FontName','Helvetica','FontSize',16);
% 
%     switch d_i
%         case 1
%             title('Total Dwelling Flow');
%             hYLabel = ylabel('Flow Weighted Error (%)');
%         case 2
%             title('Total Zonal Flows');
%             set(gca,'YTickLabel',[]);
%             hXLabel = xlabel('PRBS Sequence Period (hours)');
%             xlabh = get(gca,'XLabel');
%             get(xlabh,'Position');
%         case 3
%             title('Individual Flows');
%             set(gca,'YTickLabel',[]);
%             hleg = legend('15 bit','31 bit','63 bit','127 bit','255 bit');
%             set(hleg,'FontSize',7);
%     end
%     text((max(xlim)-min(xlim))*0.8+min(xlim), (max(ylim)-min(ylim))*0.05+min(ylim), ['n = ' separatethousands(nTotal,',',0)], 'FontSize', 7, 'VerticalAlignment' , 'middle', 'HorizontalAlignment' , 'center', 'BackgroundColor', 'White');
% 
%     set(gca, ...
%       'TickDir'     , 'out'     , ...
%       'TickLength'  , [0 0] , ...
%       'XMinorTick'  , 'off'      , ...
%       'YMinorTick'  , 'off'      , ...
%       'YGrid'       , 'on'      , ...
%       'YTick'       , -100:20:100, ...
%       'LineWidth'   , 0.5         );
%   
%   fileSaveName=['Data_' num2str(d_i) '.xlsx'];
%   xlswrite(fileSaveName,d_summary);
% end
% fileSaveName=['Test.pdf'];
% export_fig('filename', fileSaveName, '-nocrop');
% close(gcf);  

