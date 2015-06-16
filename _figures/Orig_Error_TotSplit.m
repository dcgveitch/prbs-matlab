% Group 1 (Sub plot) = Number of Zones
% Group 2 (Colour) = Sequence Length

id_perm=1; % Run identifier for drawing connecting lines between markers
id_1=2; % Page Break
id_2=3; % Sub Plot
id_series=4; % Series identifier for different colors/lines

% Percentile based statistics
id_x1=7; % Sub Plot X Axis (Processed)
id_y1=24; % Sub Plot Y Axis (Processed)
id_eL1=22; % Bounded area 1 Lower
id_eL2=23; % Bounded area 2 Lower
id_eU2=25; % Bounded area 2 Upper
id_eU1=26; % Bounded area 1 Upper

id_x2=8; % Sub Plot X Axis (Raw)
id_y2=9; % Sub Plot Y Axis (Raw)
id_w2=11; % Sub Plot Weighting (Raw)

d_reqSolve=[1];
d_reqImp=[2];
d_reqConc=[4];
d_reqFlowType=[2];

d_seqA=1;

nX=3;
nY=2;

sizeXMOff=0.05;
sizeYMOff=0.1;
sizeXM=0.9/nX;
sizeYM=0.85/nY;
sizeXm=sizeXM*0.4;
sizeYm=sizeYM/10;
sizeXmOff=sizeXM*0.5;
sizeYmOff=sizeYM*0.08;
xOffset=0.003;

cmap=[0.5 0.5 0.5; 0 0 0];

for d_solve=d_reqSolve
    for d_impulse=d_reqImp
        for d_conc=d_reqConc
            for d_flowType=d_reqFlowType
                id_group1=unique(out_indResultsProc{d_impulse,d_conc}{d_solve,d_flowType}(:,id_1));
                id_group2=unique(out_indResultsProc{d_impulse,d_conc}{d_solve,d_flowType}(:,id_2));
                id_groupPermM=unique(out_indResultsProc{d_impulse,d_conc}{d_solve,d_flowType}(:,id_series));
                seriesColour=varycolor(length(id_groupPermM));
                d_pos=1;

                for zone=1:length(id_groupPermM)
                    for i=1:length(id_group1)
                        handaxesM=[];
                        handaxesB=[];
                        count=0;
                        figure('Position', [100, 100, 1600, 800]);
                        set(gcf,'Renderer','Painters');

                        switch d_flowType
                            case 1, desFlow='Total';
                            case 2, desFlow='ZoInf';
                            case 3, desFlow='Individual';
                        end

                        switch d_conc
                            case 1, desConc='Theory';
                            case 2, desConc='1st Order';
                            case 3, desConc='Noise';
                            case 4, desConc='Sensor';
                        end

                        switch d_impulse
                            case 1, desImp='Standard';
                            case 2, desImp='Combined';
                            case 3, desImp='Direct';
                        end

                        for j=1:length(id_group2)
                            d_graphSelect=out_indResultsProc{d_impulse,d_conc}{d_solve,d_flowType}(:,id_1)==id_group1(i) & out_indResultsProc{d_impulse,d_conc}{d_solve,d_flowType}(:,id_2)==id_group2(j) & out_indResultsProc{d_impulse,d_conc}{d_solve,d_flowType}(:,id_series)==id_groupPermM(zone);
                            d_boxSelect=out_indResultsRaw{d_impulse,d_conc}{d_solve,d_flowType}(:,id_1)==id_group1(i) & out_indResultsRaw{d_impulse,d_conc}{d_solve,d_flowType}(:,id_2)==id_group2(j) & out_indResultsRaw{d_impulse,d_conc}{d_solve,d_flowType}(:,id_series)==id_groupPermM(zone);

                            if(~any(d_graphSelect))
                                continue;
                            end
                            
                            d_graph=out_indResultsProc{d_impulse,d_conc}{d_solve,d_flowType}(d_graphSelect,:);
                            d_box=out_indResultsRaw{d_impulse,d_conc}{d_solve,d_flowType}(d_boxSelect,:);
                            
                            

                            clear data*;
                            handAxesM(j) = axes('color','none');
                            dataMean(:,1)=d_graph(:,id_x1);
                            dataMean(:,2)=d_graph(:,id_y1);
                            dataMean(:,3)=d_graph(:,id_y1)-d_graph(:,id_eL1);
                            dataMean(:,4)=d_graph(:,id_eU1)-d_graph(:,id_y1);
                            dataMean(:,5)=d_graph(:,id_y1)-d_graph(:,id_eL2);
                            dataMean(:,6)=d_graph(:,id_eU2)-d_graph(:,id_y1);
                            dataMean(:,7)=d_graph(:,id_series);
                            dataMean(:,8)=d_graph(:,id_perm);
                            dataMeanSort=sortrows(dataMean,1);        
                            dataAll(:,1)=d_box(:,id_x2);
                            dataAll(:,2)=d_box(:,id_y2);
                            dataAll(:,3)=(dataAll(:,2)-dataAll(:,1))./dataAll(:,1);
                            dataAll(:,4)=d_box(:,id_w2);
                            boundedline(dataMeanSort(:,1),dataMeanSort(:,2),dataMeanSort(:,3:4),'b.',dataMeanSort(:,1),dataMeanSort(:,2),dataMeanSort(:,5:6),'k.', 'cmap', cmap);
                            hold on;
                            id_groupPerm=unique(dataMeanSort(:,8));
                            for k=length(id_groupPerm):-1:1
                                d_range=find(dataMeanSort(:,8)==id_groupPerm(k));
                                d_color=find(id_groupPermM==dataMeanSort(d_range(1),7));
                                plot(dataMeanSort(d_range,1), dataMeanSort(d_range,2), '-o', 'MarkerSize', 2, 'MarkerFaceColor', seriesColour(d_color,:), 'Color', seriesColour(d_color,:));
                            end
                            colormap(jet);
                            axis([0 400 0 400]);
                            set(handAxesM(j),'XTickLabel',[],'YTickLabel',[], 'Position', [sizeXMOff+rem(j-1,nX)*sizeXM 1-(sizeYMOff+(ceil(j/nX))*sizeYM) sizeXM sizeYM], 'color', 'none');
                            line(xlim, ylim,'Color', [0 0 0], 'LineWidth', 0.5);
                            dashline(xlim, ylim.*1.1, 1,1,1,1,'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 0.5);
                            dashline(xlim, ylim.*0.9, 1,1,1,1,'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 0.5);
                            count=size(dataMeanSort,1);
                            text((max(xlim)-min(xlim))*0.04+min(xlim), (max(ylim)-min(ylim))*0.90+min(ylim), [num2str(id_group2(j)) ' hours'], 'FontSize', 10, 'EdgeColor', 'black', 'VerticalAlignment' , 'bottom');
                            text((max(xlim)-min(xlim))*0.04+min(xlim), (max(ylim)-min(ylim))*0.89+min(ylim), ['n = ' num2str(count)], 'FontSize', 7, 'VerticalAlignment' , 'top');
                            text((max(xlim)-min(xlim))*0.04+min(xlim), (max(ylim)-min(ylim))*0.84+min(ylim), ['wmean = ' num2str(wmean(dataAll(:,3),dataAll(:,4))*100,3) ' %'], 'FontSize', 7, 'VerticalAlignment' , 'top');
                            text((max(xlim)-min(xlim))*0.04+min(xlim), (max(ylim)-min(ylim))*0.79+min(ylim), ['wsd = ' num2str(wstd(dataAll(:,3),dataAll(:,4))*100,3) ' %'], 'FontSize', 7, 'VerticalAlignment' , 'top');
                            text((max(xlim)-min(xlim))*0.04+min(xlim), (max(ylim)-min(ylim))*0.74+min(ylim), ['wmed = ' num2str(wprctile(dataAll(:,3),50,dataAll(:,4))*100,3) ' %'], 'FontSize', 7, 'VerticalAlignment' , 'top');
                            text((max(xlim)-min(xlim))*0.04+min(xlim), (max(ylim)-min(ylim))*0.69+min(ylim), ['wpersd = ' num2str(wprctile(dataAll(:,3),15.87,dataAll(:,4))*100,3) ' - ' num2str(wprctile(dataAll(:,3),84.13,dataAll(:,4))*100,3) ' %'], 'FontSize', 7, 'VerticalAlignment' , 'top');
                            text((max(xlim)-min(xlim))*0.04+min(xlim), (max(ylim)-min(ylim))*0.64+min(ylim), ['wIQR = ' num2str(wprctile(dataAll(:,3),25,dataAll(:,4))*100,3) ' - ' num2str(wprctile(dataAll(:,3),75,dataAll(:,4))*100,3) ' %'], 'FontSize', 7, 'VerticalAlignment' , 'top');
                            if (j==1)
                                text((max(xlim)-min(xlim))*0.04+min(xlim), (max(ylim)-min(ylim))*0.125+max(ylim), [num2str(id_groupPermM(zone)) ' Zone ' desFlow ' - ' desImp ' - ' desConc], 'FontSize', 13, 'VerticalAlignment' , 'bottom');
                                text((max(xlim)-min(xlim))*0.04+min(xlim), (max(ylim)-min(ylim))*0.125+max(ylim), [num2str(id_group1(i)) ' bit sequence'], 'FontSize', 11, 'VerticalAlignment' , 'top');
                            end

                            d_stats(1,1)=id_group1(i);
                            d_stats(1,2)=id_group2(j);
                            d_stats(1,3)=id_groupPermM(zone);
                            d_stats(1,4)=length(dataAll);
                            d_stats(1,5)=mean(dataAll(:,3));
                            d_stats(1,6)=std(dataAll(:,3));
                            d_stats(1,7)=prctile(dataAll(:,3),5);
                            d_stats(1,8)=prctile(dataAll(:,3),15.87);
                            d_stats(1,9)=prctile(dataAll(:,3),25);
                            d_stats(1,10)=prctile(dataAll(:,3),50);
                            d_stats(1,11)=prctile(dataAll(:,3),75);
                            d_stats(1,12)=prctile(dataAll(:,3),84.13);
                            d_stats(1,13)=prctile(dataAll(:,3),95);
                            d_stats(1,14)=length(dataAll);
                            d_stats(1,15)=wmean(dataAll(:,3),dataAll(:,4));
                            d_stats(1,16)=wstd(dataAll(:,3),dataAll(:,4));
                            d_stats(1,17)=wprctile(dataAll(:,3),5,dataAll(:,4));
                            d_stats(1,18)=wprctile(dataAll(:,3),15.87,dataAll(:,4));
                            d_stats(1,19)=wprctile(dataAll(:,3),25,dataAll(:,4));
                            d_stats(1,20)=wprctile(dataAll(:,3),50,dataAll(:,4));
                            d_stats(1,21)=wprctile(dataAll(:,3),75,dataAll(:,4));
                            d_stats(1,22)=wprctile(dataAll(:,3),84.13,dataAll(:,4));  
                            d_stats(1,23)=wprctile(dataAll(:,3),95,dataAll(:,4));   

                            out_stats{d_impulse,d_conc}{d_solve,d_flowType}(d_pos,:)=d_stats;


                            % Boxplot
                            handaxesB(i) = axes;
                            h=boxplot(dataAll(:,3), 'orientation', 'horizontal', 'width', 0.5, 'outliersize', 3.5, 'jitter', 0.5, 'symbol', 'k.');
                            % Replace standard values with weighted statistcs
                            set(h(1,1),'Xdata',d_stats([21 23]));
                            set(h(2,1),'Xdata',d_stats([17 19]));
                            set(h(3,1),'Xdata',d_stats([23 23]));
                            set(h(4,1),'Xdata',d_stats([17 17]));
                            set(h(5,1),'Xdata',d_stats([19 21 21 19 19]));
                            set(h(6,1),'Xdata',d_stats([20 20]));
                            set(h(7,1),{'Visible'},{'off'})
                            xlim([-0.5 0.5]);
                            h = findobj(handaxesB(i),'Tag','Box');
                            for x=1:length(h)
                                boxShade=patch(get(h(x),'XData'),get(h(x),'YData'), [0.7,0.7,0.7], 'linestyle', 'none');
                                uistack(boxShade,'bottom'); 
                            end
                            delete(h);
                            set(findobj(handaxesB(i),'Tag','Median'),'Color',[0 0 0]);
                            set(findobj(handaxesB(i),'Tag','Outliers'),'MarkerEdgeColor',[0 0 0]);
                            trueLine=line([0 0], ylim, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5);
                            uistack(trueLine,'bottom');   
                            set(handaxesB(i),'YTickLabel',[], 'XTick', -0.5:0.1:0.5, 'XTickLabel', {'-50%' '' '' '' '' '0' '' '' '' '' '+50%'}, 'FontSize', 7, 'Color', [0.95 0.95 0.95],'Position', [sizeXMOff+rem(j-1,nX)*sizeXM+sizeXmOff 1-(sizeYMOff+(ceil(j/nX))*sizeYM)+sizeYmOff sizeXm sizeYm]);

                            % Histogram
                            handaxesH(i) = axes;
                            % Weighted histgram
                            [N,binCenters] = histwc(dataAll(:,3),dataAll(:,4),-0.5:0.05:0.5);         
                            hBar = bar(binCenters,N,'hist');
                            index = ones(1,length(N));
                            index(1) = 2;
                            index(end) = 2;
                            set(hBar,'FaceVertexCData',index');
                    %         colormap([0 0 0.5; 1 0 0]);
                            axis tight;
                            xlim([-0.5 0.5]);
                            set(handaxesH(i),'Position', [sizeXMOff+rem(j-1,nX)*sizeXM+sizeXmOff 1-(sizeYMOff+(ceil(j/nX))*sizeYM)+sizeYmOff+sizeYm sizeXm sizeYm*1.5]);
                            axis off

                            d_pos=d_pos+1;
                        end
                        
                        if(isempty(get(gca, 'children')))
                            continue;
                        end                    

                        fileSaveName=[num2str(id_groupPermM(zone)) 'Zone ' desFlow(1:3) '_' desImp(1:3) '_' desConc(1:3) '_' num2str(id_group1(i)) ' bit.pdf'];
                        export_fig('filename', fileSaveName, '-nocrop');
                        close(gcf);          
                    end
                end 
            end
        end
    end
end