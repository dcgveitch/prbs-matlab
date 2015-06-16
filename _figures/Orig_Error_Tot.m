% Group 1 (Sub plot) = Number of Zones
% Group 2 (Colour) = Sequence Length

id_perm=1; % Run identifier for drawing connecting lines between markers
id_1=2; % Page Break - Sequence Length
id_2=3; % Sub Plot - Sequence Period
id_series=4; % Series identifier for different colors/lines

id_x1=7; % Sub Plot X Axis (Proc) - Total External Flowrate

id_y1=24; % Sub Plot Y Axis (Proc) - Calculated PRBS LTI Total External Flowrate
id_eL1=22; % Bounded area y-axis width from point
id_eL2=23; % Bounded area y-axis width from point 
id_eU2=25; % Bounded area y-axis width from point 
id_eU1=26; % Bounded area y-axis width from point 

% Weighted mean and SD
% id_y1=12; % Sub Plot Y Axis (Proc) - Calculated PRBS LTI Total External Flowrate
% id_eL1=11; % Bounded area y-axis width from point
% id_eL2=23; % Bounded area y-axis width from point 
% id_eU2=25; % Bounded area y-axis width from point 
% id_eU1=13; % Bounded area y-axis width from point 

id_x2=8; % Sub Plot X Axis (Proc) - Total External Flowrate
id_y2=9; % Sub Plot Y Axis (Raw) - Calculated PRBS LTI Total External Flowrate
id_w2=11; % Sub Plot Y Axis (Raw) - Calculated PRBS LTI Total External Flowrate
id_d2=4;

desFlow='Total Ventilation';
desFlowFN='Tot';
desTech='Sensor (Noise + 1st Order)';
desTechFN='Sensor';

d_reqSolve=[1];
d_reqImp=[1 2 3];
d_reqConc=[1 2 3 4];

d_solve=d_reqSolve;
d_seqA=1;
d_flowType=1;
d_impulse=1;
d_conc=3;

id_group1=unique(out_indResultsProc{d_impulse,d_conc}{d_solve,d_flowType}(:,id_1));
id_group2=unique(out_indResultsProc{d_impulse,d_conc}{d_solve,d_flowType}(:,id_2));
id_groupPermM=unique(out_indResultsProc{d_impulse,d_conc}{d_solve,d_flowType}(:,id_series));

seriesColour=varycolor(length(id_groupPermM));

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
stdInfo=cell(0);

for i=1:length(id_group1)
    handaxesM=[];
    handaxesB=[];
    count=0;
    figure('Position', [100, 100, 1600, 800]);
    set(gcf,'Renderer','Painters');    
    for j=1:length(id_group2)
        clear data*;
        d_graphSelect=out_indResultsProc{d_impulse,d_conc}{d_solve,d_flowType}(:,id_1)==id_group1(i) & out_indResultsProc{d_impulse,d_conc}{d_solve,d_flowType}(:,id_2)==id_group2(j);
        d_boxSelect=out_indResultsRaw{d_impulse,d_conc}{d_solve,d_flowType}(:,id_1)==id_group1(i) & out_indResultsRaw{d_impulse,d_conc}{d_solve,d_flowType}(:,id_2)==id_group2(j);
        
        d_graph=out_indResultsProc{d_impulse,d_conc}{d_solve,d_flowType}(d_graphSelect,:);
        d_box=out_indResultsRaw{d_impulse,d_conc}{d_solve,d_flowType}(d_boxSelect,:);
        
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
        dataAll(:,5)=d_box(:,id_d2);
        dataAllBox=(dataAll(:,2)-dataAll(:,1))./dataAll(:,1);
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
            text((max(xlim)-min(xlim))*0.04+min(xlim), (max(ylim)-min(ylim))*0.125+max(ylim), [desFlow ' - ' desTech], 'FontSize', 13, 'VerticalAlignment' , 'bottom');
            text((max(xlim)-min(xlim))*0.04+min(xlim), (max(ylim)-min(ylim))*0.125+max(ylim), [num2str(id_group1(i)) ' bit sequence'], 'FontSize', 11, 'VerticalAlignment' , 'top');
        end
        
        % Boxplot
        handaxesB(i) = axes;
        dataAllScreen=dataAll;
        if (d_flowType==1)
            d_screen=dataAllScreen(:,5)/2;
        elseif (d_flowType==2)
            d_screen=dataAllScreen(:,5)*2/2;
        else
            d_screen=dataAllScreen(:,5).*(dataAllScreen(:,5)+1)/2;
        end
        dataAllScreen(dataAllScreen(:,4)<1/100,:)=[];
        boxplot(dataAllScreen(:,3), 'orientation', 'horizontal', 'width', 0.5, 'outliersize', 3.5, 'jitter', 0.5, 'symbol', 'k.');
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
    end
    fileSaveName=[desFlowFN '_' desTechFN '_' num2str(id_group1(i)) ' bit.pdf'];
    export_fig('filename', fileSaveName, '-nocrop');
    close(gcf);          
end


