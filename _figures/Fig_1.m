%% Errors in discretised impulse response.
%% Use with data from '3_1_140922T1601_Impulses'

clear d_*

nX=1;
nY=1;

sizeXMOff=0.07;
sizeYMOff=0.07;
sizeXTotal=0.9;
sizeYTotal=0.78;
sizeXPad=0.01*17/14;
sizeYPad=0.02;

sizeXM=(sizeXTotal-(nX-1)*sizeXPad)/nX;
sizeYM=(sizeYTotal-(nY-1)*sizeYPad)/nY;

colours=pmkmp(10,'CubicL');

figOut=figure;
set(figOut,'Units','centimeters');
set(figOut,'Position', [5, 5, 9, 7]);
set(figOut,'Units','pixels');

set(gcf,'Renderer','Painters');

for d_i=1:length(out_prbsCrossCorr)
    crossCorr{d_i}=out_prbsCrossCorr{d_i}(:,:,1);
end
yScale=1/max(crossCorr{2}(:,1));

handAxesM(1)=axes('color','none');
d_count2=1;
for d_i=[2 3 4 5 6 8 10]
    error{d_i}=[];
    for d_flowType=1:3
        d_data=out_aFlowResults{1, d_i}{1, d_flowType};
        d_count1=1;
        for d_flow1=1:size(d_data,1)
            for d_flow2=1:size(d_data,2)
                  d_error{d_i}(d_count1,d_flowType)=(d_data{d_flow1, d_flow2}{1, 1}{1, 1}(1,3)-d_data{d_flow1, d_flow2}{1, 1}{1, 1}(1,2))/d_data{d_flow1, d_flow2}{1, 1}{1, 1}(1,2);
                  d_count1=d_count1+1;
            end
        end
    d_errorSum(d_count2,d_flowType)=mean(nonzeros(d_error{d_i}(:,d_flowType)));
    end
    d_p(d_count2)=plot(1/length(crossCorr{d_i}(:,1)):1/length(crossCorr{d_i}(:,1)):1,crossCorr{d_i}(:,1)*yScale,'color',colours(d_i,:));
    switch d_i
        case 2
            d_multipliers{d_count2}=['x1   (' num2str(d_errorSum(d_count2,3)*100,'%0.1f') ')'];
        case 3
            d_multipliers{d_count2}=['x2   (' num2str(d_errorSum(d_count2,3)*100,'%0.1f') ')'];
        case 4
            d_multipliers{d_count2}=['x3   (' num2str(d_errorSum(d_count2,3)*100,'%0.1f') ')'];
        case 5
            d_multipliers{d_count2}=['x6   (' num2str(d_errorSum(d_count2,3)*100,'%0.1f') ')'];
        case 6
            d_multipliers{d_count2}=['x12   (' num2str(d_errorSum(d_count2,3)*100,'%0.1f') ')'];
        case 8
            d_multipliers{d_count2}=['x48   (' num2str(d_errorSum(d_count2,3)*100,'%0.1f') ')'];
        case 10
            d_multipliers{d_count2}=['x192 (' num2str(d_errorSum(d_count2,3)*100,'%0.1f') ')'];
    end
    
    hold on
    d_count2=d_count2+1;    
end
set(handAxesM(1),'Position', [sizeXMOff 1-(sizeYMOff+(sizeYM+sizeYPad)) sizeXM sizeYM], 'color', 'none');
set(gca, 'layer', 'top');

xlim([0 1]);
ylim([0 1.1]);
ylab=ylabel({'Cross-Correlation' 'Output'});
set(ylab,'Position',get(ylab,'Position')+[0.065 0 0]);
xlabel('Monitoring Time (hours)');

legendflex(handAxesM(1),d_multipliers,'xscale', 0.5, 'title', {'Sampling Multiplier (Error %)'}, 'padding', [10 10 10], 'box', 'on', 'ncol', 2);

set(gca,...
    'TickLength'    , [0 0] , ...
    'XTick'         , [0:1/6:1] , ...
    'XTickLabel'    , {'0','4','8','12','16','20','24'} , ...    
    'YTickLabel'    , [] , ...    
    'YTick'         , [0:0.25:1.25]);
grid on;

fileSaveName=['1_Discretisation.pdf'];
export_fig('filename', fileSaveName, '-nocrop');
close(gcf);
