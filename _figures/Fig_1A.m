%% Errors in discretised impulse response.
%% Use with data from '3_1_140922T1601_Impulses'

clear d_*

nX=2;
nY=2;
sizeXcm=15;
sizeYcm=9;
topMargincm=0.8;
botMargincm=1.2;
lefMargincm=1.2;
rigMargincm=0.8;

sizeXMOff=lefMargincm/sizeXcm;
sizeYMOff=botMargincm/sizeYcm;
sizeXTotal=(sizeXcm-lefMargincm-rigMargincm)/sizeXcm;
sizeYTotal=(sizeYcm-botMargincm-topMargincm)/sizeYcm;
sizeXPad=0.3/sizeXcm;
sizeYPad=0.3/sizeYcm;
sizeXM=(sizeXTotal-(nX-1)*sizeXPad)/nX;
sizeYM=(sizeYTotal-(nY-1)*sizeYPad)/nY;

colours=pmkmp(6,'CubicL');

figOut=figure;
set(figOut,'Units','centimeters');
set(figOut,'Position', [5, 5, sizeXcm, sizeYcm]);
set(figOut,'Units','pixels');

set(gcf,'Renderer','Painters');
d_axis=1;
yScaleMax=[1 0.75];

for d_type=1:2
    for d_impulse=1:2
        % PRBS Cross Correlations
        yMax=0;
        d_plot=1;
        for d_i=[6 7 8 9 10]-(d_type-1)*5
            crossCorrInput=circshift(out_prbsCrossCorr{d_i}(:,d_impulse,1),1);
            crossCorr{d_plot}=crossCorrInput(1:round(10/(r_seqPeriod(d_i)/(r_seqLength(d_i)*r_seqMultiple(d_i))))+1);
            yMax=max(yMax,max(crossCorr{d_plot}));
            d_plot=d_plot+1;
        end

        % Need reference impulses
        crossCorr{6}=out_impulseSim{1}(1:10/r_stepSize(1),d_impulse);

        yScale=yScaleMax(d_impulse)/yMax;
        yScaleRef=[1.15 1];
        crossCorr{6}=(crossCorr{6}./max(crossCorr{6})).*yMax*yScaleRef(d_impulse);
        crossCorr{6}(1)=[];

        handAxesM(d_axis)=axes('color','none');
        d_plot=1;
        for d_i=1:5
            d_p(d_plot)=plot(-(1/(length(crossCorr{d_i})-1))/2:1/(length(crossCorr{d_i})-1):1,crossCorr{d_i}*yScale,'color',colours(d_plot+1,:));
            switch rem(d_i-1,5)+1
                case 1
                    d_legend{d_plot}=['15 bit'];
                case 2
                    d_legend{d_plot}=['31 bit'];
                case 3
                    d_legend{d_plot}=['63 bit'];
                case 4
                    d_legend{d_plot}=['127 bit'];
                case 5
                    d_legend{d_plot}=['255 bit'];
            end
            hold on
            d_plot=d_plot+1;    
        end
        d_p(d_plot)=plot((1/length(crossCorr{6}(:,1)))/2:1/length(crossCorr{6}(:,1)):1,crossCorr{6}(:,1)*yScale,':r');
        d_legend{d_plot}=['Reference'];
        
        set(handAxesM(d_axis),'Position', [sizeXMOff+rem(d_axis-1,nX)*(sizeXM+sizeXPad) (sizeYMOff+(floor((d_axis-1)/nX))*(sizeYM+sizeYPad)) sizeXM sizeYM]);

        xlim([0 6/10]);
        ylim([0 1.25]);
        
        x_title{1}=['Internal Input'];
        x_title{2}=['External Input'];
        y_title{1}=['Without Subsampling'];
        y_title{2}=['With Subsampling'];
        
        xa=handAxesM(d_axis);
        
        xa.TickLength=[0.01 0];
        xa.XTick=[0:2/10:1];
        xa.YTick=[0:0.25:1.25];
        
        if (d_type==1) 
            xa.XTickLabel={'0','2','4','6', '8'};
            xlabel('Time (hours)');
        else
            xa.XTickLabel=[];
            title(x_title{d_impulse});
        end
        
        if (d_impulse==1)   
           xa.YTickLabel={'0', '', '0.5', '', '1', ''};
           ylabel('Relative Response');
        else
           xa.YTickLabel=[]; 
           h=text(6.5/10,0.625,y_title{d_type});
           h.Rotation=90;
           h.HorizontalAlignment='center';
           h.FontWeight='bold';
           h.FontSize=7.7;
        end
	        
        if (d_type==2 & d_impulse==2)
            legendflex(d_legend,'xscale', 0.5, 'title', {'PRBS Sequence Length'}, 'padding', [10 10 10], 'box', 'on', 'ncol', 2);
        end
        
        grid on;
        xa.GridAlpha=0.85;        
        d_axis=d_axis+1;
    end
end

fileSaveName=['1_Discretisation.pdf'];
export_fig('filename', fileSaveName, '-nocrop');
close(gcf);
