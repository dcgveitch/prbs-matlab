clc_ccTarget=zeros(8,1);
clc_ccTarget(1:clc_nZones)=0.001; % Set target concentration of 1000ppm
clc_tracerMax=0.15; % Max release rate of tracer in kg/hr
clc_PIDpGain=100;
clc_PIDiGain=100;
clc_PIDdGain=0;
clc_PIDsample=clc_seqPeriod/clc_seqLength/clc_seqMultiple;
clc_PIDnoise=(0.001*0.05)^2;

assignin('base','sim_ccTarget',clc_ccTarget);
assignin('base','sim_tracerMax',clc_tracerMax);
assignin('base','sim_PIDpGain',clc_PIDpGain);
assignin('base','sim_PIDiGain',clc_PIDiGain);
assignin('base','sim_PIDdGain',clc_PIDdGain);
assignin('base','sim_PIDsample',clc_PIDsample);
assignin('base','sim_PIDnoise',clc_PIDnoise);

disp(['CC Sim ' datestr(now)]);
[sim_ccT,sim_ccX,sim_ccConc,sim_ccTracer,sim_ccFlow,sim_ccError]=sim('PRBS_ParaPID',clc_nDays*24);

clc_nRunSeq=(clc_nDays-setup_nDaysStab)*24/clc_seqPeriod;
sim_flowCC=sim_ccFlow(find(sim_ccT>=24),:);
sim_tracerCC=sim_ccTracer(find(sim_ccT>=24),1:clc_tZones);
sim_concCC=sim_ccConc(find(sim_ccT>=24),1:clc_tZones);
sim_Tcc=sim_ccT(find(sim_ccT>=setup_nDaysStab*24))-24;

clc_aggregate=2;

for d_time=1:(clc_nDays-setup_nDaysStab)*24/clc_aggregate
    for d_zone=1:clc_nZones
        d_timeSel=find(sim_Tcc>=(d_time-1)*clc_aggregate & sim_Tcc<d_time*clc_aggregate);
        out_ccFlow(d_time,1)=d_time*clc_aggregate;
        out_ccFlow(d_time,d_zone+1)=mean(sim_tracerCC(d_timeSel,d_zone))/(mean(sim_concCC(d_timeSel,d_zone)/1000000));
        d_refFlow=-sim_flowCC(d_timeSel,(d_zone-1)*clc_nZones+d_zone);
        for d_k=1:clc_nZones
            if (d_k~=d_zone) 
                d_refFlow=d_refFlow-sim_flowCC(d_timeSel,(d_k-1)*clc_nZones+d_zone);
            end
        end
        out_ccRef(d_time,1)=d_time*clc_aggregate;
        out_ccRef(d_time,d_zone+1)=mean(d_refFlow);
    end
end

disp(['Complete!']);