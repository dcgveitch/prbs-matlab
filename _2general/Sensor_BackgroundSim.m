setup_nModelZones=8;
clc_tZones=2;

clc_tSchedule=cell(1);
clc_tSchedule(1:setup_nModelZones)={zeros(1,2)};
clc_nSchedule=cell(1);
clc_nSchedule(1:setup_nModelZones)={zeros(1,2)};

clc_Qgain=cell(1,setup_nModelZones);
clc_Qgain(1:setup_nModelZones)={zeros(setup_nModelZones,setup_nModelZones^2)};

for d_i=1:clc_tZones
    for d_j=1:clc_tZones
        clc_Qgain{d_j}(d_i,((d_i-1)*clc_tZones)+d_j)=1;
    end    
end

clc_stepSize=1/3600; % 1 second

clc_zoneVol(1:2)=[2.492 5.314];
clc_releaseRateSim=zeros(1,8);

clc_extSens=sens_ext/1000000;
clc_ext=[0:clc_sensorStep:(length(clc_extSens)-1)*clc_sensorStep]';
clc_ext=[clc_ext clc_extSens];

clc_Qin=zeros(size(clc_flowRef,1),setup_nModelZones^2);
clc_Qin(:,1:size(clc_flowRef,2))=clc_flowRef;
clc_Q=[0:clc_sensorStep:(size(clc_Qin,1)-1)*clc_sensorStep]';
clc_Q=[clc_Q clc_Qin];

assignin('base','sim_tSchedule',clc_tSchedule);
assignin('base','sim_nSchedule',clc_nSchedule);
assignin('base','sim_Qgain',clc_Qgain);    

assignin('base','sim_stepSize',clc_stepSize);
assignin('base','sim_Q',clc_Q);
assignin('base','sim_zoneVolGain',clc_zoneVolGain);
assignin('base','sim_releaseRate',clc_releaseRateSim);
assignin('base','sim_external',clc_ext);

[sim_prbsT,sim_prbsX,sim_prbsConc,sim_prbsTracer,sim_prbsFlow]=sim('PRBS_ParaBack',max(clc_ext(:,1)));