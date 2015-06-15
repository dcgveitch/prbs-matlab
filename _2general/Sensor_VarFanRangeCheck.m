%% Read test description
clear;
%#ok<*FNDSB>
tic

[d_upperPath, d_folder, ~] = fileparts(pwd);
if d_folder(2)=='_', d_folderTS=d_folder(5:15);
else d_folderTS=d_folder(1:11); end

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'));
load(strcat(d_folderTS(1:11), '__outP1.mat'));
cd P6;
load(strcat(d_folderTS(1:11), '__outP6Single.mat'));


fan_m=[0.004353082 0.007693463 0.008192474 0.004130734];
fan_c=[-1.090344399 0.046566342 -0.946952276 -1.003501932];
fan_ul=[1250 1100 1100 1250];
fan_ll=[300 600 600 300];

for d_i=1:length(r_flowSim)
    clc_nZones=r_nZones(d_i);
    clc_zoneVol=r_zoneVol(d_i,:);
    d_flowProcess=r_flowSim{d_i};
    for d_zone=1:clc_nZones
        d_flowProcess(:,(d_zone-1)*clc_nZones+d_zone)=-sum(d_flowProcess(:,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones),2);
    end
    clc_flowRef{d_i}=d_flowProcess;
end

for d_i=1:length(clc_flowRef)
    clc_nZones=r_nZones(d_i);
    clc_zoneVol=r_zoneVol(d_i,:);
    d_flowProcess=clc_flowRef{d_i};
    d_flowFan=[];
    for d_j=1:size(d_flowProcess,2)
        d_flowFan(:,d_j)=(d_flowProcess(:,d_j)-fan_c(d_j))/fan_m(d_j);
    end
    clc_flowRefFan{d_i}=d_flowFan;
end

for d_i=1:length(clc_flowRefFan)
    d_flowProcess=clc_flowRefFan{d_i};
    d_flowFit=[];
    for d_j=1:size(d_flowProcess,2)
        d_flowFit(:,d_j)=d_flowProcess(:,d_j)>fan_ll(d_j) & d_flowProcess(:,d_j)<fan_ul(d_j);
    end
    clc_flowFit{d_i}=d_flowFit;
end

for d_i=1:length(clc_flowFit)
    d_flowProcess=clc_flowFit{d_i};
    for d_j=1:size(d_flowProcess,2)
        clc_flowFitSummary(d_i,d_j)=sum(d_flowProcess(:,d_j))/length(d_flowProcess(:,d_j));
    end
end

clc_flowFitSummary(:,end+1)=mean(clc_flowFitSummary,2);


for d_i=1:4
    clc_comparison{d_i}(:,1)=1-clc_flowFitSummary(:,end);
    clc_comparison{d_i}(:,2)=out_resultsCombinedSummary{d_i, 6}{5, 2}(:,8);
    clc_comparison{d_i}(:,3)=out_resultsCombinedSummary{d_i, 6}{5, 2}(:,9)-out_resultsCombinedSummary{d_i, 6}{5, 2}(:,7);
    clc_comparison{d_i}(:,4)=out_resultsCombinedSummary{d_i, 6}{5, 2}(:,10)-out_resultsCombinedSummary{d_i, 6}{5, 2}(:,6);
    clc_comparison{d_i}(:,5)=out_resultsCombinedSummary{d_i, 6}{5, 2}(:,11)-out_resultsCombinedSummary{d_i, 6}{5, 2}(:,5);
end
    
    



    