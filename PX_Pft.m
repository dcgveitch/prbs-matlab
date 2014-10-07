% Process Simulink results

%% Read test description
clear;
%#ok<*FNDSB>
tic

[d_upperPath, d_folder, ~] = fileparts(pwd);
if d_folder(2)=='_', d_folderTS=d_folder(5:15);
else d_folderTS=d_folder(1:11); end

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'))
mat_outP3=matfile(strcat(d_folderTS, '__outP3.mat'),'Writable',true);

for d_i=1:ceil(length(r_flowSim)/setup_nMC)
    for d_j=1:min(length(r_flowSim)-(d_i-1)*setup_nMC,setup_nMC)
        d_k=(d_i-1)*setup_nMC+d_j;
        d_flow=reshape(r_flowSim{d_k}(1:r_nZones(d_k)^2),r_nZones(d_k),r_nZones(d_k)).*-1;
        d_tracerT=r_releaseRateT{d_k}(1:r_nZones(d_k))'/2;
        out_pftConcZT{d_i}(d_j,1:r_nZones(d_k))=inv(d_flow)*d_tracerT*1000000; 
        d_tracer=r_releaseRate{d_k}(1:r_nZones(d_k))'/2;
        out_pftConcZ{d_i}(d_j,1:r_nZones(d_k))=inv(d_flow)*d_tracer*1000000; % Use real tracer release
        out_pftConcZ{d_i}(d_j,1:r_nZones(d_k))=out_pftConcZ{d_i}(d_j,1:r_nZones(d_k)).*(randn(1,r_nZones(d_k))*0.05+1); % Add 5% uncertainty
        
        % Calculate PFT flow error
        out_pftConc{d_i}(d_j,1)=sum(out_pftConcZT{d_i}(d_j,1:r_nZones(d_k)).*r_zoneVol(d_k,1:r_nZones(d_k)))/sum(r_zoneVol(d_k,1:r_nZones(d_k))); % Precise volume weighted
        out_pftConc{d_i}(d_j,2)=sum(out_pftConcZ{d_i}(d_j,1:r_nZones(d_k)).*r_zoneVol(d_k,1:r_nZones(d_k)))/sum(r_zoneVol(d_k,1:r_nZones(d_k))); % Precise volume weighted
        out_pftFlow{d_i}(d_j,1)=-sum(r_flowSim{d_k}(1:r_nZones(d_k)^2)); % Reference
        d_meanRatio=out_pftConcZT{d_i}(d_j,1:r_nZones(d_k))./mean(out_pftConcZT{d_i}(d_j,1:r_nZones(d_k)));
        out_pftFlow{d_i}(d_j,2)=range(d_meanRatio); % Range of zone concentrations
        out_pftFlow{d_i}(d_j,3)=sum(d_tracerT)/(out_pftConc{d_i}(d_j,1)/1000000); % Theoretical flow
        out_pftFlow{d_i}(d_j,4)=sum(d_tracer)/(out_pftConc{d_i}(d_j,2)/1000000); % Includes release and conc errors
        out_pftFlow{d_i}(d_j,5)=out_pftFlow{d_i}(d_j,3)/out_pftFlow{d_i}(d_j,1);        
        out_pftFlow{d_i}(d_j,6)=out_pftFlow{d_i}(d_j,4)/out_pftFlow{d_i}(d_j,1);        
        out_pftFlowSummary(d_k,:)=out_pftFlow{d_i}(d_j,1:6);
        out_cRMS(d_k,1)=out_pftConc{d_i}(d_j,1);
        out_cRMS(d_k,2)=sqrt(sum((out_pftConcZT{d_i}(d_j,1:r_nZones(d_k))-out_cRMS(d_k,1)).^2.*(r_zoneVol(d_k,1:r_nZones(d_k))/sum(r_zoneVol(d_k,1:r_nZones(d_k))))));
        out_cRMS(d_k,3)=out_cRMS(d_k,2)/out_cRMS(d_k,1);
        out_cRMS(d_k,4)=out_cRMS(d_k,3)*sqrt(r_nZones(d_k));
    end
end

d_procTime=toc
mat_outP3.d_procTimePFT=d_procTime;