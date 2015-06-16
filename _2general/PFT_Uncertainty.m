in_mcNZones=[2 3 5 8];
in_mcZoneVol=[60 52 40 30];

setup_nModelZones=8;
setup_mFlowExtLogSD = 0.5;
setup_mFlowIntLogSD = 0.5;
setup_nZones=[2 3 5 8];

r_afRefs=[0.5 1.5];

setup_nPerm=5;
setup_nMC=1000;

for ref_zone=1:4;
    for ref_perm=1:setup_nPerm
        d_sim=(ref_zone-1)*setup_nPerm+ref_perm;
        clc_nZones=in_mcNZones(ref_zone);
    
        for d_zone=1:setup_nModelZones
            if (d_zone<=clc_nZones)
                r_zoneVol(d_sim,d_zone)=0.5*in_mcZoneVol(ref_zone)+rand*in_mcZoneVol(ref_zone);
                r_zoneVolGain(d_sim,d_zone)=1/r_zoneVol(d_sim,d_zone);
            else
                r_zoneVol(d_sim,d_zone)=0;
                r_zoneVolGain(d_sim,d_zone)=0;
            end
        end

        d_rateSD=[0.02 0.02];
        for d_zone=1:r_nZones(d_sim)
            d_releaseT(d_zone)=in_mcZoneVol(ref_zone)*0.5*1000/1000000*2; % Approximately 1000ppm at 0.5ACH
        end
        d_releaseBias=d_releaseT*d_rateSD(1)*randn; % Add bias (consistent across zones)
        d_release=d_releaseT+d_releaseBias+d_releaseT*d_rateSD(2).*randn(1,length(d_releaseT)); % Add random uncertainty

        r_releaseRate{d_sim} = d_release; % Actual release rate
        r_releaseRateT{d_sim}= d_releaseT; % Theoretical release rate for calculations

        d_failFlag=1;
        d_failReselect=-1;
        while d_failFlag>0
            d_ext=random('logn',log(r_afRefs(1)),setup_mFlowExtLogSD,1,1);
            d_int=random('logn',log(r_afRefs(2)),setup_mFlowIntLogSD,1,1);

            d_totalVol=sum(r_zoneVol(d_sim,:))*d_ext;
            d_split=sort(rand(r_nZones(d_sim)-1,1));
            d_split=[0; d_split; 1];
            d_split=(d_split(2:end)-d_split(1:end-1))'.*r_zoneVol(d_sim,1:r_nZones(d_sim));
            d_split=d_split*d_totalVol/sum(d_split);
            r_flowExt(d_sim,1:r_nZones(d_sim))=d_split./r_zoneVol(d_sim,1:r_nZones(d_sim));

            d_failReselect=d_failReselect+1;
            d_failCount=-1;
            while d_failFlag>0 && d_failCount<50000
                d_failFlag=0;
                d_flowTest=[];
                d_failCount=d_failCount+1;

                % New weighted assignment of internal flows
                d_totalVol=sum(r_zoneVol(d_sim,:))*d_int;
                d_split=sort(rand(r_nZones(d_sim)*(r_nZones(d_sim)-1)-1,1));
                d_split=[0; d_split; 1];
                d_zoneWeight=[];
                for d_zone1=1:r_nZones(d_sim)
                    for d_zone2=1:r_nZones(d_sim)
                        if (d_zone1~=d_zone2)
                            d_zoneWeight(end+1)=r_zoneVol(d_sim,d_zone1)*r_zoneVol(d_sim,d_zone2);
                        end
                    end
                end
                d_split=(d_split(2:end)-d_split(1:end-1))'.*d_zoneWeight;
                d_split=d_split*d_totalVol/sum(d_split);

                d_count=1;                              
                for d_zone1=1:r_nZones(d_sim)
                    for d_zone2=1:r_nZones(d_sim)
                        if (d_zone1==d_zone2)
                            d_flowTest(d_zone1,d_zone2) = r_flowExt(d_sim,d_zone1)*r_zoneVol(d_sim,d_zone1);
                        else
                            d_flowTest(d_zone1,d_zone2) = d_split(d_count);
                            r_flowIntSplit{d_sim}(d_zone1,d_zone2)=d_split(d_count)/r_zoneVol(d_sim,d_zone1);
                            d_count=d_count+1;
                        end
                    end
                end

                for d_zone=1:r_nZones(d_sim) % Check internal flows are not driving extra exfiltration
                    if (sum(d_flowTest(d_zone,:),2)<(sum(d_flowTest(:,d_zone),1)-d_flowTest(d_zone,d_zone)))
                        d_failFlag=1;
                    end
                end
            end
        end
        r_failReselect(d_sim)=d_failReselect;
        d_flow=d_flowTest;

for d_zone=1:r_tZones(d_sim)
    d_flow(d_zone,d_zone,:)=-(max(sum(d_flow(d_zone,:,:),2),(sum(d_flow(:,d_zone,:),1)-d_flow(d_zone,d_zone,:))));
end

% Linear vector of flowrates for simulation input
for d_zone1=1:r_tZones(d_sim)
    for d_zone2=1:r_tZones(d_sim)
        d_i = (d_zone1-1)*r_tZones(d_sim) + d_zone2;
        r_flowSim{d_sim}(:,d_i) = squeeze(d_flow(d_zone1,d_zone2,:));
    end
end


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