for d_perm=1:length(out_prbsFlow)
    clc_nZones=r_nZones(d_perm);
    %% Total External flow
    d_3DFlow=[];                
    d_3DFlow(:,1)=out_simFlowTimeFull{d_perm};
    d_3DFlow(:,2)=-sum(out_prbsFlow{d_perm},2);

    out_varFlowSummary{1,d_perm}{1,1}=d_3DFlow;

    %% Zone flow - Exfiltration
    for d_zone=1:clc_nZones
        d_3DFlow=[];
        d_3DFlow(:,1)=out_simFlowTimeFull{d_perm};
        d_3DFlow(:,2)=-sum(out_prbsFlow{d_perm}(:,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones),2);

        out_varFlowSummary{2,d_perm}{1,d_zone}=d_3DFlow;
    end 

    %% Zone flow - Infiltration
    for d_zone=1:clc_nZones
        d_3DFlow=[];
        d_3DFlow(:,1)=out_simFlowTimeFull{d_perm};
        d_3DFlow(:,2)=-out_prbsFlow{d_perm}(:,(d_zone-1)*clc_nZones+d_zone);
        for d_k=1:clc_nZones
            if (d_k~=d_zone) 
                d_3DFlow(:,2)=d_3DFlow(:,2)-out_prbsFlow{d_perm}(:,(d_k-1)*clc_nZones+d_zone);
            end
        end

        out_varFlowSummary{2,d_perm}{2,d_zone}=d_3DFlow;
    end 

    %% Individual flows
    for d_zone1=1:clc_nZones+1
        for d_zone2=1:clc_nZones
            d_3DFlow=[];
            d_3DFlow(:,1)=out_simFlowTimeFull{d_perm};                
            if (d_zone1==d_zone2) % Zone exfiltration
                d_3DFlow(:,2)=out_varFlowSummary{2,d_perm}{1,d_zone2}(:,2);
            elseif (d_zone1==clc_nZones+1) % Zone infiltration
                d_3DFlow(:,2)=out_varFlowSummary{2,d_perm}{2,d_zone2}(:,2);
            else
                d_3DFlow(:,2)=out_prbsFlow{d_perm}(:,(d_zone1-1)*clc_nZones+d_zone2);
            end
            out_varFlowSummary{3,d_perm}{d_zone1,d_zone2}=d_3DFlow;
        end
    end                                
end
