for ref_bPerm = 1:length(r_flowSim)
    clc_nZones=r_nZones(ref_bPerm);
    clc_flowResults=[];

    clc_simFlow=r_flowSim{ref_bPerm};
    clc_flow=out_refFlow(ref_bPerm,:);

    
    for d_refFlow=1:6
        %% External - Total
        d_3DFlow=[];
        d_3DFlow1=[];
        
        d_3DFlow(:,1)=-sum(clc_simFlow);
        d_3DFlow1(:,1)=d_3DFlow;
        d_nNoise=size(clc_flow{1,d_refFlow},2);
        for d_noise=1:d_nNoise
            switch d_refFlow
                case {1,4}; d_3DFlow1(:,1+d_noise)=clc_flow{1,d_refFlow}(1,d_noise);
                case {2,5}; d_3DFlow1(:,1+d_noise)=-sum(clc_flow{1,d_refFlow}(:,d_noise));
                case {3,6}; d_3DFlow1(:,1+d_noise)=sum(clc_flow{1,d_refFlow}(:,d_noise));
            end
        end

        clc_flowResults{1}{1,1}{d_refFlow}=d_3DFlow1;
        
        %% External - Exfiltration
        if any(d_refFlow==[2 3 5 6])
            for d_zone=1:clc_nZones
                d_3DFlow=[];
                d_3DFlow1=[];

                d_3DFlow(:,1)=-sum(clc_simFlow((d_zone-1)*clc_nZones+1:d_zone*clc_nZones));
                d_3DFlow1(:,1)=d_3DFlow;
                d_nNoise=size(clc_flow{1,d_refFlow},2);
                for d_noise=1:d_nNoise
                    switch d_refFlow
                        case {2,5}; d_3DFlow1(:,1+d_noise)=-sum(clc_flow{1,d_refFlow}((d_zone-1)*clc_nZones+1:d_zone*clc_nZones,d_noise));
                        case {3,6}; d_3DFlow1(:,1+d_noise)=clc_flow{1,d_refFlow}(d_zone,d_noise);
                    end
                end
                clc_flowResults{2}{1,d_zone}{d_refFlow}=d_3DFlow1;
            end
        end 
    end
    outB_aFlowResults{ref_bPerm}=clc_flowResults;
    
    %% External - Exfiltration
%     for d_zone=1:clc_nZones
%         for d_solve=d_reqSolve
%              for d_seqA=clc_nSeqAverage
%                 if (clc_afType=='S' || clc_afType=='F')
%                     d_seqVlim=clc_nRunSeq-d_seqA;
%                 else
%                     d_seqVlim=1;
%                 end
%                 d_3DFlow=[];
%                 d_3DFlow1=[];
%                 d_3DFlow2=[];
%                 for d_imp=d_reqImp
%                     if (d_imp==3 && d_seqA>1)
%                         continue;
%                     end
%                     for d_conc=d_reqConc
%                         d_3DFlow(:,1)=clc_simFlowTime{d_seqA}(1:d_seqVlim,:);
%                         d_3DFlow(:,2)=-sum(clc_simFlow{d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones),2);
%                         d_3DFlow1{d_imp,d_conc}(:,1:2)=d_3DFlow;
%                         d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA},3);
%                         for d_noise=1:d_nNoise
%                             d_3DFlow1{d_imp,d_conc}(:,2+d_noise)=-sum(clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones,d_noise),2);
%                         end
% 
%                         if (d_seqA>1)
%                             d_3DFlow2{d_imp,d_conc}(:,1:2)=d_3DFlow;
%                             d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA},3);
%                             for d_noise=1:d_nNoise
%                                 d_3DFlow2{d_imp,d_conc}(:,2+d_noise)=-sum(clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones,d_noise),2);
%                             end
%                         end
% 
%                     end
%                 end
% 
%                 clc_flowResults{d_solve,2}{1,d_zone}{1,find(clc_nSeqAverage==d_seqA)}=d_3DFlow1;
%                 if (d_seqA>1)
%                     clc_flowResults{d_solve,2}{1,d_zone}{2,find(clc_nSeqAverage==d_seqA)}=d_3DFlow2;
%                 end
%             end
%         end
%     end
    
    
    
end