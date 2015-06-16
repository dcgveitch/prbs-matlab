d_reqTSeqA=[1 2];

%% External - Total
for d_sens=1:4    
    for d_solve=d_reqSolve
        for d_tSeqA=d_reqTSeqA
            for d_nSeqA=1:clc_nSeqAverage
                for d_imp=d_reqImp
                    d_seqVlim=clc_nRunSeq-d_nSeqA;
                    clc_flowResults{d_test}{d_solve,1}{1}{d_tSeqA,d_nSeqA}{d_imp,d_sens}(:,1)=clc_simFlowTime{d_nSeqA}(1:d_seqVlim,:);
                    clc_flowResults{d_test}{d_solve,1}{1}{d_tSeqA,d_nSeqA}{d_imp,d_sens}(:,2)=-sum(clc_simFlow{d_nSeqA}(1:d_seqVlim,:),2);
                    clc_flowResults{d_test}{d_solve,1}{1}{d_tSeqA,d_nSeqA}{d_imp,d_sens}(:,3)=-sum(clc_flow{d_sens,d_solve}{d_imp,d_conc}{d_tSeqA,d_nSeqA}(1:d_seqVlim,:),2);
                end
            end
        end
    end 
end

%% External - Exfiltration
for d_zone=1:clc_nZones
    for d_sens=1:4    
        for d_solve=d_reqSolve
            for d_tSeqA=d_reqTSeqA
                for d_nSeqA=1:clc_nSeqAverage
                    for d_imp=d_reqImp
                        d_seqVlim=clc_nRunSeq-d_nSeqA;
                        clc_flowResults{d_test}{d_solve,2}{1,d_zone}{d_tSeqA,d_nSeqA}{d_imp,d_sens}(:,1)=clc_simFlowTime{d_nSeqA}(1:d_seqVlim,:);
                        clc_flowResults{d_test}{d_solve,2}{1,d_zone}{d_tSeqA,d_nSeqA}{d_imp,d_sens}(:,2)=-sum(clc_simFlow{d_nSeqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones),2);
                        clc_flowResults{d_test}{d_solve,2}{1,d_zone}{d_tSeqA,d_nSeqA}{d_imp,d_sens}(:,3)=-sum(clc_flow{d_sens,d_solve}{d_imp,d_conc}{d_tSeqA,d_nSeqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones),2);
                    end
                end
            end
        end 
    end
end
 
%% External - Infiltration
for d_zone=1:clc_nZones
    for d_sens=1:4    
        for d_solve=d_reqSolve
            for d_tSeqA=d_reqTSeqA
                for d_nSeqA=1:clc_nSeqAverage
                    for d_imp=d_reqImp
                        d_seqVlim=clc_nRunSeq-d_nSeqA;
                        clc_flowResults{d_test}{d_solve,2}{2,d_zone}{d_tSeqA,d_nSeqA}{d_imp,d_sens}(:,1)=clc_simFlowTime{d_nSeqA}(1:d_seqVlim,:);
                        d_3DFlow=-clc_simFlow{d_nSeqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+d_zone);
                        for d_k=1:clc_nZones
                            if (d_k~=d_zone) 
                                d_3DFlow=d_3DFlow-clc_simFlow{d_nSeqA}(1:d_seqVlim,(d_k-1)*clc_nZones+d_zone);
                            end
                        end
                        clc_flowResults{d_test}{d_solve,2}{2,d_zone}{d_tSeqA,d_nSeqA}{d_imp,d_sens}(:,2)=d_3DFlow;
                        d_3DFlow=-clc_flow{d_sens,d_solve}{d_imp,d_conc}{d_tSeqA,d_nSeqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+d_zone);
                        for d_k=1:clc_nZones
                            if (d_k~=d_zone) 
                                d_3DFlow=d_3DFlow-clc_flow{d_sens,d_solve}{d_imp,d_conc}{d_tSeqA,d_nSeqA}(1:d_seqVlim,(d_k-1)*clc_nZones+d_zone);
                            end
                        end
                        clc_flowResults{d_test}{d_solve,2}{2,d_zone}{d_tSeqA,d_nSeqA}{d_imp,d_sens}(:,3)=d_3DFlow;
                    end
                end
            end
        end 
    end
end

%% Internal - Total
for d_sens=1:4    
    for d_solve=d_reqSolve
        for d_tSeqA=d_reqTSeqA
            for d_nSeqA=1:clc_nSeqAverage
                for d_imp=d_reqImp
                    d_seqVlim=clc_nRunSeq-d_nSeqA;
                    clc_flowResults{d_test}{d_solve,3}{1}{d_tSeqA,d_nSeqA}{d_imp,d_sens}(:,1)=clc_simFlowTime{d_nSeqA}(1:d_seqVlim,:);
                    d_3DFlow=sum(clc_simFlow{d_nSeqA}(1:d_seqVlim,:),2);
                    for d_k=1:clc_nZones 
                        d_3DFlow=d_3DFlow-clc_simFlow{d_nSeqA}(1:d_seqVlim,(d_k-1)*clc_nZones+d_k);
                    end
                    clc_flowResults{d_test}{d_solve,3}{1}{d_tSeqA,d_nSeqA}{d_imp,d_sens}(:,2)=d_3DFlow;
                    d_3DFlow=sum(clc_flow{d_sens,d_solve}{d_imp,d_conc}{d_tSeqA,d_nSeqA}(1:d_seqVlim,:),2);
                    for d_k=1:clc_nZones 
                        d_3DFlow=d_3DFlow-clc_flow{d_sens,d_solve}{d_imp,d_conc}{d_tSeqA,d_nSeqA}(1:d_seqVlim,(d_k-1)*clc_nZones+d_k);
                    end
                    clc_flowResults{d_test}{d_solve,3}{1}{d_tSeqA,d_nSeqA}{d_imp,d_sens}(:,3)=d_3DFlow;
                end
            end
        end
    end 
end

%% Internal - Individual
for d_zone1=1:clc_nZones
    for d_zone2=1:clc_nZones-1
        if (d_zone2>=d_zone1)
            d_zone2n=d_zone2+1;
        else
            d_zone2n=d_zone2;
        end
        for d_sens=1:4    
            for d_solve=d_reqSolve
                for d_tSeqA=d_reqTSeqA
                    for d_nSeqA=1:clc_nSeqAverage
                        for d_imp=d_reqImp
                            d_seqVlim=clc_nRunSeq-d_nSeqA;
                            clc_flowResults{d_test}{d_solve,4}{d_zone1,d_zone2}{d_tSeqA,d_nSeqA}{d_imp,d_sens}(:,1)=clc_simFlowTime{d_nSeqA}(1:d_seqVlim,:);
                            clc_flowResults{d_test}{d_solve,4}{d_zone1,d_zone2}{d_tSeqA,d_nSeqA}{d_imp,d_sens}(:,2)=clc_simFlow{d_nSeqA}(1:d_seqVlim,(d_zone1-1)*clc_nZones+d_zone2n);
                            clc_flowResults{d_test}{d_solve,4}{d_zone1,d_zone2}{d_tSeqA,d_nSeqA}{d_imp,d_sens}(:,3)=clc_flow{d_sens,d_solve}{d_imp,d_conc}{d_tSeqA,d_nSeqA}(1:d_seqVlim,(d_zone1-1)*clc_nZones+d_zone2n);
                        end
                    end
                end
            end 
        end
    end
end


d_count=ones(4,8,5,2,2,clc_nSeqAverage);
for d_sens=1:4    
    for d_solve=d_reqSolve
        for d_tSeqA=d_reqTSeqA
            for d_nSeqA=1:clc_nSeqAverage
                for d_imp=d_reqImp
                    d_countMin=d_count;
                    for d_flowType=1:4
                        try
                            d_flowProcess=clc_flowResults{d_test}{d_solve,d_flowType};
                            d_flowProcess{1,1}{d_tSeqA,d_nSeqA};
                        catch
                            continue;
                        end

                        d_nFlow1=size(d_flowProcess,1);
                        d_nFlow2=size(d_flowProcess,2);
                        d_nFlows=d_nFlow1*d_nFlow2;
                        d_count2=ones(d_nFlows,1);
                        d_ndt=size(d_flowProcess{1,1}{d_tSeqA,d_nSeqA}{d_imp,d_sens},1);

                        d_flowTotal=zeros(d_ndt,1);
                        % Total flow for airflow type
                        for d_flow1=1:d_nFlow1
                            for d_flow2=1:d_nFlow2
                                d_flowCount=(d_flow1-1)*d_nFlow2+d_flow2;
                                d_flowProcessP=d_flowProcess{d_flow1,d_flow2}{d_tSeqA,d_nSeqA}{d_imp,d_sens};
                                d_flowTotal=d_flowTotal+abs(d_flowProcessP(:,2));
                            end
                        end

                        d_resultsCombined=[];
                        d_resultsCombinedSum=[];
                        % Output airflow and weight
                        for d_dt=1:d_ndt
                            for d_flow1=1:d_nFlow1
                                for d_flow2=1:d_nFlow2
                                    d_flowCount=(d_flow1-1)*d_nFlow2+d_flow2;                        
                                    d_flowProcessP=d_flowProcess{d_flow1,d_flow2}{d_tSeqA,d_nSeqA}{d_imp,d_sens}';
                                    d_countL=d_count(d_flowType,d_sens,d_solve,d_imp,d_tSeqA,d_nSeqA);

                                    % Individual raw results with references
                                    out_resultsCombined{d_test}{d_flowType,d_sens}{d_solve,d_imp}{d_tSeqA,d_nSeqA}(d_countL,1)=d_flowProcessP(3:end,d_dt);
                                    out_resultsCombined{d_test}{d_flowType,d_sens}{d_solve,d_imp}{d_tSeqA,d_nSeqA}(d_countL,2)=d_flowProcessP(2,d_dt);
                                    out_resultsCombined{d_test}{d_flowType,d_sens}{d_solve,d_imp}{d_tSeqA,d_nSeqA}(d_countL,3)=(d_flowProcessP(3:end,d_dt)-d_flowProcessP(2,d_dt))/d_flowProcessP(2,d_dt); % Only processing single version of flow
                                    out_resultsCombined{d_test}{d_flowType,d_sens}{d_solve,d_imp}{d_tSeqA,d_nSeqA}(d_countL,4)=abs(d_flowProcessP(2,d_dt))/d_flowTotal(d_dt)/d_ndt; 
                                    out_resultsCombined{d_test}{d_flowType,d_sens}{d_solve,d_imp}{d_tSeqA,d_nSeqA}(d_countL,5)=clc_nZones;
                                    out_resultsCombined{d_test}{d_flowType,d_sens}{d_solve,d_imp}{d_tSeqA,d_nSeqA}(d_countL,6)=1;
                                    out_resultsCombined{d_test}{d_flowType,d_sens}{d_solve,d_imp}{d_tSeqA,d_nSeqA}(d_countL,7)=d_flowCount;
                                    d_count(d_flowType,d_sens,d_solve,d_imp,d_tSeqA,d_nSeqA)=d_countL+1;
                                end
                            end
                        end
                    end
                end
            end
        end 
    end
end




