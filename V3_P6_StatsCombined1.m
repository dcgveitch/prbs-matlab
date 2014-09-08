% Process Simulink results

% %% Read test description
clear;
%#ok<*FNDSB>
tic

[d_upperPath, d_folderTS, ~] = fileparts(pwd);

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'), '-regexp', '^(?!r_flowSim)...')
mat_outP3=matfile(strcat(d_folderTS(1:11), '__outP3.mat'),'Writable',true);
mat_outP6=matfile(strcat(d_folderTS(1:11), '__outP6.mat'),'Writable',true);

out_aFlowResults=mat_outP3.out_aFlowResults;

d_reqSolve=[1 4 5];
d_reqImp=[1 2];
d_reqConc=[1 2 3 4];
d_reqNSeqA=1;
d_reqTSeqA=[1];

ind_1_Group=unique(r_seqLength);
ind_2_Group=unique(r_seqPeriod);
ind_3_Group=unique(r_nZones);

setup_nSim=38;

for d_perm=1:setup_nSim
    if (rem(d_perm,10)==0)
        disp(['Processing Test ' num2str(d_perm) '/' num2str(setup_nSim)]);
    end
    ind_1=find(ind_1_Group==r_seqLength(d_perm));
    ind_2=find(ind_2_Group==r_seqPeriod(d_perm));
    ind_3=find(ind_3_Group==r_nZones(d_perm));
    
    for d_solve=d_reqSolve
        for d_impulse=d_reqImp
            for d_nSeqA=d_reqNSeqA
                for d_tSeqA=d_reqTSeqA
                    for d_conc=d_reqConc
                        for d_flowType=1:3
                            try
                                switch d_flowType
                                    case 1
                                        d_flowProcess=out_aFlowResults{d_perm}{d_solve,d_flowType}(1,1);
                                    case 2
                                        d_flowProcess=out_aFlowResults{d_perm}{d_solve,d_flowType}(2,:);
                                    case 3
                                        d_flowProcess=out_aFlowResults{d_perm}{d_solve,d_flowType}(1:end,:);
                                end
                                d_flowProcess{1,1}{d_tSeqA,d_nSeqA};
                            catch
                                continue;
                            end

                            d_nFlow1=size(d_flowProcess,1);
                            d_nFlow2=size(d_flowProcess,2);
                            d_nFlows=d_nFlow1*d_nFlow2;
                            d_count2=ones(d_nFlows,1);
                            d_nNoise=size(d_flowProcess{1,1}{d_tSeqA,d_nSeqA}{d_impulse,d_conc},2)-2;
                            d_ndt=size(d_flowProcess{1,1}{d_tSeqA,d_nSeqA}{d_impulse,d_conc},1);
                            
                            try
                                out_resultsCombined{ind_1,ind_2,ind_3,d_solve,d_impulse,d_nSeqA,d_tSeqA,d_conc,d_flowType};
                            catch
                                out_resultsCombined{ind_1,ind_2,ind_3,d_solve,d_impulse,d_nSeqA,d_tSeqA,d_conc,d_flowType}=[];
                            end
                          
                            d_flowTotal=zeros(d_ndt,1);
                            % Total flow for airflow type
                            for d_flow1=1:d_nFlow1
                                for d_flow2=1:d_nFlow2
                                    d_flowCount=(d_flow1-1)*d_nFlow2+d_flow2;
                                    d_flowProcessP=d_flowProcess{d_flow1,d_flow2}{d_tSeqA,d_nSeqA}{d_impulse,d_conc};
                                    d_flowTotal=d_flowTotal+abs(d_flowProcessP(:,2));
                                end
                            end

                            % Output airflow and weight
                            for d_dt=1:d_ndt
                                for d_flow1=1:d_nFlow1
                                    for d_flow2=1:d_nFlow2
                                        d_flowCount=(d_flow1-1)*d_nFlow2+d_flow2;                        
                                        d_flowProcessP=d_flowProcess{d_flow1,d_flow2}{d_tSeqA,d_nSeqA}{d_impulse,d_conc}';
                                        d_resultsRaw=[];

                                        % Individual raw results with references
                                        d_resultsRaw(1:d_nNoise,1)=d_perm;
                                        d_resultsRaw(1:d_nNoise,2)=d_flowCount;
                                        d_resultsRaw(1:d_nNoise,3)=d_dt;
                                        d_resultsRaw(1:d_nNoise,4)=1:d_nNoise;
                                        d_resultsRaw(1:d_nNoise,5)=d_flowProcessP(2,d_dt);
                                        d_resultsRaw(1:d_nNoise,6)=d_flowProcessP(3:end,d_dt);
                                        d_resultsRaw(1:d_nNoise,7)=(d_flowProcessP(3:end,d_dt)-d_flowProcessP(2,d_dt))/d_flowProcessP(2,d_dt);
                                        d_resultsRaw(1:d_nNoise,8)=abs(d_flowProcessP(2,d_dt))/d_flowTotal(d_dt)/d_ndt/d_nNoise;
                                        out_resultsCombined{ind_1,ind_2,ind_3,d_solve,d_impulse,d_nSeqA,d_tSeqA,d_conc,d_flowType}(end+1:end+d_nNoise,1:8)=d_resultsRaw;
                                    end
                                end
                            end   
                        end          
                    end
                end
            end 
        end
    end
end
               
mat_outP6.out_resultsCombined=out_resultsCombined;
            
cd ..;
toc