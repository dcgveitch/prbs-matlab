% Process Simulink results

% %% Read test description
clear;
%#ok<*FNDSB>
tic

[d_upperPath, d_folderTS, ~] = fileparts(pwd);

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'), '-regexp', '^(?!r_flowSim)...')
mat_outP3=matfile(strcat(d_folderTS(1:11), '__outP3.mat'),'Writable',true);
mat_outP5=matfile(strcat(d_folderTS(1:11), '__outP5.mat'),'Writable',true);

out_aFlowResults=mat_outP3.out_aFlowResults;

d_reqSolve=[1 2 3];
d_reqImp=[1 2];
d_reqConc=[1 2];

d_seqA=1;


for d_impulse=d_reqImp
    for d_conc=d_reqConc
        outB_graphR=cell(length(d_reqSolve),3);
        outB_graphP=cell(length(d_reqSolve),3);
        for d_solve=d_reqSolve
            for d_flowType=1:3
                disp(['Type ' num2str(d_impulse) ':' num2str(d_conc) ':' num2str(d_flowType)]);
                d_count1=1;
                for d_perm=1:setup_nSim
                    d_flowInput=cell(1,3);
                    d_flowOutput=cell(1,3);
                    d_errorSummary=[];
                    d_count3=1;                
                    if (rem(d_perm,100)==0)
                        disp(['Processing Test ' num2str(d_perm) '/' num2str(setup_nSim)]);
                    end
                    switch d_flowType
                        case 1
                            d_flowProcess=out_aFlowResults{d_perm}{d_solve,d_flowType}(1,1);
                        case 2
                            d_flowProcess=out_aFlowResults{d_perm}{d_solve,d_flowType}(2,:);
                        case 3
                            d_flowProcess=out_aFlowResults{d_perm}{d_solve,d_flowType}(1:end,:);
                    end

                    d_nFlow1=size(d_flowProcess,1);
                    d_nFlow2=size(d_flowProcess,2);
                    d_nFlows=d_nFlow1*d_nFlow2;
                    d_count2=ones(d_nFlows,1);
                    d_nNoise=size(d_flowProcess{1,1}{1,d_seqA}{d_impulse,d_conc},2)-2;
                    d_ndt=size(d_flowProcess{1,1}{1,d_seqA}{d_impulse,d_conc},1);

                    d_graphR=zeros(d_ndt*d_nNoise*d_nFlows,11);

                    for d_dt=1:d_ndt
                        d_flowTotal=0;
                        % Total flow for airflow type
                        for d_flow1=1:d_nFlow1
                            for d_flow2=1:d_nFlow2
                                d_flowCount=(d_flow1-1)*d_nFlow2+d_flow2;
                                d_flowProcessP=d_flowProcess{d_flow1,d_flow2}{1,d_seqA}{d_impulse,d_conc};
                                d_flowTotal=d_flowTotal+abs(d_flowProcessP(d_dt,2));
                            end
                        end

                        % Input airflow and weight 
                        for d_flow1=1:d_nFlow1
                            for d_flow2=1:d_nFlow2
                                d_flowCount=(d_flow1-1)*d_nFlow2+d_flow2;
                                d_flowProcessP=d_flowProcess{d_flow1,d_flow2}{1,d_seqA}{d_impulse,d_conc};
                                d_flowInput{d_flowCount}(d_dt,1)=d_flowProcessP(d_dt,2);
                                d_flowInput{d_flowCount}(d_dt,2)=abs(d_flowProcessP(d_dt,2))/d_flowTotal/d_ndt;
                            end
                        end

                        % Output airflow and weight
                        for d_noise=1:d_nNoise
                            for d_flow1=1:d_nFlow1
                                for d_flow2=1:d_nFlow2
                                    d_flowCount=(d_flow1-1)*d_nFlow2+d_flow2;                        
                                    d_flowProcessP=d_flowProcess{d_flow1,d_flow2}{1,d_seqA}{d_impulse,d_conc};
                                    d_flowOutput{d_flowCount}(d_count2(d_flowCount),1)=d_flowProcessP(d_dt,d_noise+2);
                                    d_flowOutput{d_flowCount}(d_count2(d_flowCount),2)=(d_flowProcessP(d_dt,d_noise+2)-d_flowProcessP(d_dt,2))/d_flowProcessP(d_dt,2);
                                    d_flowOutput{d_flowCount}(d_count2(d_flowCount),3)=abs(d_flowProcessP(d_dt,2))/d_flowTotal/d_ndt/d_nNoise;
                                    d_count2(d_flowCount)=d_count2(d_flowCount)+1;

                                    % Individual raw results with references
                                    d_graphR(d_count3,1)=d_perm;
                                    d_graphR(d_count3,2)=r_seqLength(d_perm);
                                    d_graphR(d_count3,3)=r_seqPeriod(d_perm);
                                    d_graphR(d_count3,4)=r_nZones(d_perm);
                                    d_graphR(d_count3,5)=d_flowCount;
                                    d_graphR(d_count3,6)=d_dt;
                                    d_graphR(d_count3,7)=d_noise;
                                    d_graphR(d_count3,8)=d_flowProcessP(d_dt,2);
                                    d_graphR(d_count3,9)=d_flowProcessP(d_dt,d_noise+2);
                                    d_graphR(d_count3,10)=(d_flowProcessP(d_dt,d_noise+2)-d_flowProcessP(d_dt,2))/d_flowProcessP(d_dt,2);
                                    d_graphR(d_count3,11)=abs(d_flowProcessP(d_dt,2))/d_flowTotal/d_ndt/d_nNoise;
                                    d_count3=d_count3+1;
                                end
                            end
                        end   
                    end

                    outB_graphR{d_solve,d_flowType}=[outB_graphR{d_solve,d_flowType}; d_graphR];

                    for d_flowCount=1:d_nFlows
                        % Individual processed results for each flow within each case
                        outB_graphP{d_solve,d_flowType}(d_count1,1)=d_perm;
                        outB_graphP{d_solve,d_flowType}(d_count1,2)=r_seqLength(d_perm);
                        outB_graphP{d_solve,d_flowType}(d_count1,3)=r_seqPeriod(d_perm);
                        outB_graphP{d_solve,d_flowType}(d_count1,4)=r_nZones(d_perm);
                        outB_graphP{d_solve,d_flowType}(d_count1,5)=d_flowCount;
                        outB_graphP{d_solve,d_flowType}(d_count1,6)=size(d_flowOutput{d_flowCount},1);
                        outB_graphP{d_solve,d_flowType}(d_count1,7)=wmean(d_flowInput{d_flowCount}(:,1),d_flowInput{d_flowCount}(:,2));
                        outB_graphP{d_solve,d_flowType}(d_count1,9)=wmean(d_flowOutput{d_flowCount}(:,2),d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,8)=outB_graphP{d_solve,d_flowType}(d_count1,9)-wstd(d_flowOutput{d_flowCount}(:,2),d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,10)=outB_graphP{d_solve,d_flowType}(d_count1,9)+wstd(d_flowOutput{d_flowCount}(:,2),d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,12)=wmean(d_flowOutput{d_flowCount}(:,1),d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,11)=outB_graphP{d_solve,d_flowType}(d_count1,12)-wstd(d_flowOutput{d_flowCount}(:,1),d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,13)=outB_graphP{d_solve,d_flowType}(d_count1,12)+wstd(d_flowOutput{d_flowCount}(:,1),d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,14)=wprctile(d_flowOutput{d_flowCount}(:,2),5,d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,15)=wprctile(d_flowOutput{d_flowCount}(:,2),15.87,d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,16)=wprctile(d_flowOutput{d_flowCount}(:,2),25,d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,17)=wprctile(d_flowOutput{d_flowCount}(:,2),50,d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,18)=wprctile(d_flowOutput{d_flowCount}(:,2),75,d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,19)=wprctile(d_flowOutput{d_flowCount}(:,2),84.13,d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,20)=wprctile(d_flowOutput{d_flowCount}(:,2),95,d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,21)=wprctile(d_flowOutput{d_flowCount}(:,1),5,d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,22)=wprctile(d_flowOutput{d_flowCount}(:,1),15.87,d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,23)=wprctile(d_flowOutput{d_flowCount}(:,1),25,d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,24)=wprctile(d_flowOutput{d_flowCount}(:,1),50,d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,25)=wprctile(d_flowOutput{d_flowCount}(:,1),75,d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,26)=wprctile(d_flowOutput{d_flowCount}(:,1),84.13,d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,27)=wprctile(d_flowOutput{d_flowCount}(:,1),95,d_flowOutput{d_flowCount}(:,3));
                        outB_graphP{d_solve,d_flowType}(d_count1,28)=sum(d_flowInput{d_flowCount}(:,2));
                        d_count1=d_count1+1;
                    end
                end          
            end
        end
        mat_outP5.out_graphP(d_impulse,d_conc)={outB_graphP};
        mat_outP5.out_graphR(d_impulse,d_conc)={outB_graphR};
        clear outB_*;  
    end
end
                
cd ..;

