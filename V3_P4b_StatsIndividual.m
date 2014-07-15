% Process Simulink results

% Read test description
clear;
[d_upperPath, d_folderTS, ~] = fileparts(pwd);

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'));
load(strcat(d_folderTS(1:11), '__results.mat'));
mat_prbsConc=matfile(strcat(d_folderTS(1:11), '__prbsConc.mat'),'Writable',true);
outB_prbsConcDisc=mat_prbsConc.out_prbsConcDisc;

clear out_aFlowError*;
    

d_nSolve=[1 2 3];
d_impulseType=[1 2];
d_concType=[1];

for d_perm=1:setup_nSim
    disp(['Processing Test ' num2str(d_perm) '/' num2str(setup_nSim)]);
    d_seqAList=1:length(r_nSeqAverage{d_perm});
    d_count1=1;
    for d_solve=d_nSolve
        for d_conc=d_concType
            for d_impulse=d_impulseType
                for d_seqA=d_seqAList
                    d_errorSummary=[];
                    for d_flowType=1:3
                        switch d_flowType
                            case 1
                                d_flowProcess=out_aFlowResults{d_perm}{d_solve,d_flowType}(1,1);
                            case 2
                                d_flowProcess=out_aFlowResults{d_perm}{d_solve,d_flowType}(1,:);
                            case 3
                                d_flowProcess=out_aFlowResults{d_perm}{d_solve,d_flowType}(1:end,:);
                            case 4
                                d_flowProcess=out_aFlowResults{d_perm}{d_solve,d_flowType}(1,1);
                            case 5
                                d_flowProcess=out_aFlowResults{d_perm}{d_solve,d_flowType-1}(1,1);
                            case 6
                                d_flowProcess=out_aFlowResults{d_perm}{d_solve,d_flowType-1}(1,:);
                        end
                        d_nFlow1=size(d_flowProcess,1);
                        d_nFlow2=size(d_flowProcess,2);
                        d_nNoise=size(d_flowProcess{1,1}{1,d_seqA}{d_impulse,d_conc},2)-2;
                        d_ndt=size(d_flowProcess{1,1}{1,d_seqA}{d_impulse,d_conc},1);

                        for d_noise=1:d_nNoise
                            d_count2=1;
                            d_flowError=[];
                            for d_dt=1:d_ndt
                                d_flowTotal=0;
                                for d_flow1=1:d_nFlow1
                                    for d_flow2=1:d_nFlow2
                                        d_flowTotal=d_flowTotal+abs(d_flowProcess{d_flow1,d_flow2}{1,d_seqA}{d_impulse,d_conc}(d_dt,2));
                                    end
                                end               

                                for d_flow1=1:d_nFlow1
                                    for d_flow2=1:d_nFlow2
                                        d_flowProcessP=d_flowProcess{d_flow1,d_flow2}{1,d_seqA}{d_impulse,d_conc};
                                        d_flowError(d_count2,1)=(d_flowProcessP(d_dt,d_noise+2)-d_flowProcessP(d_dt,2))/d_flowProcessP(d_dt,2);
                                        d_flowError(d_count2,2)=abs(d_flowProcessP(d_dt,2))/d_flowTotal/d_nNoise;
                                        d_count2=d_count2+1;
                                    end
                                end
                            end

                            % Weighted statistics
                            d_errorSummary(1,(d_flowType-1)*7+1)=size(d_flowError,1);
                            d_errorSummary(1,(d_flowType-1)*7+2)=wmean(d_flowError(:,1),d_flowError(:,2));
                            d_errorSummary(1,(d_flowType-1)*7+3)=wprctile(d_flowError(:,1),5,d_flowError(:,2));
                            d_errorSummary(1,(d_flowType-1)*7+4)=wprctile(d_flowError(:,1),25,d_flowError(:,2));
                            d_errorSummary(1,(d_flowType-1)*7+5)=wprctile(d_flowError(:,1),50,d_flowError(:,2));
                            d_errorSummary(1,(d_flowType-1)*7+6)=wprctile(d_flowError(:,1),75,d_flowError(:,2));
                            d_errorSummary(1,(d_flowType-1)*7+7)=wprctile(d_flowError(:,1),95,d_flowError(:,2));

                            % Non-weighted statistics
                            d_errorSummary(2,(d_flowType-1)*7+1)=size(d_flowError,1);
                            d_errorSummary(2,(d_flowType-1)*7+2)=mean(d_flowError(:,1));
                            d_errorSummary(2,(d_flowType-1)*7+3)=prctile(d_flowError(:,1),5);
                            d_errorSummary(2,(d_flowType-1)*7+4)=prctile(d_flowError(:,1),25);
                            d_errorSummary(2,(d_flowType-1)*7+5)=prctile(d_flowError(:,1),50);
                            d_errorSummary(2,(d_flowType-1)*7+6)=prctile(d_flowError(:,1),75);
                            d_errorSummary(2,(d_flowType-1)*7+7)=prctile(d_flowError(:,1),95);    
                        end
                    end
                    out_aFlowErrorW{d_perm}(d_count1,1)=d_solve;
                    out_aFlowErrorW{d_perm}(d_count1,2)=d_conc;
                    out_aFlowErrorW{d_perm}(d_count1,3)=d_impulse;
                    out_aFlowErrorW{d_perm}(d_count1,4)=d_seqA;
                    out_aFlowErrorW{d_perm}(d_count1,5)=r_nZones(d_perm);
                    out_aFlowErrorW{d_perm}(d_count1,6)=sum(r_zoneVol(d_perm,:));
                    out_aFlowErrorW{d_perm}(d_count1,7)=r_seqLength(d_perm);
                    out_aFlowErrorW{d_perm}(d_count1,8)=r_seqPeriod(d_perm);
                    out_aFlowErrorW{d_perm}(d_count1,9:29)=d_errorSummary(1,:);

                    out_aFlowErrorNW{d_perm}(d_count1,1)=d_solve;
                    out_aFlowErrorNW{d_perm}(d_count1,2)=d_conc;
                    out_aFlowErrorNW{d_perm}(d_count1,3)=d_impulse;
                    out_aFlowErrorNW{d_perm}(d_count1,4)=d_seqA;
                    out_aFlowErrorNW{d_perm}(d_count1,5)=r_nZones(d_perm);
                    out_aFlowErrorNW{d_perm}(d_count1,6)=sum(r_zoneVol(d_perm,:));
                    out_aFlowErrorNW{d_perm}(d_count1,7)=r_seqLength(d_perm);
                    out_aFlowErrorNW{d_perm}(d_count1,8)=r_seqPeriod(d_perm);
                    out_aFlowErrorNW{d_perm}(d_count1,9:29)=d_errorSummary(2,:);
                    d_count1=d_count1+1;
                end
            end
        end
    end
end

%% Save output
save(strcat(d_folderTS(1:11), '__Results.mat'), 'out_*','-v7.3');
cd ..;

