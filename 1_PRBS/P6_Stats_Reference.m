% Process Simulink results

% %% Read test description
clear;
%#ok<*FNDSB>
tic

[d_upperPath, d_folder, ~] = fileparts(pwd);
if d_folder(2)=='_', d_folderTS=d_folder(5:15);
else d_folderTS=d_folder(1:11); end

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'), '-regexp', '^(?!r_flowSim)...')
mat_outP3=matfile(strcat(d_folderTS(1:11), '__outP3.mat'),'Writable',true);
mkdir('P6');
cd P6;
mat_outP6=matfile(strcat(d_folderTS(1:11), '__outP6Ref.mat'),'Writable',true);

d_reqTSeqA=[1 2];

ind_1_Group=unique(r_seqLength);
ind_2_Group=unique(r_seqPeriod);
ind_3_Group=unique(r_nZones);

setup_batchSize=setup_nMC;
setup_batchProc=100;
setup_batchTrim=setup_nMC;
d_batchRef=[];

for d_i=1:ceil(setup_nSim/setup_batchSize)
    for d_j=1:setup_batchTrim
        d_batchRef=[d_batchRef (d_i-1)*setup_batchSize+d_j];
    end
end

d_batchRef(d_batchRef>setup_nSim)=[];

for d_batch=1:ceil(length(d_batchRef)/setup_batchProc)
    d_batchL=(d_batch-1)*setup_batchProc+1;
    d_batchH=min(d_batchL+setup_batchProc-1,length(d_batchRef));
    d_batchRun=d_batchRef(d_batchL:d_batchH);
    d_batchSize=length(d_batchRun);
    
    rB_nZones=r_nZones(d_batchRun);
    
    rB_seqLength=r_seqLength(d_batchRun);
    rB_seqPeriod=r_seqPeriod(d_batchRun);
    rB_seqMultiple=r_seqMultiple(d_batchRun);
    rB_stepSize=r_stepSize(d_batchRun);
    rB_nSeqAverage=r_nSeqAverage(d_batchRun);
    rB_releaseRate=r_releaseRate(d_batchRun);
    rB_releaseRateT=r_releaseRateT(d_batchRun);
    rB_zoneVol=r_zoneVol(d_batchRun,:);
    rB_afType=r_afType(d_batchRun,:);
    
    rB_flowResults=mat_outP3.out_aFlowRef(1,d_batchRun);

    disp(['Processing Batch ' num2str(d_batch) '/' num2str(ceil(length(d_batchRef)/setup_batchProc))]);
    d_batchList=[];
    
    for ref_bPerm = 1:d_batchSize
        d_count=0;
        disp([' -Run ' num2str(ref_bPerm) '/' num2str(d_batchSize)]);
        
        clc_nZones=rB_nZones(ref_bPerm);
        clc_seqLength=rB_seqLength(ref_bPerm);
        clc_seqPeriod=rB_seqPeriod(ref_bPerm);
        clc_flowResults=rB_flowResults{ref_bPerm};
        clc_nSeqAverage=rB_nSeqAverage{ref_bPerm};
        
        ind_3=find(ind_3_Group==clc_nZones);
        % Overrides for Reference flows vs PRBS
        ind_1=1;
        ind_2=1;
        d_reqSolve=3; % Repurposed to cover the reference flow methods
        d_reqNSeqA=1;
        d_reqTSeqA=1;
        d_reqImp=1;
        d_reqConc=[1 7];
        d_reqFlowType={[1] [1 2 3 4] [1 2]}; % Restrict flowtypes to what is available for each reference method

        for d_solve=1:d_reqSolve
            for d_impulse=d_reqImp
                for d_nSeqA=d_reqNSeqA
                    for d_tSeqA=d_reqTSeqA
                        for d_conc=d_reqConc
                            for d_flowType=d_reqFlowType{d_solve}
                                try
                                    d_flowProcess=clc_flowResults{d_solve,d_flowType};
                                    d_flowProcess{1,1}{d_tSeqA,d_nSeqA};
                                catch
                                    continue;
                                end
                                
                                d_batchList(end+1,:)=[ind_1,ind_2,ind_3,d_solve,d_impulse,d_nSeqA,d_tSeqA,d_conc,d_flowType];

                                d_nFlow1=size(d_flowProcess,1);
                                d_nFlow2=size(d_flowProcess,2);
                                d_nFlows=d_nFlow1*d_nFlow2;
                                d_count2=ones(d_nFlows,1);
                                d_nNoise=size(d_flowProcess{1,1}{d_tSeqA,d_nSeqA}{d_impulse,d_conc},2)-2;
                                d_ndt=size(d_flowProcess{1,1}{d_tSeqA,d_nSeqA}{d_impulse,d_conc},1);

                                try
                                    outB_resultsCombined{ind_1,ind_2,ind_3,d_solve,d_impulse,d_nSeqA,d_tSeqA,d_conc,d_flowType};
                                catch
                                    outB_resultsCombined{ind_1,ind_2,ind_3,d_solve,d_impulse,d_nSeqA,d_tSeqA,d_conc,d_flowType}=[];
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

                                            % Individual raw results with references
                                            outB_resultsCombined{ind_1,ind_2,ind_3,d_solve,d_impulse,d_nSeqA,d_tSeqA,d_conc,d_flowType}{end+1,1}=[d_batchRun(ref_bPerm) clc_nZones d_flowCount d_dt d_flowProcessP(2,d_dt) d_flowProcessP(2,d_dt)/d_flowTotal(d_dt)/d_ndt/d_nNoise]; 
                                            outB_resultsCombined{ind_1,ind_2,ind_3,d_solve,d_impulse,d_nSeqA,d_tSeqA,d_conc,d_flowType}{end,2}=(d_flowProcessP(3:end,d_dt)-d_flowProcessP(2,d_dt))/d_flowProcessP(2,d_dt);
                                            d_count=d_count+d_nNoise;
                                        end
                                    end
                                end   
                            end          
                        end
                    end
                end 
            end
        end
        disp([' --Count ' num2str(d_count)]);
    end
    
    if(~isempty(d_batchList))
        d_bL=unique(d_batchList,'rows');
        for d_i=1:size(d_bL,1)
            filename=['R' num2str(d_bL(d_i,1)) '_' num2str(d_bL(d_i,2)) '_' num2str(d_bL(d_i,3)) '_' num2str(d_bL(d_i,4)) '_' num2str(d_bL(d_i,5)) '_' num2str(d_bL(d_i,6)) '_' num2str(d_bL(d_i,7)) '_' num2str(d_bL(d_i,8)) '_' num2str(d_bL(d_i,9))];
            out_results=outB_resultsCombined{d_bL(d_i,1),d_bL(d_i,2),d_bL(d_i,3),d_bL(d_i,4),d_bL(d_i,5),d_bL(d_i,6),d_bL(d_i,7),d_bL(d_i,8),d_bL(d_i,9)};
            save(filename,'out_results','-v6');
            clear out_results;
        end
    end
    
    clear outB_* d_in;
end

d_procTime=toc
mat_outP6.d_procTime=d_procTime;
cd ../..;



