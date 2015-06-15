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
mat_outP6=matfile(strcat(d_folderTS(1:11), '__outP6Single.mat'),'Writable',true);

d_reqTSeqA=[1 2];

ind_1_Group=unique(r_seqLength);
ind_2_Group=unique(r_seqPeriod);
ind_3_Group=unique(r_nZones);

setup_batchSize=setup_nMC;
setup_batchProc=95;
setup_batchTrim=setup_nMC;
d_batchRef=[];

d_nZonesUnique=[2 3 5 8];

for d_i=1:ceil(setup_nSim/setup_batchSize)
    for d_j=1:setup_batchTrim
        d_batchRef=[d_batchRef (d_i-1)*setup_batchSize+d_j];
    end
end

d_batchRef(d_batchRef>setup_nSim)=[];
d_count=ones(4,8,5,2);

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
    
    rB_flowResults=mat_outP3.out_aFlowResults(1,d_batchRun);

    disp(['Processing Batch ' num2str(d_batch) '/' num2str(ceil(length(d_batchRef)/setup_batchProc))]);
    d_batchList=[];
    
    for ref_bPerm = 1:d_batchSize
        disp([' -Run ' num2str(ref_bPerm) '/' num2str(d_batchSize)]);
        
        clc_nZones=rB_nZones(ref_bPerm);
        clc_seqLength=rB_seqLength(ref_bPerm);
        clc_seqPeriod=rB_seqPeriod(ref_bPerm);
        clc_flowResults=rB_flowResults{ref_bPerm};
        clc_nSeqAverage=rB_nSeqAverage{ref_bPerm};
        
        ind_1=find(ind_1_Group==clc_seqLength);
        ind_2=find(ind_2_Group==clc_seqPeriod);
        ind_3=find(ind_3_Group==clc_nZones);
        d_reqNSeqA=clc_nSeqAverage;

        for d_solve=d_reqSolve
            for d_impulse=d_reqImp
                for d_nSeqA=d_reqNSeqA
                    for d_tSeqA=d_reqTSeqA
                        for d_conc=d_reqConc
                            d_countMin=d_count;
                            for d_flowType=1:4
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

                                d_flowTotal=zeros(d_ndt,1);
                                % Total flow for airflow type
                                for d_flow1=1:d_nFlow1
                                    for d_flow2=1:d_nFlow2
                                        d_flowCount=(d_flow1-1)*d_nFlow2+d_flow2;
                                        d_flowProcessP=d_flowProcess{d_flow1,d_flow2}{d_tSeqA,d_nSeqA}{d_impulse,d_conc};
                                        d_flowTotal=d_flowTotal+abs(d_flowProcessP(:,2));
                                    end
                                end
                                out_flowTotal{d_batchRun(ref_bPerm)}(d_flowType,:)=d_flowTotal;
                                
                                d_resultsCombined=[];
                                d_resultsCombinedSum=[];
                                % Output airflow and weight
                                for d_dt=1:d_ndt
                                    for d_flow1=1:d_nFlow1
                                        for d_flow2=1:d_nFlow2
                                            d_flowCount=(d_flow1-1)*d_nFlow2+d_flow2;                        
                                            d_flowProcessP=d_flowProcess{d_flow1,d_flow2}{d_tSeqA,d_nSeqA}{d_impulse,d_conc}';
                                            d_countL=d_count(d_flowType,d_conc,d_solve,d_impulse);

                                            % Individual raw results with references
                                            out_resultsCombined{d_flowType,d_conc}{d_solve,d_impulse}(d_countL:d_countL+d_nNoise-1,1)=repmat(abs(d_flowProcessP(2,d_dt))/d_flowTotal(d_dt)/d_ndt/d_nNoise,d_nNoise,1); 
                                            out_resultsCombined{d_flowType,d_conc}{d_solve,d_impulse}(d_countL:d_countL+d_nNoise-1,2)=(d_flowProcessP(3:end,d_dt)-d_flowProcessP(2,d_dt))/d_flowProcessP(2,d_dt); % Only processing single version of flow
                                            out_resultsCombined{d_flowType,d_conc}{d_solve,d_impulse}(d_countL:d_countL+d_nNoise-1,3)=repmat(find(d_nZonesUnique==clc_nZones),d_nNoise,1);
                                            out_resultsCombined{d_flowType,d_conc}{d_solve,d_impulse}(d_countL:d_countL+d_nNoise-1,4)=repmat(d_batchRun(ref_bPerm),d_nNoise,1);
                                            out_resultsCombined{d_flowType,d_conc}{d_solve,d_impulse}(d_countL:d_countL+d_nNoise-1,5)=repmat(d_flowCount,d_nNoise,1);
                                            d_count(d_flowType,d_conc,d_solve,d_impulse)=d_countL+d_nNoise;
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
    
    clear outB_* d_in;
end

d_procTime=toc
mat_outP6.d_procTime=d_procTime;
mat_outP6.out_resultsCombined=out_resultsCombined;
mat_outP6.out_flowTotal=out_flowTotal;
cd ../..;



