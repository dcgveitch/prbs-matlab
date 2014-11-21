% Process Simulink results

%% Read test description
clear;
%#ok<*FNDSB>
tic

[d_upperPath, d_folder, ~] = fileparts(pwd);
if d_folder(2)=='_', d_folderTS=d_folder(5:15);
else d_folderTS=d_folder(1:11); end

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'), '-regexp', '^(?!r_flowSim)...')
mat_outP2=matfile(strcat(d_folderTS(1:11), '__outP2.mat'),'Writable',true);
mat_outP3=matfile(strcat(d_folderTS(1:11), '__outP3.mat'),'Writable',true);

setup_batchSize=setup_nMC;
setup_batchProc=10;
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
    rB_nDays=r_nDays(d_batchRun);
    rB_tZones=r_tZones(d_batchRun);
    rB_seqLength=r_seqLength(d_batchRun);
    rB_seqPeriod=r_seqPeriod(d_batchRun);
    rB_seqMultiple=r_seqMultiple(d_batchRun);
    rB_stepSize=r_stepSize(d_batchRun);
    rB_nSeqAverage=r_nSeqAverage(d_batchRun);
    rB_releaseRate=r_releaseRate(d_batchRun);
    rB_releaseRateT=r_releaseRateT(d_batchRun);
    rB_zoneVol=r_zoneVol(d_batchRun,:);
    rB_afType=r_afType(d_batchRun,:);
    
    rB_simFlow=mat_outP2.out_simFlow(1,d_batchRun);
    rB_simFlowTime=mat_outP2.out_simFlowTime(1,d_batchRun);
    rB_simFlowTimeFull=mat_outP2.out_simFlowTimeFull(1,d_batchRun);
    rB_flow=mat_outP2.out_flow(1,d_batchRun);
    
    disp(['Processing Batch ' num2str(d_batch) '/' num2str(ceil(length(d_batchRef)/setup_batchProc))]);
    
    parfor ref_bPerm = 1:d_batchSize
        disp([' -Run ' num2str(ref_bPerm) '/' num2str(d_batchSize)]);

        %% Assign temporary variables
        clc_nZones=rB_nZones(ref_bPerm);
        clc_nDays=rB_nDays(ref_bPerm);
        clc_tZones=rB_tZones(ref_bPerm);
        clc_seqLength=rB_seqLength(ref_bPerm);
        clc_seqPeriod=rB_seqPeriod(ref_bPerm);
        clc_seqMultiple=rB_seqMultiple(ref_bPerm);
        clc_stepSize=rB_stepSize(ref_bPerm);
        clc_nSeqAverage=rB_nSeqAverage{ref_bPerm};
        clc_releaseRate=rB_releaseRate{ref_bPerm};
        clc_releaseRateT=rB_releaseRateT{ref_bPerm};
        clc_zoneVol=rB_zoneVol(ref_bPerm,:);
        clc_afType=rB_afType(ref_bPerm,:);

        clc_cSeqLength=clc_seqLength*clc_seqMultiple;
        clc_nSeq = 24/clc_seqPeriod;
        clc_ndt = clc_nSeq * clc_seqLength * clc_seqMultiple;
        clc_dt = clc_seqPeriod*60*60/(clc_seqLength*clc_seqMultiple); % Seconds
        clc_dth = clc_seqPeriod/(clc_seqLength*clc_seqMultiple); % Hours
        clc_nRunSeq=(clc_nDays-setup_nDaysStab)*24/clc_seqPeriod;

        clc_flowResults=[];
        clc_flowFullRef=[];

        clc_simFlow=rB_simFlow{ref_bPerm};
        clc_simFlowTime=rB_simFlowTime{ref_bPerm};
%         clc_prbsFlow=rB_prbsFlow{ref_bPerm};
        clc_simFlowTimeFull=rB_simFlowTimeFull{ref_bPerm};
        clc_flow=rB_flow{ref_bPerm};

        %% Total External flow
        for d_solve=d_reqSolve
            for d_seqA=clc_nSeqAverage
                if (clc_afType=='S' || clc_afType=='F')
                    d_seqVlim=clc_nRunSeq-d_seqA;
                else
                    d_seqVlim=1;
                end
                d_3DFlow=[];
                d_3DFlow1=[];
                d_3DFlow2=[];
                for d_imp=d_reqImp
                    if (d_imp==3 && d_seqA>1)
                        continue;
                    end
                    for d_conc=d_reqConc
                        d_3DFlow(:,1)=clc_simFlowTime{d_seqA}(1:d_seqVlim,:);
                        d_3DFlow(:,2)=-sum(clc_simFlow{d_seqA}(1:d_seqVlim,:),2);
                        d_3DFlow1{d_imp,d_conc}(:,1:2)=d_3DFlow;
                        d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA},3);
                        for d_noise=1:d_nNoise
                            d_3DFlow1{d_imp,d_conc}(:,2+d_noise)=-sum(clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA}(1:d_seqVlim,:,d_noise),2);
                        end

                        if (d_seqA>1)
                            d_3DFlow2{d_imp,d_conc}(:,1:2)=d_3DFlow;
                            d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA},3);
                            for d_noise=1:d_nNoise
                                d_3DFlow2{d_imp,d_conc}(:,2+d_noise)=-sum(clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA}(1:d_seqVlim,:,d_noise),2);
                            end
                        end

                    end
                end

                clc_flowResults{d_solve,1}{1}{1,find(clc_nSeqAverage==d_seqA)}=d_3DFlow1;
                if (d_seqA>1)
                    clc_flowResults{d_solve,1}{1}{2,find(clc_nSeqAverage==d_seqA)}=d_3DFlow2;
                end
            end
        end

        %% Zone flow - Exfiltration
        for d_zone=1:clc_nZones
            for d_solve=d_reqSolve
                 for d_seqA=clc_nSeqAverage
                    if (clc_afType=='S' || clc_afType=='F')
                        d_seqVlim=clc_nRunSeq-d_seqA;
                    else
                        d_seqVlim=1;
                    end
                    d_3DFlow=[];
                    d_3DFlow1=[];
                    d_3DFlow2=[];
                    for d_imp=d_reqImp
                        if (d_imp==3 && d_seqA>1)
                            continue;
                        end
                        for d_conc=d_reqConc
                            d_3DFlow(:,1)=clc_simFlowTime{d_seqA}(1:d_seqVlim,:);
                            d_3DFlow(:,2)=-sum(clc_simFlow{d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones),2);
                            d_3DFlow1{d_imp,d_conc}(:,1:2)=d_3DFlow;
                            d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA},3);
                            for d_noise=1:d_nNoise
                                d_3DFlow1{d_imp,d_conc}(:,2+d_noise)=-sum(clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones,d_noise),2);
                            end

                            if (d_seqA>1)
                                d_3DFlow2{d_imp,d_conc}(:,1:2)=d_3DFlow;
                                d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA},3);
                                for d_noise=1:d_nNoise
                                    d_3DFlow2{d_imp,d_conc}(:,2+d_noise)=-sum(clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones,d_noise),2);
                                end
                            end

                        end
                    end

                    clc_flowResults{d_solve,2}{1,d_zone}{1,find(clc_nSeqAverage==d_seqA)}=d_3DFlow1;
                    if (d_seqA>1)
                        clc_flowResults{d_solve,2}{1,d_zone}{2,find(clc_nSeqAverage==d_seqA)}=d_3DFlow2;
                    end
                end
            end
        end

        %% Zone flow - Infiltration
        for d_zone=1:clc_nZones
            for d_solve=d_reqSolve
                 for d_seqA=clc_nSeqAverage
                    if (clc_afType=='S' || clc_afType=='F')
                        d_seqVlim=clc_nRunSeq-d_seqA;
                    else
                        d_seqVlim=1;
                    end
                    d_3DFlow=[];
                    d_3DFlow1=[];
                    d_3DFlow2=[];
                    for d_imp=d_reqImp
                        if (d_imp==3 && d_seqA>1)
                            continue;
                        end
                        for d_conc=d_reqConc
                            d_3DFlow(:,1)=clc_simFlowTime{d_seqA}(1:d_seqVlim,:);
                            d_3DFlow(:,2)=-clc_simFlow{d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+d_zone);
                            for d_k=1:clc_nZones
                                if (d_k~=d_zone) 
                                    d_3DFlow(:,2)=d_3DFlow(:,2)-clc_simFlow{d_seqA}(1:d_seqVlim,(d_k-1)*clc_nZones+d_zone);
                                end
                            end
                            d_3DFlow1{d_imp,d_conc}(:,1:2)=d_3DFlow;
                            d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA},3);
                            for d_noise=1:d_nNoise
                                d_3DFlow1{d_imp,d_conc}(:,2+d_noise)=-clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+d_zone,d_noise);
                                for d_k=1:clc_nZones
                                    if (d_k~=d_zone) 
                                        d_3DFlow1{d_imp,d_conc}(:,2+d_noise)=d_3DFlow1{d_imp,d_conc}(:,2+d_noise)-clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA}(1:d_seqVlim,(d_k-1)*clc_nZones+d_zone,d_noise);
                                    end
                                end
                            end

                            if (d_seqA>1)
                                d_3DFlow2{d_imp,d_conc}(:,1:2)=d_3DFlow;
                                d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA},3);
                                for d_noise=1:d_nNoise
                                    d_3DFlow2{d_imp,d_conc}(:,2+d_noise)=-clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+d_zone,d_noise);
                                    for d_k=1:clc_nZones
                                        if (d_k~=d_zone) 
                                            d_3DFlow2{d_imp,d_conc}(:,2+d_noise)=d_3DFlow2{d_imp,d_conc}(:,2+d_noise)-clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA}(1:d_seqVlim,(d_k-1)*clc_nZones+d_zone,d_noise);
                                        end
                                    end
                                end
                            end
                        end
                    end

                    clc_flowResults{d_solve,2}{2,d_zone}{1,find(clc_nSeqAverage==d_seqA)}=d_3DFlow1;
                    if (d_seqA>1)
                        clc_flowResults{d_solve,2}{2,d_zone}{2,find(clc_nSeqAverage==d_seqA)}=d_3DFlow2;
                    end
                end
            end
        end

        %% Individual flows
        for d_zone1=1:clc_nZones+1
            for d_zone2=1:clc_nZones
                for d_solve=d_reqSolve
                    for d_seqA=clc_nSeqAverage
                        if (clc_afType=='S' || clc_afType=='F')
                            d_seqVlim=clc_nRunSeq-d_seqA;
                        else
                            d_seqVlim=1;
                        end
                        d_3DFlow=[];
                        d_3DFlow1=[];
                        d_3DFlow2=[];
                        for d_imp=d_reqImp
                            if (d_imp==3 && d_seqA>1)
                                continue;
                            end
                            for d_conc=d_reqConc
                                d_3DFlow(:,1)=clc_simFlowTime{d_seqA}(1:d_seqVlim,:);
                                if (d_zone1==d_zone2) % Zone exfiltration
                                    d_3DFlow(:,2)=clc_flowResults{d_solve,2}{1,d_zone2}{1,find(clc_nSeqAverage==d_seqA)}{d_imp,d_conc}(:,2);
                                elseif (d_zone1==clc_nZones+1) % Zone infiltration
                                    d_3DFlow(:,2)=clc_flowResults{d_solve,2}{2,d_zone2}{1,find(clc_nSeqAverage==d_seqA)}{d_imp,d_conc}(:,2);
                                else
                                    d_3DFlow(:,2)=clc_simFlow{d_seqA}(1:d_seqVlim,(d_zone1-1)*clc_nZones+d_zone2);
                                end
                                d_3DFlow1{d_imp,d_conc}(:,1:2)=d_3DFlow;
                                d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA},3);
                                for d_noise=1:d_nNoise
                                    if (d_zone1==d_zone2) % Zone exfiltration
                                        d_3DFlow1{d_imp,d_conc}(:,2+d_noise)=clc_flowResults{d_solve,2}{1,d_zone2}{1,find(clc_nSeqAverage==d_seqA)}{d_imp,d_conc}(:,2+d_noise);
                                    elseif (d_zone1==clc_nZones+1) % Zone infiltration
                                        d_3DFlow1{d_imp,d_conc}(:,2+d_noise)=clc_flowResults{d_solve,2}{2,d_zone2}{1,find(clc_nSeqAverage==d_seqA)}{d_imp,d_conc}(:,2+d_noise);
                                    else
                                        d_3DFlow1{d_imp,d_conc}(:,2+d_noise)=clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA}(1:d_seqVlim,(d_zone1-1)*clc_nZones+d_zone2,d_noise);
                                    end
                                end

                                if (d_seqA>1)
                                    d_3DFlow2{d_imp,d_conc}(:,1:2)=d_3DFlow;
                                    d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA},3);
                                    for d_noise=1:d_nNoise
                                        if (d_zone1==d_zone2) % Zone exfiltration
                                            d_3DFlow2{d_imp,d_conc}(:,2+d_noise)=clc_flowResults{d_solve,2}{1,d_zone2}{2,find(clc_nSeqAverage==d_seqA)}{d_imp,d_conc}(:,2+d_noise);
                                        elseif (d_zone1==clc_nZones+1) % Zone infiltration
                                            d_3DFlow2{d_imp,d_conc}(:,2+d_noise)=clc_flowResults{d_solve,2}{2,d_zone2}{2,find(clc_nSeqAverage==d_seqA)}{d_imp,d_conc}(:,2+d_noise);
                                        else
                                            d_3DFlow2{d_imp,d_conc}(:,2+d_noise)=clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA}(1:d_seqVlim,(d_zone1-1)*clc_nZones+d_zone2,d_noise);
                                        end
                                    end
                                end
                            end
                        end
                        clc_flowResults{d_solve,3}{d_zone1,d_zone2}{1,find(clc_nSeqAverage==d_seqA)}=d_3DFlow1;
                        if (d_seqA>1)
                            clc_flowResults{d_solve,3}{d_zone1,d_zone2}{2,find(clc_nSeqAverage==d_seqA)}=d_3DFlow2;
                        end
                    end
                end
            end
        end
        
        outB_aFlowResults{ref_bPerm}=clc_flowResults;
    end
    
    try
        mat_outP3.out_aFlowResults(1,d_batchRun)=outB_aFlowResults;
    catch
        outB_aFlowResults{d_batchSize}=[];
        mat_outP3.out_aFlowResults(1,d_batchRun)=outB_aFlowResults;
    end
    
    clear outB_* rB_*;
end

%% Save output
cd ..;

d_procTime=toc
mat_outP3.d_procTime=d_procTime;