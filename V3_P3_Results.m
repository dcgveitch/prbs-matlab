% Process Simulink results

%% Read test description
clear;
%#ok<*FNDSB>
tic

[d_upperPath, d_folderTS, ~] = fileparts(pwd);

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'), '-regexp', '^(?!r_flowSim)...')

% mat_setup=matfile(strcat(d_folderTS(1:11), '_setup.mat'));
mat_outP1=matfile(strcat(d_folderTS(1:11), '__outP1.mat'),'Writable',true);
mat_outP2=matfile(strcat(d_folderTS(1:11), '__outP2.mat'),'Writable',true);
mat_outP3=matfile(strcat(d_folderTS(1:11), '__outP3.mat'),'Writable',true);

mat_outP1Info=whos(mat_outP1);
if (ismember('out_pftConc', {mat_outP1Info.name}))
    d_pft=1;
else
    d_pft=0;
end  

d_reqSolve=[1];
d_reqImp=[1 2 3];
d_reqConc=[1 2 3 4];

d_nBatch=ceil(setup_nSim/setup_batchSize);
setup_batchTrim=setup_batchSize;

for d_batch=1:d_nBatch
    d_batchL=(d_batch-1)*setup_batchSize+1;
    d_batchH=min(d_batchL+setup_batchTrim-1,setup_nSim);
    d_batchSize=d_batchH-d_batchL+1;
    
    rB_nZones=r_nZones(d_batchL:d_batchH);
    rB_nDays=r_nDays(d_batchL:d_batchH);
    rB_tZones=r_tZones(d_batchL:d_batchH);
    rB_seqLength=r_seqLength(d_batchL:d_batchH);
    rB_seqPeriod=r_seqPeriod(d_batchL:d_batchH);
    rB_seqMultiple=r_seqMultiple(d_batchL:d_batchH);
    rB_stepSize=r_stepSize(d_batchL:d_batchH);
    rB_nSeqAverage=r_nSeqAverage(d_batchL:d_batchH);
    rB_releaseRate=r_releaseRate(d_batchL:d_batchH);
    rB_releaseRateT=r_releaseRateT(d_batchL:d_batchH);
    rB_zoneVol=r_zoneVol(d_batchL:d_batchH,:);
    rB_afType=r_afType(d_batchL:d_batchH,:);
    
    rB_simFlow=mat_outP2.out_simFlow(1,d_batchL:d_batchH);
    rB_simFlowTime=mat_outP2.out_simFlowTime(1,d_batchL:d_batchH);
    rB_prbsFlow=mat_outP1.out_prbsFlow(1,d_batchL:d_batchH);
    rB_simFlowTimeFull=mat_outP2.out_simFlowTimeFull(1,d_batchL:d_batchH);
    rB_flow=mat_outP2.out_flow(1,d_batchL:d_batchH);
    if (d_pft==1)
        rB_pftConc=mat_outP1.out_pftConc(1,d_batchL:d_batchH);
        rB_pftTracer=mat_outP1.out_pftTracer(1,d_batchL:d_batchH);
    else
        rB_pftConc=[];
        rB_pftTracer=[];
    end
    
    disp(['Processing Batch ' num2str(d_batch) '/' num2str(d_nBatch)]);
    
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
        clc_prbsFlow=rB_prbsFlow{ref_bPerm};
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


        if (clc_afType=='S' || clc_afType=='F')
            clc_flowFullRef{1}{1}(:,1)=clc_simFlowTimeFull;
            clc_flowFullRef{1}{1}(:,2)=-sum(clc_prbsFlow,2);
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
            if (clc_afType=='S' || clc_afType=='F')
                clc_flowFullRef{2}{1,d_zone}(:,1)=clc_simFlowTimeFull;
                clc_flowFullRef{2}{1,d_zone}(:,2)=-sum(clc_prbsFlow(:,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones),2);
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

            if (clc_afType=='S' || clc_afType=='F')
                clc_flowFullRef{2}{2,d_zone}(:,1)=clc_simFlowTimeFull;
                clc_flowFullRef{2}{2,d_zone}(:,2)=-clc_prbsFlow(:,(d_zone-1)*clc_nZones+d_zone);
                for d_k=1:clc_nZones
                    if (d_k~=d_zone) 
                        clc_flowFullRef{2}{2,d_zone}(:,2)=clc_flowFullRef{2}{1,d_zone}(:,2)-clc_prbsFlow(:,(d_k-1)*clc_nZones+d_zone);
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

                if (clc_afType=='S' || clc_afType=='F')
                    clc_flowFullRef{3}{d_zone1,d_zone2}(:,1)=clc_simFlowTimeFull;
                    if (d_zone1==d_zone2) % Zone exfiltration
                        clc_flowFullRef{3}{d_zone1,d_zone2}(:,2)=clc_flowFullRef{2}{1,d_zone2}(:,2);
                    elseif (d_zone1==clc_nZones+1) % Zone infiltration
                        clc_flowFullRef{3}{d_zone1,d_zone2}(:,2)=clc_flowFullRef{2}{2,d_zone2}(:,2);
                    else
                        clc_flowFullRef{3}{d_zone1,d_zone2}(:,2)=clc_prbsFlow(:,(d_zone1-1)*clc_nZones+d_zone2);
                    end
                end
            end
        end

%         clc_pftFlow=[];
%         %% Calculate PFT flowrate and potential bias
%         if ((clc_afType=='S' || clc_afType=='F'))
%             clc_pftConc=rB_pftConc{ref_bPerm};
%             clc_pftTracer=rB_pftTracer{ref_bPerm};
%             clc_pftFlow(1)=-mean(sum(clc_prbsFlow(1:end,:),2));
%             clc_pftFlow(2)=sum(mean(clc_pftConc(round(setup_nDaysStab*24/clc_stepSize):end,:),1).*clc_zoneVol(1:clc_nZones))/sum(clc_zoneVol(1:clc_nZones));
%             clc_pftFlow(3)=sum(clc_pftTracer(1,:));
%             clc_pftFlow(4)=clc_pftFlow(3)/(clc_pftFlow(2)/1000000);
%             clc_pftFlow(5)=clc_pftFlow(4)/clc_pftFlow(1);
%             outB_aPftResults{ref_bPerm}=clc_pftFlow;
%         end

        outB_aFlowResults{ref_bPerm}=clc_flowResults;
%         if (clc_afType=='S' || clc_afType=='F')
%             outB_aFlowFullRef{ref_bPerm}=clc_flowFullRef;
%         end
    end

    mat_outP3.out_aFlowResults(1,d_batchL:d_batchH)=outB_aFlowResults;
%     if (d_pft==1)
%         mat_outP3.out_aFlowFullRef(1,d_batchL:d_batchH)=outB_aFlowFullRef;
%         mat_outP3.out_aPftResults(1,d_batchL:d_batchH)=outB_aPftResults;
%     end
    
    clear outB_* rB_*;
end

%% Save output
cd ..;

toc