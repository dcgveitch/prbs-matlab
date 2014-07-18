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

d_reqSolve=[1 2 3];
d_reqImp=[1 2];
d_reqConc=[1 2];

mat_outP1Info=whos(mat_outP1);
if (ismember('out_pftConc', {mat_outP1Info.name}))
    out_pftConc=mat_outP1.out_pftConc(1,1:14);
    out_pftTracer=mat_outP1.out_pftTracer(1,1:14);
end    

for ref_perm = 1:14
    disp(['Processing Test ' num2str(ref_perm) '/' num2str(setup_nSim)]);
    
    %% Assign temporary variables
    clc_nZones=r_nZones(ref_perm);
    clc_nDays=r_nDays(ref_perm);
    clc_tZones=r_tZones(ref_perm);
    clc_seqLength=r_seqLength(ref_perm);
    clc_seqPeriod=r_seqPeriod(ref_perm);
    clc_seqMultiple=r_seqMultiple(ref_perm);
    clc_nSeqAverage=r_nSeqAverage{ref_perm};
    clc_stepSize=r_stepSize(ref_perm);
    clc_releaseRate=r_releaseRate{ref_perm};
    clc_releaseRateT=r_releaseRateT{ref_perm};
    clc_zoneVol=r_zoneVol(ref_perm,:);
    clc_afType=r_afType(ref_perm,:);
    
    clc_cSeqLength=clc_seqLength*clc_seqMultiple;
    clc_nSeq = 24/clc_seqPeriod;
    clc_ndt = clc_nSeq * clc_seqLength * clc_seqMultiple;
    clc_dt = clc_seqPeriod*60*60/(clc_seqLength*clc_seqMultiple); % Seconds
    clc_dth = clc_seqPeriod/(clc_seqLength*clc_seqMultiple); % Hours
    clc_nRunSeq=(clc_nDays-setup_nDaysStab)*24/clc_seqPeriod;

    clc_flowResults=[];
    clc_flowFullRef=[];
    
    clc_simFlow=mat_outP2.out_simFlow(1,ref_perm);
    clc_simFlowTime=mat_outP2.out_simFlowTime(1,ref_perm);
    clc_simFlow=clc_simFlow{1};
    clc_simFlowTime=clc_simFlowTime{1};
    clc_prbsFlow=cell2mat(mat_outP1.out_prbsFlow(1,ref_perm));
    clc_simFlowTimeFull=cell2mat(mat_outP2.out_simFlowTimeFull(1,ref_perm));

    clc_flow=mat_outP2.out_flow(1,ref_perm);
    clc_flow=clc_flow{1};
    
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
    
    %% Calculate PFT flowrate and potential bias
    if ((clc_afType=='S' || clc_afType=='F'))
        clc_pftConc=mat_outP1.out_pftConc(1,ref_perm);
        clc_pftTracer=mat_outP1.out_pftTracer(1,ref_perm);
        clc_pftConc=clc_pftConc{1};
        clc_pftTracer=clc_pftTracer{1};
        clc_pftFlow(1)=-mean(sum(clc_prbsFlow(1:end,:),2));
        clc_pftFlow(2)=sum(mean(clc_pftConc(setup_nDaysStab*24/clc_stepSize:end,:),1).*clc_zoneVol(1:clc_nZones))/sum(clc_zoneVol(1:clc_nZones));
        clc_pftFlow(3)=sum(clc_pftTracer(1,:));
        clc_pftFlow(4)=clc_pftFlow(3)/(clc_pftFlow(2)/1000000);
        clc_pftFlow(5)=clc_pftFlow(4)/clc_pftFlow(1);
        out_aPftResults(ref_perm,:)=clc_pftFlow;
    end
        
    out_aFlowResults{ref_perm}=clc_flowResults;
    if (clc_afType=='S' || clc_afType=='F')
        out_aFlowFullRef{ref_perm}=clc_flowFullRef;
    end
end

mat_outP3.out_aFlowResults=out_aFlowResults;
if (clc_afType=='S' || clc_afType=='F')
    mat_outP3.out_aFlowFullRef=out_aFlowFullRef;
    mat_outP3.out_aPftResults=out_aPftResults;
end

%% Save output
cd ..;

toc