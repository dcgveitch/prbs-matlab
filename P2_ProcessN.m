% Process Simulink results

%% Read test description
clear;
%#ok<*FNDSB>
tic

[d_upperPath, d_folderTS, ~] = fileparts(pwd);

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'), '-regexp', '^(?!r_flowSim)...')
mat_outP1=matfile(strcat(d_folderTS(1:11), '__outP1.mat'),'Writable',true);
mat_outP2=matfile(strcat(d_folderTS(1:11), '__outP2.mat'),'Writable',true);

d_reqSolve=[1 4 5];
d_reqImp=[1 2];
d_reqConc=[1 5 6 7 8];

save(strcat(d_folderTS(1:11), '_setup.mat'), 'd_req*', '-append')

setup_batchSize=10;
setup_batchProc=35;
setup_batchTrim=10;
d_batchRef=[];

for d_i=1:ceil(setup_nSim/setup_batchSize)
    for d_j=1:setup_batchTrim
        d_batchRef=[d_batchRef (d_i-1)*setup_batchSize+d_j];
    end
end

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
    rB_zoneVolGain=r_zoneVolGain(d_batchRun,:);
    rB_mixModel=r_mixModel(d_batchRun);
    rB_sensorSpecType=r_sensorSpecType(d_batchRun);
    rB_sensorSpecRefs=r_sensorSpecRefs(d_batchRun);
    rB_sensorResp=r_sensorResp(d_batchRun);
    rB_noiseRefsStr=r_noiseRefsStr(d_batchRun);
    rB_afType=r_afType(d_batchRun,:);
    rB_randSeeds=r_randSeeds(d_batchRun,:);
    
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
        clc_zoneVolGain=rB_zoneVolGain(ref_bPerm,:);
        clc_mixModel=rB_mixModel(ref_bPerm);
        clc_sensorSpecType=rB_sensorSpecType(ref_bPerm);
        clc_sensorSpecRefs=rB_sensorSpecRefs{ref_bPerm};
        clc_sensorResp=rB_sensorResp(ref_bPerm);
        clc_noiseRefsStr=rB_noiseRefsStr{ref_bPerm};
        clc_afType=rB_afType(ref_bPerm,:);
        clc_randSeeds=rB_randSeeds(ref_bPerm,:);

        clc_cSeqLength=clc_seqLength*clc_seqMultiple;
        clc_dSeqLength=floor(clc_seqLength/clc_nZones)*clc_seqMultiple;
        clc_nSeq = 24/clc_seqPeriod;
        clc_ndt = clc_nSeq * clc_seqLength * clc_seqMultiple;
        clc_dt = clc_seqPeriod*60*60/(clc_seqLength*clc_seqMultiple); % Seconds
        clc_dth = clc_seqPeriod/(clc_seqLength*clc_seqMultiple); % Hours
        clc_nRunSeq=(clc_nDays-setup_nDaysStab)*24/clc_seqPeriod;
        
        sim_inputFull=cell2mat(mat_outP1.out_input(1,d_batchRun(ref_bPerm)));
        sim_prbsConc=cell2mat(mat_outP1.out_prbsConc(1,d_batchRun(ref_bPerm)));
        sim_prbsFlow=cell2mat(mat_outP1.out_prbsFlow(1,d_batchRun(ref_bPerm)));
        
        clc_flow=[];
        clc_crossCorr=[];
        clc_crossCorrD=[];
        
        %% Form full set of input traces from theoretical            
        d_a=clc_stepSize/clc_sensorResp;
        if d_a~=inf
            sim_prbsSensConc=filter(d_a, [1 d_a-1], sim_prbsConc);
        else
            sim_prbsSensConc=sim_prbsConc;
        end

        sim_conc=[];
        sim_impulse=[];
        clc_sensorSpec=[];

        d_count=1;
        for d_i=1:clc_cSeqLength*clc_nRunSeq;
            while round(clc_stepSize*(d_count-1)*10^7) < round((setup_nDaysStab*24+(d_i*clc_dth))*10^7);
                d_count=d_count+1;
            end
            sim_conc{1}(d_i,1:clc_tZones)=sim_prbsConc(d_count,1:clc_tZones);
            sim_conc{5}(d_i,1:clc_tZones)=sim_prbsSensConc(d_count,1:clc_tZones);
        end

        sim_prbsConc=[];
        sim_prbsSensConc=[];
        outB_prbsConcDisc{ref_bPerm}{1}=sim_conc{1};
        outB_prbsConcDisc{ref_bPerm}{2}=sim_conc{2};

        % Sensor Spec MC
        if clc_randSeeds(4)==0
            rng('shuffle');
        else
            rng(clc_randSeeds(4));
        end

        if (clc_sensorSpecType==1)
            for d_zone=1:clc_nZones
                for d_noiseR=1:in_noiseAveNum
                    clc_sensorSpec(d_zone,d_noiseR,1) = clc_sensorSpecRefs(1)+clc_sensorSpecRefs(2)*randn(1);
                    clc_sensorSpec(d_zone,d_noiseR,2) = clc_sensorSpecRefs(3)+clc_sensorSpecRefs(4)*randn(1);
                    clc_sensorSpec(d_zone,d_noiseR,3) = clc_sensorSpecRefs(5)+clc_sensorSpecRefs(6)*randn(1);
                end
            end
        else
            for d_noiseR=1:in_noiseAveNum
                d_count=1;
                for d_zone=1:clc_nZones
                    for d_i=1:3
                        clc_sensorSpec(d_zone,d_noiseR,d_i) = clc_sensorSpecRefs(d_count);
                        d_count=d_count+1;
                    end
                end
            end
        end

        if clc_randSeeds(5)==0
            rng('shuffle');
        else
            rng(clc_randSeeds(5));
        end

        for d_i=1:clc_nZones
            for d_noise=1:in_noiseAveNum
                d_concRand=randn(length(sim_conc{1}(:,d_i)),1);
                sim_conc{2}(:,d_i,d_noise)=sim_conc{1}(:,d_i)*clc_sensorSpec(d_i,d_noise,1); % Span only
                sim_conc{3}(:,d_i,d_noise)=sim_conc{1}(:,d_i)+clc_sensorSpec(d_i,d_noise,2); % Offset only                
                sim_conc{4}(:,d_i,d_noise)=sim_conc{1}(:,d_i)+(clc_sensorSpec(d_i,d_noise,3)*(1+sim_conc{1}(:,d_i)./setup_sensRNSDrange).*d_concRand); % Random only
                sim_conc{6}(:,d_i,d_noise)=((sim_conc{5}(:,d_i)*clc_sensorSpec(d_i,d_noise,1))+clc_sensorSpec(d_i,d_noise,2)+(clc_sensorSpec(d_i,d_noise,3)*(1+sim_conc{5}(:,d_i)./setup_sensRNSDrange).*d_concRand));
                sim_conc{7}(:,d_i,d_noise)=((sim_conc{1}(:,d_i)*clc_sensorSpec(d_i,d_noise,1))+clc_sensorSpec(d_i,d_noise,2)+(clc_sensorSpec(d_i,d_noise,3)*(1+sim_conc{1}(:,d_i)./setup_sensRNSDrange).*d_concRand));
            end
        end  

        if (ismember(8,d_reqConc))
            d_a=(clc_dt/3600)/(clc_sensorResp);
            if d_a~=inf
                sim_conc{8}=filter([1 d_a-1], d_a, sim_conc{6});
            else
                sim_conc{8}=sim_conc{6};
            end
        end

        % Direct impulse traces
        if (ismember(3,d_reqImp))
            sim_impulseRaw=cell2mat(mat_outP1.out_impulseSim(1,d_batchRun(ref_bPerm)))/1000;
            d_a=clc_stepSize/clc_sensorResp;
            if d_a~=inf
                sim_impulseSens=filter(d_a, [1 d_a-1], sim_impulseRaw);
            else
                sim_impulseSens=sim_impulseRaw;
            end

            d_count=1;
            for d_i=1:clc_cSeqLength
                while round(clc_stepSize*(d_count-1)*10^7) < round((d_i*clc_dth)*10^7);
                    d_count=d_count+1;
                end
                sim_impulse{1}(d_i,1:clc_tZones^2)=sim_impulseRaw(d_count,1:clc_tZones^2);
                sim_impulse{2}(d_i,1:clc_tZones^2)=sim_impulseSens(d_count,1:clc_tZones^2);
            end

            sim_impulseRaw=[];
            sim_impulseSens=[];

            for d_i=1:clc_nZones
                for d_j=1:clc_nZones
                    d_count=(d_i-1)*clc_nZones+d_j;
                    clc_crossCorrD{3,1}(:,d_i,d_j)=sim_impulse{1}(:,d_count);
                    clc_crossCorrD{3,5}(:,d_i,d_j)=sim_impulse{2}(:,d_count);
                    for d_noise=1:in_noiseAveNum
                        d_concRand=randn(length(sim_impulse{1}(:,d_count)),1);
                        clc_crossCorrD{3,2}(:,d_i,d_j,d_noise)=sim_impulse{1}(:,d_count)*clc_sensorSpec(d_j,d_noise,1); % Span only
                        clc_crossCorrD{3,3}(:,d_i,d_j,d_noise)=sim_impulse{1}(:,d_count)+clc_sensorSpec(d_j,d_noise,2); % Offset only
                        clc_crossCorrD{3,4}(:,d_i,d_j,d_noise)=sim_impulse{1}(:,d_count)+(clc_sensorSpec(d_j,d_noise,3)*(1+sim_impulse{1}(:,d_count)./setup_sensRNSDrange).*d_concRand); % Random only
                        clc_crossCorrD{3,6}(:,d_i,d_j,d_noise)=((sim_impulse{2}(:,d_count)*clc_sensorSpec(d_j,d_noise,1))+clc_sensorSpec(d_j,d_noise,2)+(clc_sensorSpec(d_j,d_noise,3)*(1+sim_impulse{2}(:,d_count)./setup_sensRNSDrange).*d_concRand));
                        clc_crossCorrD{3,7}(:,d_i,d_j,d_noise)=((sim_impulse{1}(:,d_count)*clc_sensorSpec(d_j,d_noise,1))+clc_sensorSpec(d_j,d_noise,2)+(clc_sensorSpec(d_j,d_noise,3)*(1+sim_impulse{1}(:,d_count)./setup_sensRNSDrange).*d_concRand));
                    end
                end
            end

            if (ismember(8,d_reqConc))
                d_a=(clc_dt/3600)/(clc_sensorResp);
                if d_a~=inf
                    for d_i=1:clc_nZones
                        for d_j=1:clc_nZones
                            d_count=(d_i-1)*clc_nZones+d_j;
                            for d_noise=1:in_noiseAveNum
                                clc_crossCorrD{3,8}(:,d_i,d_j,d_noise)=filter([1 d_a-1], d_a, clc_crossCorrD{3,6}(:,d_i,d_j,d_noise));
                            end
                        end
                    end
                else
                    clc_crossCorrD{3,8}=clc_crossCorrD{3,6};
                end
            end
        end

        %% Reference Information
        clc_simFlowTime=[];
        clc_simFlow=[];
        
        for d_seqA=clc_nSeqAverage
            %%%
            if (clc_afType=='S' || clc_afType=='F')
                d_seqVlim=clc_nRunSeq-d_seqA;
            else
                d_seqVlim=1;
            end
            clc_simFlowTime{d_seqA}=[((d_seqA+1)*clc_seqPeriod)/2:clc_seqPeriod:(((d_seqA+1)*clc_seqPeriod)/2)+(d_seqVlim-1)*clc_seqPeriod]';
            for d_seqV=1:d_seqVlim
                %%%
                clc_simFlow{d_seqA}(d_seqV,:)=mean(sim_prbsFlow(((d_seqV-1)*clc_cSeqLength)+1:(d_seqV+d_seqA)*clc_cSeqLength,:));
            end
        end
        outB_simFlowTimeFull{ref_bPerm}=[0:clc_dth:(size(sim_prbsFlow,1)-1)*clc_dth]';
        outB_simFlowTime{ref_bPerm}=clc_simFlowTime;
        outB_simFlow{ref_bPerm}=clc_simFlow;        

        %% Prepare inputs for cross correlations
        % Cross-correlation input variables
        clc_crossCorrInput=[];    
        clc_a2=[];
        clc_seqPosA=[];
        clc_seqPosAi=[];

        % Static cross correlation variables
        for d_i=1:clc_nZones
            clc_crossCorrInput(1:((max(clc_nSeqAverage)+1)*clc_cSeqLength),d_i) = sim_inputFull(1:((max(clc_nSeqAverage)+1)*clc_cSeqLength),d_i)*clc_releaseRateT(d_i)-(clc_releaseRateT(d_i)/2);
            clc_a2(d_i)=(clc_releaseRateT(d_i)/2)^2;
            clc_seqPosA(d_i)=round(d_i*clc_seqLength/clc_nZones)*clc_seqMultiple;
        end

        % Adjustment for varying space between sequence repeats in zones
        clc_seqPosA=[0 clc_seqPosA max(clc_seqPosA)+clc_seqPosA];

        for d_i=1:clc_nZones
            d_temp=clc_seqPosA-clc_seqPosA(d_i);
            clc_seqPosAi(d_i,1)=1;
            for d_j=1:clc_nZones
                clc_seqPosAi(d_i,d_j+1)=d_temp(d_i+d_j)+1;
            end
        end

        clc_crossCorrInputMat=[];
        cll_crossCorrSum=[];
        cll_crossCorrCalc=[];
        cll_dt=[];

        for d_relZone=1:clc_nZones
            clc_crossCorrInputMat{d_relZone}=full(spdiags(repmat(clc_crossCorrInput(1:clc_cSeqLength,d_relZone)'/clc_cSeqLength,clc_cSeqLength,1),0:-1:-clc_cSeqLength+1,2*clc_cSeqLength,clc_cSeqLength));
        end

        %% Calculate flows
        for d_conc=d_reqConc
            %%%
            for d_noise=1:size(sim_conc{d_conc},3)
                % disp([' --Run ' num2str(ref_bPerm) ' Noise ' num2str(d_noise) '/' num2str(size(sim_conc{d_conc},3))]);
                %%%
                if (clc_afType=='S' || clc_afType=='F')
                    d_seqVlim=clc_nRunSeq-1;
                else
                    d_seqVlim=max(clc_nSeqAverage);
                end
                clc_crossCorr=[];
                % Calculate cross correlations and flows for all sequence positions required
                d_seqA=1;
                for d_seqV=1:d_seqVlim
                    %%%
                    if (ismember(1,d_reqImp) || ismember(2,d_reqImp))
                        for d_relZone=1:clc_nZones
                            for d_concZone=1:clc_nZones
                                clc_crossCorrOutput = sim_conc{d_conc}(((d_seqV-1)*clc_cSeqLength)+1:(d_seqV+1)*clc_cSeqLength,d_concZone,d_noise)/1000000;
                                clc_concMean = mean(clc_crossCorrOutput);
                                clc_gain = clc_concMean/(clc_releaseRateT(d_relZone)/2);
                                clc_crossCorrOutput = clc_crossCorrOutput-clc_concMean;
                                clc_crossCorrOutputMat=clc_crossCorrOutput';                  
                                d_crossCorr=clc_crossCorrOutputMat*clc_crossCorrInputMat{d_relZone};
                                clc_crossCorr{1}(:,d_relZone,d_concZone,d_seqV)=((clc_seqLength*d_crossCorr)+(clc_a2(d_relZone)*clc_gain))/((clc_seqLength+1)*clc_dth*clc_seqMultiple*clc_a2(d_relZone));
                            end
                        end
                    end
                    

                    if (ismember(2,d_reqImp))
                        for d_relZone=1:clc_nZones
                            for d_concZone=1:clc_nZones
                                d_offset=circshift(clc_nZones:-1:1,[0 d_relZone]);
                                clc_crossCorr{2}(:,d_relZone,d_concZone,d_seqV)=zeros(clc_dSeqLength,1);
                                for d_i=1:clc_nZones
                                    clc_crossCorr{2}(:,d_relZone,d_concZone,d_seqV)=clc_crossCorr{2}(:,d_relZone,d_concZone,d_seqV)+clc_crossCorr{1}(clc_seqPosAi(d_offset(d_i),d_i):clc_seqPosAi(d_offset(d_i),d_i)+clc_dSeqLength-1,d_offset(d_i),d_concZone,d_seqV);
                                end
                                clc_crossCorr{2}(:,d_relZone,d_concZone,d_seqV)=clc_crossCorr{2}(:,d_relZone,d_concZone,d_seqV)./clc_nZones;
                            end
                        end
                    end

                    % Direct impulses
                    if (ismember(3,d_reqImp) && (d_seqV<2))
                        clc_crossCorr{3}=clc_crossCorrD{3,d_conc}(:,:,:,d_noise);
                    end

                    for d_imp=d_reqImp
                        if (d_imp==3 && d_seqV>1)
                            continue;
                        end
                        cll_crossCorrSum=[];
                        cll_D=[];
                        cll_E=[];
                        cll_X=[];
                        cll_Xach=[];
                        cll_Xflow=[];
                        cll_dt=clc_dth;

                        % Calculate flowrates
                        if (d_imp==3)
                            cll_crossCorrCalc=clc_crossCorr{d_imp}(1:(floor(clc_seqLength/clc_nZones)*clc_seqMultiple),:,:,d_seqV);
                            cll_sumStart = 1;
                            cll_sumEnd = size(cll_crossCorrCalc,1);
                        else
                            cll_crossCorrCalc=clc_crossCorr{d_imp}(1:(floor(clc_seqLength/clc_nZones)*clc_seqMultiple),:,:,d_seqV);
                            cll_sumStart = clc_seqMultiple+1;
                            cll_sumEnd = size(cll_crossCorrCalc,1)-clc_seqMultiple+1;
                        end
                        cll_crossCorrCalc=cll_crossCorrCalc(cll_sumStart:cll_sumEnd,:,:);

                        for d_relZone=1:clc_nZones
                            for d_concZone=1:clc_nZones
                                cll_crossCorrSum(:,d_relZone,d_concZone,1)=cll_crossCorrCalc(2:end,d_relZone,d_concZone)-cll_crossCorrCalc(1:end-1,d_relZone,d_concZone);
                                cll_crossCorrSum(:,d_relZone,d_concZone,2)=cll_crossCorrCalc(2:end,d_relZone,d_concZone)+cll_crossCorrCalc(1:end-1,d_relZone,d_concZone);
                            end
                        end

                        for d_zone=1:clc_nZones
                            for d_i=1:clc_nZones
                                for d_j=1:clc_nZones
                                    cll_D(d_i,d_j,d_zone)=sum(cll_crossCorrSum(:,d_i,d_j,2).*cll_crossCorrSum(:,d_i,d_zone,2));
                                end
                                cll_E(d_i,d_zone)=sum(cll_crossCorrSum(:,d_i,d_zone,1).*cll_crossCorrSum(:,d_i,d_zone,2));
                            end
                            cll_X(:,d_zone)=cll_D(:,:,d_zone)\cll_E(:,d_zone);
                            cll_Xach(:,d_zone)=(cll_D(:,:,d_zone)\cll_E(:,d_zone))*2/cll_dt;
                            cll_Xflow(:,d_zone)=(cll_D(:,:,d_zone)\cll_E(:,d_zone))*2/cll_dt*clc_zoneVol(d_zone);
                        end

                        for d_i=1:clc_nZones
                            clc_flow{1}{d_imp,d_conc}{1,d_seqA}(d_seqV,(d_i-1)*clc_nZones+1:d_i*clc_nZones,d_noise)=cll_Xflow(d_i,:);    
                        end

                        if (ismember(2,d_reqSolve) || ismember(3,d_reqSolve) || ismember(4,d_reqSolve) || ismember(5,d_reqSolve))
                            % Prepare inputs for additional solvers
                            if (d_imp==3)
                                cll_crossCorrCalc=reshape(permute(clc_crossCorr{d_imp}(1:floor(clc_seqLength/clc_nZones)*clc_seqMultiple,:,:,d_seqV),[1 3 2]),floor(clc_seqLength/clc_nZones)*clc_seqMultiple,clc_nZones^2);
                                cll_sumStart = 1;
                                cll_sumEnd = size(cll_crossCorrCalc,1);
                            else
                                cll_crossCorrCalc=reshape(permute(clc_crossCorr{d_imp}(1:floor(clc_seqLength/clc_nZones)*clc_seqMultiple,:,:,d_seqV),[1 3 2]),floor(clc_seqLength/clc_nZones)*clc_seqMultiple,clc_nZones^2);
                                cll_sumStart = clc_seqMultiple+1;
                                cll_sumEnd = size(cll_crossCorrCalc,1)-clc_seqMultiple+1;
                            end
                            cll_crossCorrCalc=cll_crossCorrCalc(cll_sumStart:cll_sumEnd,:);

                            clc_impAve = (cll_crossCorrCalc(1:end-1,:)+cll_crossCorrCalc(2:end,:))*1000/2;
                            clc_impdt = (cll_crossCorrCalc(2:end,:)-cll_crossCorrCalc(1:end-1,:))*1000/clc_dth;

                            if (ismember(2,d_reqSolve))
                                % Non-linear least squares - SPLIT EQUATIONS      
                                d_flowSort=[];
                                for d_zone=1:clc_nZones
                                    x0=ones(clc_nZones,1)*100;
                                    lb=zeros(clc_nZones,1);
                                    ub=ones(clc_nZones,1)*500;
                                    f=@(x)vecObj_Split(x,d_zone,clc_nZones,clc_zoneVol,clc_impAve,clc_impdt);
                                    opts = optimoptions(@lsqnonlin,'Display', 'off');
                                    d_flow=lsqnonlin(f,x0,lb,ub,opts)';
                                    d_flow(d_zone)=-d_flow(d_zone);
                                    d_flowSort(d_zone,:)=d_flow;
                                end
                                clc_flow{2}{d_imp,d_conc}{1,d_seqA}(d_seqV,:,d_noise)=reshape(d_flowSort,clc_nZones^2,1);
                            end
                            
                            if (ismember(4,d_reqSolve))
                                % Linear least squares - SPLIT EQUATIONS      
                                d_flowSort=[];
                                for d_zone=1:clc_nZones
                                    % Prepare matrices
                                    Ccell=cell(clc_nZones,clc_nZones);
                                    Ccell(:,:)={zeros(length(clc_impAve(:,1)),1)};
                                    Dcell=cell(clc_nZones,1);
                                    Dcell(:,:)={zeros(length(clc_impAve(:,1)),1)};

                                    for i=1:clc_nZones
                                        Ccell{i,d_zone}=-clc_impAve(:,(i-1)*clc_nZones+d_zone)/clc_zoneVol(d_zone);
                                        for k=1:clc_nZones
                                            if (k~=d_zone)
                                                Ccell{i,k}=clc_impAve(:,(i-1)*clc_nZones+k)/clc_zoneVol(d_zone);
                                            end
                                        end
                                        Dcell{i,1}=clc_impdt(:,(i-1)*clc_nZones+d_zone);
                                    end

                                    C=cell2mat(Ccell);
                                    D=cell2mat(Dcell);   
                                    opts = optimset('Display', 'off');
                                    d_flow=lsqnonneg(C,D,opts)';
                                    d_flow(d_zone)=-d_flow(d_zone);
                                    d_flowSort(d_zone,:)=d_flow;
                                end
                                clc_flow{4}{d_imp,d_conc}{1,d_seqA}(d_seqV,:,d_noise)=reshape(d_flowSort,clc_nZones^2,1);
                            end

                            if (ismember(3,d_reqSolve))
                                % Non-linear least squares - SINGLE EQUATION
                                d_flowSort=[];
                                x0=ones(clc_nZones*(clc_nZones+1),1)*100;
                                lb=zeros(clc_nZones*(clc_nZones+1),1);
                                ub=ones(clc_nZones*(clc_nZones+1),1)*500;
                                f=@(x)vecObj(x,clc_nZones,clc_zoneVol,clc_impAve,clc_impdt);
                                opts = optimoptions(@lsqnonlin,'Display', 'off');
                                d_flow=lsqnonlin(f,x0,lb,ub,opts)';
                                d_flowSort=d_flow(clc_nZones+1:end);
                                for d_i=1:clc_nZones
                                    d_flowSort((d_i-1)*clc_nZones+d_i)=-sum(d_flowSort((d_i-1)*clc_nZones+1:d_i*clc_nZones));
                                end
                                clc_flow{3}{d_imp,d_conc}{1,d_seqA}(d_seqV,:,d_noise)=d_flowSort;
                            end
                            
                            if (ismember(5,d_reqSolve))
                                % Linear least squares - SINGLE EQUATION
                                d_flowSort=[];
                                
                                % Prepare matrices
                                Ccell=cell(clc_nZones^2,clc_nZones^2);
                                Ccell(:,:)={zeros(length(clc_impAve(:,1)),1)};
                                Dcell=cell(clc_nZones^2,1);
                                Dcell(:,:)={zeros(length(clc_impAve(:,1)),1)};

                                for i=1:clc_nZones
                                    for j=1:clc_nZones
                                            Ccell((i-1)*clc_nZones+j,(i-1)*clc_nZones+1:i*clc_nZones)={-clc_impAve(:,(j-1)*clc_nZones+i)/clc_zoneVol(i)};
                                            Dcell{(i-1)*clc_nZones+j,1}=clc_impdt(:,(j-1)*clc_nZones+i);
                                            for k=1:clc_nZones
                                                if (k~=i)
                                                    Ccell{(i-1)*clc_nZones+j,(k-1)*clc_nZones+i}=clc_impAve(:,(j-1)*clc_nZones+k)/clc_zoneVol(i);
                                                end
                                            end
                                    end
                                end

                                C=cell2mat(Ccell);
                                D=cell2mat(Dcell);

                                C=[zeros(size(C,1),clc_nZones) C];
                                D=[D; zeros(clc_nZones,1)];

                                Cext=zeros(clc_nZones,clc_nZones*(clc_nZones+1));

                                for i=1:clc_nZones
                                    Cext(i,i)=-1;
                                    Cext(i,i*clc_nZones+1:(i+1)*clc_nZones)=1;
                                    for j=1:clc_nZones
                                        if (j~=i)
                                            Cext(i,j*clc_nZones+i)=-1;
                                        end
                                    end
                                end

                                C=[C; Cext];      
                                
                                opts = optimset('Display', 'off');
                                d_flow=lsqnonneg(C,D,opts)';
                                d_flowSort=d_flow(clc_nZones+1:end);
                                for d_i=1:clc_nZones
                                    d_flowSort((d_i-1)*clc_nZones+d_i)=-sum(d_flowSort((d_i-1)*clc_nZones+1:d_i*clc_nZones));
                                end
                                clc_flow{5}{d_imp,d_conc}{1,d_seqA}(d_seqV,:,d_noise)=d_flowSort;
                            end
                        end
                    end
                end

                for d_seqA=clc_nSeqAverage(1:end)
                    %%%
                    if (clc_afType=='S' || clc_afType=='F')
                        d_seqVlim=clc_nRunSeq-d_seqA;
                    else
                        d_seqVlim=1;
                    end
                    for d_seqV=1:d_seqVlim
                        %%%
                        d_reqImpA=d_reqImp;
                        d_reqImpA(d_reqImpA==3)=[];
                        for d_imp=d_reqImpA
                            %%%
                            % Calculate averaged flows through combining results
                            for d_solve=d_reqSolve
                                clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA}(d_seqV,:,d_noise)=mean(clc_flow{d_solve}{d_imp,d_conc}{1,1}(d_seqV:2:d_seqV+d_seqA-1,:,d_noise),1);
                            end

                            % Calculate averaged flows through combining inputs
                            clc_crossCorrA=[];
                            for d_relZone=1:clc_nZones
                                for d_concZone=1:clc_nZones
                                    clc_crossCorrA(:,d_relZone,d_concZone)=mean(clc_crossCorr{d_imp}(:,d_relZone,d_concZone,d_seqV:d_seqV+d_seqA-1),4);
                                end
                            end

                            cll_crossCorrSum=[];
                            cll_D=[];
                            cll_E=[];
                            cll_X=[];
                            cll_Xach=[];
                            cll_Xflow=[];
                            cll_dt=clc_dth;

                            % Calculate flowrates
                            cll_crossCorrCalc=clc_crossCorrA(1:(floor(clc_seqLength/clc_nZones)*clc_seqMultiple),:,:);
                            cll_sumStart = clc_seqMultiple+1;
                            cll_sumEnd = size(cll_crossCorrCalc,1)-clc_seqMultiple+1;
                            cll_crossCorrCalc=cll_crossCorrCalc(cll_sumStart:cll_sumEnd,:,:);

                            for d_relZone=1:clc_nZones
                                for d_concZone=1:clc_nZones
                                    cll_crossCorrSum(:,d_relZone,d_concZone,1)=cll_crossCorrCalc(2:end,d_relZone,d_concZone)-cll_crossCorrCalc(1:end-1,d_relZone,d_concZone);
                                    cll_crossCorrSum(:,d_relZone,d_concZone,2)=cll_crossCorrCalc(2:end,d_relZone,d_concZone)+cll_crossCorrCalc(1:end-1,d_relZone,d_concZone);
                                end
                            end

                            for d_zone=1:clc_nZones
                                for d_i=1:clc_nZones
                                    for d_j=1:clc_nZones
                                        cll_D(d_i,d_j,d_zone)=sum(cll_crossCorrSum(:,d_i,d_j,2).*cll_crossCorrSum(:,d_i,d_zone,2));
                                    end
                                    cll_E(d_i,d_zone)=sum(cll_crossCorrSum(:,d_i,d_zone,1).*cll_crossCorrSum(:,d_i,d_zone,2));
                                end
                                cll_X(:,d_zone)=cll_D(:,:,d_zone)\cll_E(:,d_zone);
                                cll_Xach(:,d_zone)=(cll_D(:,:,d_zone)\cll_E(:,d_zone))*2/cll_dt;
                                cll_Xflow(:,d_zone)=(cll_D(:,:,d_zone)\cll_E(:,d_zone))*2/cll_dt*clc_zoneVol(d_zone);
                            end

                            for d_i=1:clc_nZones
                                clc_flow{1}{d_imp,d_conc}{2,d_seqA}(d_seqV,(d_i-1)*clc_nZones+1:d_i*clc_nZones,d_noise)=cll_Xflow(d_i,:);    
                            end
                            
                            % Additional solvers
                            if (ismember(2,d_reqSolve) || ismember(3,d_reqSolve) || ismember(4,d_reqSolve) || ismember(5,d_reqSolve))
                                % Prepare inputs for additional solvers
                                cll_crossCorrCalc=reshape(permute(clc_crossCorrA(1:floor(clc_seqLength/clc_nZones)*clc_seqMultiple,:,:),[1 3 2]),floor(clc_seqLength/clc_nZones)*clc_seqMultiple,clc_nZones^2);
                                cll_sumStart = clc_seqMultiple+1;
                                cll_sumEnd = size(cll_crossCorrCalc,1)-clc_seqMultiple+1;
                                cll_crossCorrCalc=cll_crossCorrCalc(cll_sumStart:cll_sumEnd,:);

                                clc_impAve = (cll_crossCorrCalc(1:end-1,:)+cll_crossCorrCalc(2:end,:))*1000/2;
                                clc_impdt = (cll_crossCorrCalc(2:end,:)-cll_crossCorrCalc(1:end-1,:))*1000/clc_dth;

                                if (ismember(2,d_reqSolve))
                                    % Non-linear least squares - SPLIT EQUATIONS      
                                    d_flowSort=[];
                                    for d_zone=1:clc_nZones
                                        x0=ones(clc_nZones,1)*100;
                                        lb=zeros(clc_nZones,1);
                                        ub=ones(clc_nZones,1)*500;
                                        f=@(x)vecObj_Split(x,d_zone,clc_nZones,clc_zoneVol,clc_impAve,clc_impdt);
                                        opts = optimoptions(@lsqnonlin,'Display', 'off');
                                        d_flow=lsqnonlin(f,x0,lb,ub,opts)';
                                        d_flow(d_zone)=-d_flow(d_zone);
                                        d_flowSort(d_zone,:)=d_flow;
                                    end
                                    clc_flow{2}{d_imp,d_conc}{2,d_seqA}(d_seqV,:,d_noise)=reshape(d_flowSort,clc_nZones^2,1);
                                end

                                if (ismember(4,d_reqSolve))
                                    % Linear least squares - SPLIT EQUATIONS      
                                    d_flowSort=[];
                                    for d_zone=1:clc_nZones
                                        % Prepare matrices
                                        Ccell=cell(clc_nZones,clc_nZones);
                                        Ccell(:,:)={zeros(length(clc_impAve(:,1)),1)};
                                        Dcell=cell(clc_nZones,1);
                                        Dcell(:,:)={zeros(length(clc_impAve(:,1)),1)};

                                        for i=1:clc_nZones
                                            Ccell{i,d_zone}=-clc_impAve(:,(i-1)*clc_nZones+d_zone)/clc_zoneVol(d_zone);
                                            for k=1:clc_nZones
                                                if (k~=d_zone)
                                                    Ccell{i,k}=clc_impAve(:,(i-1)*clc_nZones+k)/clc_zoneVol(d_zone);
                                                end
                                            end
                                            Dcell{i,1}=clc_impdt(:,(i-1)*clc_nZones+d_zone);
                                        end

                                        C=cell2mat(Ccell);
                                        D=cell2mat(Dcell);   
                                        opts = optimset('Display', 'off');
                                        d_flow=lsqnonneg(C,D,opts)';
                                        d_flow(d_zone)=-d_flow(d_zone);
                                        d_flowSort(d_zone,:)=d_flow;
                                    end
                                    clc_flow{4}{d_imp,d_conc}{2,d_seqA}(d_seqV,:,d_noise)=reshape(d_flowSort,clc_nZones^2,1);
                                end

                                if (ismember(3,d_reqSolve))
                                    % Non-linear least squares - SINGLE EQUATION
                                    d_flowSort=[];
                                    x0=ones(clc_nZones*(clc_nZones+1),1)*100;
                                    lb=zeros(clc_nZones*(clc_nZones+1),1);
                                    ub=ones(clc_nZones*(clc_nZones+1),1)*500;
                                    f=@(x)vecObj(x,clc_nZones,clc_zoneVol,clc_impAve,clc_impdt);
                                    opts = optimoptions(@lsqnonlin,'Display', 'off');
                                    d_flow=lsqnonlin(f,x0,lb,ub,opts)';
                                    d_flowSort=d_flow(clc_nZones+1:end);
                                    for d_i=1:clc_nZones
                                        d_flowSort((d_i-1)*clc_nZones+d_i)=-sum(d_flowSort((d_i-1)*clc_nZones+1:d_i*clc_nZones));
                                    end
                                    clc_flow{3}{d_imp,d_conc}{2,d_seqA}(d_seqV,:,d_noise)=d_flowSort;
                                end

                                if (ismember(5,d_reqSolve))
                                    % Linear least squares - SINGLE EQUATION
                                    d_flowSort=[];

                                    % Prepare matrices
                                    Ccell=cell(clc_nZones^2,clc_nZones^2);
                                    Ccell(:,:)={zeros(length(clc_impAve(:,1)),1)};
                                    Dcell=cell(clc_nZones^2,1);
                                    Dcell(:,:)={zeros(length(clc_impAve(:,1)),1)};

                                    for i=1:clc_nZones
                                        for j=1:clc_nZones
                                                Ccell((i-1)*clc_nZones+j,(i-1)*clc_nZones+1:i*clc_nZones)={-clc_impAve(:,(j-1)*clc_nZones+i)/clc_zoneVol(i)};
                                                Dcell{(i-1)*clc_nZones+j,1}=clc_impdt(:,(j-1)*clc_nZones+i);
                                                for k=1:clc_nZones
                                                    if (k~=i)
                                                        Ccell{(i-1)*clc_nZones+j,(k-1)*clc_nZones+i}=clc_impAve(:,(j-1)*clc_nZones+k)/clc_zoneVol(i);
                                                    end
                                                end
                                        end
                                    end

                                    C=cell2mat(Ccell);
                                    D=cell2mat(Dcell);

                                    C=[zeros(size(C,1),clc_nZones) C];
                                    D=[D; zeros(clc_nZones,1)];

                                    Cext=zeros(clc_nZones,clc_nZones*(clc_nZones+1));

                                    for i=1:clc_nZones
                                        Cext(i,i)=-1;
                                        Cext(i,i*clc_nZones+1:(i+1)*clc_nZones)=1;
                                        for j=1:clc_nZones
                                            if (j~=i)
                                                Cext(i,j*clc_nZones+i)=-1;
                                            end
                                        end
                                    end

                                    C=[C; Cext];      

                                    opts = optimset('Display', 'off');
                                    d_flow=lsqnonneg(C,D,opts)';
                                    d_flowSort=d_flow(clc_nZones+1:end);
                                    for d_i=1:clc_nZones
                                        d_flowSort((d_i-1)*clc_nZones+d_i)=-sum(d_flowSort((d_i-1)*clc_nZones+1:d_i*clc_nZones));
                                    end
                                    clc_flow{5}{d_imp,d_conc}{2,d_seqA}(d_seqV,:,d_noise)=d_flowSort;
                                end
                            end
                        end
                    end
                end
            end
        end
        outB_flow{ref_bPerm}=clc_flow;
    end
    
    mat_outP2.out_prbsConcDisc(1,d_batchRun)=outB_prbsConcDisc;
    mat_outP2.out_simFlowTimeFull(1,d_batchRun)=outB_simFlowTimeFull;
    mat_outP2.out_simFlowTime(1,d_batchRun)=outB_simFlowTime;
    mat_outP2.out_simFlow(1,d_batchRun)=outB_simFlow;    
    mat_outP2.out_flow(1,d_batchRun)=outB_flow;
    
    clear outB_* rB_*;
end

%% Save output
cd ..;

toc