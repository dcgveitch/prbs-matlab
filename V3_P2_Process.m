% Process Simulink results

%% Read test description
clear;
%#ok<*FNDSB>
tic

[d_upperPath, d_folderTS, ~] = fileparts(pwd);

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'));
mat_outBig=matfile(strcat(d_folderTS(1:11), '__outBig.mat'),'Writable',true);
mat_impulse=matfile(strcat(d_folderTS(1:11), '__impulse.mat'),'Writable',true);

d_solve=[1];
d_reqImp=[1 2 3];
d_reqConc=[1 2 3 4];
d_impSave=0;

d_nBatch=ceil(setup_nSim/setup_batchSize);
setup_batchTrim=setup_batchSize;

if (exist('out_impulseSimDisc','var'))
    d_impulseSim=1;
else
    d_impulseSim=0;
end

for d_batch=1:d_nBatch
    d_batchL=(d_batch-1)*setup_batchSize+1;
    d_batchH=min(d_batchL+setup_batchTrim-1,setup_nSim);
    d_batchSize=d_batchH-d_batchL+1;
    
    outB_input(1,:)=mat_outBig.out_input(1,d_batchL:d_batchH);    
    outB_prbsFlow(1,:)=mat_outBig.out_prbsFlow(1,d_batchL:d_batchH);
    outB_prbsConc(1,:)=mat_outBig.out_prbsConc(1,d_batchL:d_batchH);
    try
        rB_impulse=mat_impulse.out_impulse(1,d_batchL:d_batchH);
        rB_impulseLoad=1;
    catch
        rB_impulse=zeros(1,d_batchSize);
        rB_impulseLoad=0;
    end
    
    % Try loading direct simulated impulse responses
    try
        rB_impulseSim=mat_outBig.out_impulseSim(d_batchL:d_batchH);
        rB_impulseSimLoad=1;
    catch
        rB_impulseSim=zeros(1,d_batchSize);
        rB_impulseSimLoad=0;
    end
    
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
    rB_zoneVolGain=r_zoneVolGain(d_batchL:d_batchH,:);
    rB_mixModel=r_mixModel(d_batchL:d_batchH);
    rB_sensorSpecType=r_sensorSpecType(d_batchL:d_batchH);
    rB_sensorSpecRefs=r_sensorSpecRefs(d_batchL:d_batchH);
    rB_sensorResp=r_sensorResp(d_batchL:d_batchH);
    rB_noiseRefsStr=r_noiseRefsStr(d_batchL:d_batchH);
    rB_afType=r_afType(d_batchL:d_batchH,:);
    rB_randSeeds=r_randSeeds(d_batchL:d_batchH,:);
    
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

        clc_inputFull=outB_input{ref_bPerm};
        sim_prbsConc=outB_prbsConc{ref_bPerm};
        
        clc_flow=[];
        clc_impulse=[];
        clc_impulseFull=[];
        clc_crossCorr=[];
        rB_impulseGen=0;
        
        if (ismember(1,d_solve))
            % Form full set of concentration traces from theoretical            
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
                sim_conc{2}(d_i,1:clc_tZones)=sim_prbsSensConc(d_count,1:clc_tZones);
            end
            sim_prbsConc(round(setup_nDaysStab*24/clc_stepSize)+1:end,1:clc_tZones);

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

            % Sensor concentration traces - attempting to reverse engineer 
            if clc_randSeeds(5)==0
                rng('shuffle');
            else
                rng(clc_randSeeds(5));
            end

            for d_i=1:clc_nZones
                for d_noise=1:in_noiseAveNum
                    sim_conc{3}(:,d_i,d_noise)=((sim_conc{1}(:,d_i)*clc_sensorSpec(d_i,d_noise,1))+clc_sensorSpec(d_i,d_noise,2)+(clc_sensorSpec(d_i,d_noise,3)*(1+sim_conc{1}(:,d_i)./setup_sensRNSDrange).*randn(length(sim_conc{1}(:,d_i)),1)));
                    sim_conc{4}(:,d_i,d_noise)=((sim_conc{2}(:,d_i)*clc_sensorSpec(d_i,d_noise,1))+clc_sensorSpec(d_i,d_noise,2)+(clc_sensorSpec(d_i,d_noise,3)*(1+sim_conc{2}(:,d_i)./setup_sensRNSDrange).*randn(length(sim_conc{2}(:,d_i)),1)));
                end
            end  

            d_a=(clc_dt/3600)/(clc_sensorResp);
            if d_a~=inf
                sim_conc{5}=filter([1 d_a-1], d_a, sim_conc{4});
            else
                sim_conc{5}=sim_conc{4};
            end
            
            % Direct impulse traces
            if (rB_impulseSimLoad==1)
                sim_impulseRaw=rB_impulseSim{ref_bPerm}/1000; % Scaling Factor
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
                
                for d_i=1:clc_nZones
                    for d_j=1:clc_nZones
                        d_count=(d_i-1)*clc_nZones+d_j;
                        clc_crossCorr{3,1}(:,d_i,d_j,1,1,1)=sim_impulse{1}(:,d_count);
                        clc_crossCorr{3,2}(:,d_i,d_j,1,1,1)=sim_impulse{2}(:,d_count);
                        clc_impulse{3,1}{1,1}(1:clc_cSeqLength,d_count,1)=clc_crossCorr{3,1}(1:clc_cSeqLength,d_i,d_j,1,1,1);
                        clc_impulse{3,2}{1,1}(1:clc_cSeqLength,d_count,1)=clc_crossCorr{3,2}(1:clc_cSeqLength,d_i,d_j,1,1,1);
                        for d_noise=1:in_noiseAveNum
                            clc_crossCorr{3,3}(:,d_i,d_j,1,1,d_noise)=((sim_impulse{1}(:,d_count)*clc_sensorSpec(d_j,d_noise,1))+clc_sensorSpec(d_j,d_noise,2)+(clc_sensorSpec(d_j,d_noise,3)*(1+sim_impulse{1}(:,d_count)./setup_sensRNSDrange).*randn(length(sim_impulse{1}(:,d_count)),1)));
                            clc_crossCorr{3,4}(:,d_i,d_j,1,1,d_noise)=((sim_impulse{2}(:,d_count)*clc_sensorSpec(d_j,d_noise,1))+clc_sensorSpec(d_j,d_noise,2)+(clc_sensorSpec(d_j,d_noise,3)*(1+sim_impulse{2}(:,d_count)./setup_sensRNSDrange).*randn(length(sim_impulse{2}(:,d_count)),1)));
                            clc_impulse{3,3}{1,1}(1:clc_cSeqLength,d_count,d_noise)=clc_crossCorr{3,3}(1:clc_cSeqLength,d_i,d_j,1,1,d_noise);
                            clc_impulse{3,4}{1,1}(1:clc_cSeqLength,d_count,d_noise)=clc_crossCorr{3,4}(1:clc_cSeqLength,d_i,d_j,1,1,d_noise);
                        end
                    end
                end
                
                if (ismember(5,d_reqConc))
                    d_a=(clc_dt/3600)/(clc_sensorResp);
                    if d_a~=inf
                        for d_i=1:clc_nZones
                            for d_j=1:clc_nZones
                                d_count=(d_i-1)*clc_nZones+d_j;
                                for d_noise=1:in_noiseAveNum
                                    clc_crossCorr{3,5}(:,d_i,d_j,1,1,d_noise)=filter([1 d_a-1], d_a, clc_crossCorr{3,4}(:,d_i,d_j,1,1,d_noise));
                                    clc_impulse{3,5}{1,1}(1:clc_dSeqLength,d_count,d_noise)=clc_crossCorr{3,5}(1:clc_dSeqLength,d_i,d_j,1,1,d_noise);
                                end
                            end
                        end
                    else
                        clc_crossCorr{3,5}=clc_crossCorr{3,4};
                        clc_impulse{3,5}=clc_crossCorr{3,5};
                    end
                end
            end

            %% Airflows
            outB_simFlowTimeFull{ref_bPerm}=[0:clc_dth:(size(outB_prbsFlow{ref_bPerm},1)-1)*clc_dth]';    

            %% Calculate cross correlations
            % Cross-correlation variables
            clc_crossCorrInput=[];    
            clc_a2=[];
            clc_seqPosA=[];
            clc_seqPosAi=[];

            % Static cross correlation variables
            for d_i=1:clc_nZones
                clc_crossCorrInput(1:((max(clc_nSeqAverage)+1)*clc_cSeqLength),d_i) = clc_inputFull(1:((max(clc_nSeqAverage)+1)*clc_cSeqLength),d_i)*clc_releaseRateT(d_i)-(clc_releaseRateT(d_i)/2);
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

            for d_seqA=clc_nSeqAverage
                if (clc_afType=='S' || clc_afType=='F')
                    d_seqVlim=clc_nRunSeq-d_seqA;
                else
                    if (d_seqA==1)
                        d_seqVlim=max(clc_nSeqAverage);
                    else
                        d_seqVlim=1;
                    end
                end
                for d_seqV=1:d_seqVlim
                    % Clear large variables
                    clc_concMean=[];
                    clc_crossCorrOutput=[];
                    clc_gain=[];
                    clc_crossCorrOutputMat=[];

                    for d_conc=d_reqConc % Required types of concentration
                        % Calculate cross correlations
                        d_count=1;     
                        for d_relZone=1:clc_nZones
                            clc_crossCorrInputMat=spdiags(repmat(clc_crossCorrInput(1:(d_seqA*clc_cSeqLength),d_relZone)',clc_cSeqLength,1),0:-1:-(d_seqA*clc_cSeqLength)+1,(d_seqA+1)*clc_cSeqLength,clc_cSeqLength);
                            for d_concZone=1:clc_nZones
                                for d_noise=1:size(sim_conc{d_conc},3)
                                    clc_crossCorrOutput = sim_conc{d_conc}(((d_seqV-1)*clc_cSeqLength)+1:(d_seqV+d_seqA)*clc_cSeqLength,d_concZone,d_noise)/1000000;
                                    clc_concMean = mean(clc_crossCorrOutput);
                                    clc_gain = clc_concMean/(clc_releaseRateT(d_relZone)/2);
                                    clc_crossCorrOutput = clc_crossCorrOutput-clc_concMean;
                                    clc_crossCorrOutputMat=clc_crossCorrOutput';                  
                                    d_crossCorr=(clc_crossCorrOutputMat*clc_crossCorrInputMat)/(d_seqA*clc_cSeqLength);
                                    clc_crossCorr{1,d_conc}(:,d_relZone,d_concZone,d_seqV,d_seqA,d_noise)=((clc_seqLength*d_crossCorr)+(clc_a2(d_relZone)*clc_gain))/((clc_seqLength+1)*clc_dth*clc_seqMultiple*clc_a2(d_relZone));
                                    clc_impulse{1,d_conc}{d_seqV,d_seqA}(1:clc_dSeqLength,d_count,d_noise)=clc_crossCorr{1,d_conc}(1:clc_dSeqLength,d_relZone,d_concZone,d_seqV,d_seqA,d_noise);
                                    clc_impulseFull{1,d_conc}{d_seqV,d_seqA}(:,d_count,d_noise)=clc_crossCorr{1,d_conc}(:,d_relZone,d_concZone,d_seqV,d_seqA,d_noise);
                                end
                                d_count=d_count+1;
                            end
                        end

                        % Average cross correlations
                        d_count=1;
                        for d_relZone=1:clc_nZones
                            for d_concZone=1:clc_nZones
                                d_offset=circshift(clc_nZones:-1:1,[0 d_relZone]);
                                for d_noise=1:size(sim_conc{d_conc},3)
                                    clc_crossCorr{2,d_conc}(:,d_relZone,d_concZone,d_seqV,d_seqA,d_noise)=zeros(clc_dSeqLength,1);
                                    for d_i=1:clc_nZones
                                        clc_crossCorr{2,d_conc}(:,d_relZone,d_concZone,d_seqV,d_seqA,d_noise)=clc_crossCorr{2,d_conc}(:,d_relZone,d_concZone,d_seqV,d_seqA,d_noise)+clc_crossCorr{1,d_conc}(clc_seqPosAi(d_offset(d_i),d_i):clc_seqPosAi(d_offset(d_i),d_i)+clc_dSeqLength-1,d_offset(d_i),d_concZone,d_seqV,d_seqA,d_noise);
                                    end
                                    clc_crossCorr{2,d_conc}(:,d_relZone,d_concZone,d_seqV,d_seqA,d_noise)=clc_crossCorr{2,d_conc}(:,d_relZone,d_concZone,d_seqV,d_seqA,d_noise)./clc_nZones;
                                    clc_impulse{2,d_conc}{d_seqV,d_seqA}(1:clc_dSeqLength,d_count,d_noise)=clc_crossCorr{2,d_conc}(1:clc_dSeqLength,d_relZone,d_concZone,d_seqV,d_seqA,d_noise);
                                end
                                d_count=d_count+1;
                            end
                        end
                    end
                end
            end

            % Calculate flowrates
            clc_simFlowTime=[];
            clc_simFlow=[];
            cll_crossCorrSum=[];
            cll_crossCorrCalc=[];
            cll_dt=[];

            for d_seqA=clc_nSeqAverage
                if (clc_afType=='S' || clc_afType=='F')
                    d_seqVlim=clc_nRunSeq-d_seqA;
                else
                    if (d_seqA==1)
                        d_seqVlim=max(clc_nSeqAverage);
                    else
                        d_seqVlim=1;
                    end
                end
                clc_simFlowTime{d_seqA}=[((d_seqA+1)*clc_seqPeriod)/2:clc_seqPeriod:(((d_seqA+1)*clc_seqPeriod)/2)+(d_seqVlim-1)*clc_seqPeriod]';
                for d_seqV=1:d_seqVlim
                    clc_simFlow{d_seqA}(d_seqV,:)=mean(outB_prbsFlow{ref_bPerm}(((d_seqV-1)*clc_cSeqLength)+1:(d_seqV+d_seqA)*clc_cSeqLength,:));
                    for d_imp=d_reqImp
                        if (d_imp==3 && (d_seqA+d_seqV>2))
                            continue;
                        end
                        for d_conc=d_reqConc
                            for d_noise=1:size(sim_conc{d_conc},3)
                                cll_crossCorrSum=[];
                                cll_D=[];
                                cll_E=[];
                                cll_X=[];
                                cll_Xach=[];
                                cll_Xflow=[];
                                cll_dt=clc_dth;

                                if (d_imp==3)
                                    cll_crossCorrCalc=clc_crossCorr{d_imp,d_conc}(1:clc_cSeqLength,:,:,d_seqV,d_seqA,d_noise);
                                    cll_sumStart = 1;
                                    cll_sumEnd = size(cll_crossCorrCalc,1);
                                else
                                    cll_crossCorrCalc=clc_crossCorr{d_imp,d_conc}(1:(floor(clc_seqLength/clc_nZones)*clc_seqMultiple),:,:,d_seqV,d_seqA,d_noise);
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
                                    clc_flow{1}{d_imp,d_conc}{d_seqA}(d_seqV,(d_i-1)*clc_nZones+1:d_i*clc_nZones,d_noise)=cll_Xflow(d_i,:);    
                                end
                            end
                        end
                    end
                end
            end
            outB_impulse{ref_bPerm}=clc_impulse;
            outB_impulseFull{ref_bPerm}=clc_impulseFull;
            outB_simFlowTime{ref_bPerm}=clc_simFlowTime;
            outB_simFlow{ref_bPerm}=clc_simFlow;
            rB_impulseGen=1;
        end
            
        if (ismember(2,d_solve))
            if (rB_impulseGen==0 && rB_impulseLoad==1)
                clc_impulse=rB_impulse{ref_bPerm};
            end            
            
            %% Non-linear least squares - SPLIT EQUATIONS
            for d_seqA=clc_nSeqAverage
                if (clc_afType=='S' || clc_afType=='F')
                    d_seqVlim=clc_nRunSeq-d_seqA;
                else
                    if (d_seqA==1)
                        d_seqVlim=max(clc_nSeqAverage);
                    else
                        d_seqVlim=1;
                    end
                end
                for d_seqV=1:d_seqVlim
                    for d_imp=d_reqImp
                        if (d_imp==3 && (d_seqA+d_seqV>2))
                            continue;
                        end
                        for d_conc=d_reqConc
                            clc_impBasis=clc_impulse{d_imp,d_conc}{d_seqV,d_seqA};
                            d_nNoise=size(clc_impBasis,3);
                            if (d_imp==3)
                                cll_sumStart = 2;
                                cll_sumEnd = size(clc_impBasis,1);
                            else
                                cll_sumStart = clc_seqMultiple+1;
                                cll_sumEnd = size(clc_impBasis,1)-clc_seqMultiple+1;
                            end
                            for d_noise=1:d_nNoise
                                d_flowSort=[];
                                clc_impBasisTrim=clc_impBasis(:,:,d_noise);
                                clc_impAve = (clc_impBasisTrim(cll_sumStart:cll_sumEnd-1,:)+clc_impBasisTrim(cll_sumStart+1:cll_sumEnd,:))*1000/2;
                                clc_impdt = (clc_impBasisTrim(cll_sumStart+1:cll_sumEnd,:)-clc_impBasisTrim(cll_sumStart:cll_sumEnd-1,:))*1000/clc_dth;
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
                                d_flowSort=reshape(d_flowSort,clc_nZones^2,1);
                                clc_flow{2}{d_imp,d_conc}{d_seqA}(d_seqV,:,d_noise)=d_flowSort;
                            end
                        end
                    end  
                end
            end
        end
        
        if (ismember(3,d_solve))
            if (rB_impulseGen==0 && rB_impulseLoad==1)
                clc_impulse=rB_impulse{ref_bPerm};
            end
            
            %% Non-linear least squares - SINGLE EQUATION
            for d_seqA=clc_nSeqAverage
                if (clc_afType=='S' || clc_afType=='F')
                    d_seqVlim=clc_nRunSeq-d_seqA;
                else
                    if (d_seqA==1)
                        d_seqVlim=max(clc_nSeqAverage);
                    else
                        d_seqVlim=1;
                    end
                end
                for d_seqV=1:d_seqVlim
                    for d_imp=d_reqImp
                        if (d_imp==3 && (d_seqA+d_seqV>2))
                            continue;
                        end
                        for d_conc=d_reqConc
                            clc_impBasis=clc_impulse{d_imp,d_conc}{d_seqV,d_seqA};
                            d_nNoise=size(clc_impBasis,3);
                            cll_sumStart = clc_seqMultiple+1;
                            cll_sumEnd = size(clc_impBasis,1)-clc_seqMultiple+1;
                            for d_noise=1:d_nNoise
                                clc_impBasisTrim=clc_impBasis(:,:,d_noise);
                                clc_impAve = (clc_impBasisTrim(cll_sumStart:cll_sumEnd-1,:)+clc_impBasisTrim(cll_sumStart+1:cll_sumEnd,:))*1000/2;
                                clc_impdt = (clc_impBasisTrim(cll_sumStart+1:cll_sumEnd,:)-clc_impBasisTrim(cll_sumStart:cll_sumEnd-1,:))*1000/clc_dth;
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
                                clc_flow{3}{d_imp,d_conc}{d_seqA}(d_seqV,:,d_noise)=d_flowSort;
                            end
                        end
                    end  
                end
            end
        end
        outB_flow{ref_bPerm}=clc_flow;
    end
    
    if (ismember(1,d_solve))
        mat_outBig.out_simFlowTimeFull(1,d_batchL:d_batchH)=outB_simFlowTimeFull;
        mat_outBig.out_simFlowTime(1,d_batchL:d_batchH)=outB_simFlowTime;
        mat_outBig.out_simFlow(1,d_batchL:d_batchH)=outB_simFlow;
        if (d_impSave==1)
            mat_impulse.out_impulse(1,d_batchL:d_batchH)=outB_impulse;
            mat_impulse.out_impulseFull(1,d_batchL:d_batchH)=outB_impulseFull;
        end
        mat_outBig.out_prbsConcDisc(1,d_batchL:d_batchH)=outB_prbsConcDisc;
    end
    
    if (rB_impulseLoad==1)
        % Already have results so append new entries
        outBT_flow(1,:)=mat_outBig.out_flow(1,d_batchL:d_batchH);    
        for ref_bPerm = 1:d_batchSize
            outBT_flow{ref_bPerm}(find(~cellfun(@isempty,outB_flow{ref_bPerm})))=outB_flow{ref_bPerm}(find(~cellfun(@isempty,outB_flow{ref_bPerm})));
        end
        mat_outBig.out_flow(1,d_batchL:d_batchH)=outBT_flow;
    else
        mat_outBig.out_flow(1,d_batchL:d_batchH)=outB_flow;
    end
    
    clear outB_* rB_*;
end

%% Save output
cd ..;

toc