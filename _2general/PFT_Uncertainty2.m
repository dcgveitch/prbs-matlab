% Run batch of Simulink runs for PRBS Tracer Method
clear;
%#ok<*FNDSB>
tic

%% Read permuations from Excel file
d_testN=0;
setup_excelInput='InputParametersPFT.xlsx';
setup_timeStamp=datestr(now,'yymmddTHHMM');
setup_nModelZones=8;
setup_ranSeedOffset=44332211;

% Read in sheets
[in_numSetup,in_txtSetup,in_inputSetup] = xlsread(setup_excelInput, 'Setup');
[in_numMain,in_txtMain,in_inputMain] = xlsread(setup_excelInput, 'Permutations');
[in_numSchFlow,in_txtSchFlow,in_inputSchFlow] = xlsread(setup_excelInput, 'Sch_Flow');
[in_numSchSpect,in_txtSchSpect,in_inputSchSpect] = xlsread(setup_excelInput, 'Sch_Spectrum');
[in_numSchdT,in_txtSchdT,in_inputSchdT] = xlsread(setup_excelInput, 'Sch_dT');
[in_numSchNoise,in_txtSchNoise,in_inputSchNoise] = xlsread(setup_excelInput, 'Sch_Noise');
[in_numGeo,in_txtGeo,in_inputGeo] = xlsread(setup_excelInput, 'Geometry');

% General Setup Variables
setup_runType = in_inputSetup{find(strcmp('runType',in_inputSetup(:,1))),2};
setup_runDesc = in_inputSetup{find(strcmp('runDesc',in_inputSetup(:,1))),2}; 
setup_nDaysStab = in_inputSetup{find(strcmp('nDaysStab',in_inputSetup(:,1))),2};
setup_nDaysRun = in_inputSetup{find(strcmp('nDaysRun',in_inputSetup(:,1))),2};
setup_dT = in_inputSetup{find(strcmp('dT',in_inputSetup(:,1))),2};
setup_stepTime = setup_dT/3600;
setup_mFlowExtLogMean = in_inputSetup{find(strcmp('mFlowExtLogMean',in_inputSetup(:,1))),2};
setup_mFlowExtLogSD = in_inputSetup{find(strcmp('mFlowExtLogSD',in_inputSetup(:,1))),2};
setup_mFlowIntLogMean = in_inputSetup{find(strcmp('mFlowIntLogMean',in_inputSetup(:,1))),2};
setup_mFlowIntLogSD = in_inputSetup{find(strcmp('mFlowIntLogSD',in_inputSetup(:,1))),2};
setup_sensRNSDrange = in_inputSetup{find(strcmp('sensRNSDrange',in_inputSetup(:,1))),2};
setup_batchSize=in_inputSetup{find(strcmp('batchSize',in_inputSetup(:,1))),2};
in_noiseAveNum = in_inputSetup{find(strcmp('noiseAveNum',in_inputSetup(:,1))),2};
in_mcSeqLength=[15 31 63 127 255];
in_mcSeqPeriod=[2 3 4 6 8 12 24];
in_mcNZones=[2 3 5 8];
in_mcZoneVol=[60 52 40 30];

% Set up directory to store run files in
setup_runFolder=[setup_timeStamp '_' setup_runDesc];
mkdir('Runs',setup_runFolder);
copyfile('_INPUT',strcat('Runs/', setup_runFolder, '/Input'));
mkdir(strcat('Runs/', setup_runFolder, '/Results'));
addpath(genpath(strcat('Runs/', setup_runFolder)), '-begin');
cd(strcat('Runs/', setup_runFolder));
cd Results;

mat_outP1=matfile(strcat(setup_timeStamp, '__outP1.mat'));

% Calculate number of permutations (either Monte Carlo or Scheduled)
if (strcmp(setup_runType(1),'M'))
    % If Monte Carlo, number of runs included after MC
    setup_nPerm = size(in_numMain,1);
    setup_nMC = str2double(setup_runType(3:end));
else
    % Else each individual line in the permutation sheet is a run
    setup_nPerm = size(in_numMain,1);
    setup_nMC = 1;
end

setup_nSim =setup_nPerm*setup_nMC;

% Extract individual matrices of variables
ind_nZones=cell2mat(in_inputMain(2:end,find(strcmp('nZones',in_inputMain(1,:))))); 
ind_seqLength=cell2mat(in_inputMain(2:end,find(strcmp('seqLength', in_inputMain(1,:)))));
ind_seqPeriod=cell2mat(in_inputMain(2:end,find(strcmp('seqPeriod',in_inputMain(1,:)))));
ind_seqMultipleStr=in_inputMain(2:end,find(strcmp('seqMultiple',in_inputMain(1,:))));
ind_nSeqAverage=in_inputMain(2:end,find(strcmp('nSeqAverage',in_inputMain(1,:))));
ind_releaseRate=in_inputMain(2:end,find(strcmp('TRelRate',in_inputMain(1,:))));

ind_sensorSpecStr=in_inputMain(2:end,find(strcmp('sensorSpec',in_inputMain(1,:))));
ind_sensorResp=cell2mat(in_inputMain(2:end,find(strcmp('sensorResp',in_inputMain(1,:)))));
ind_noiseRefsStr=in_inputMain(2:end,find(strcmp('schNoiseRefs',in_inputMain(1,:))));
ind_afType=cell2mat(in_inputMain(2:end,find(strcmp('afType',in_inputMain(1,:)))));
ind_afRefsStr=in_inputMain(2:end,find(strcmp('schAfRefs',in_inputMain(1,:))));
ind_afVAmpStr=in_inputMain(2:end,find(strcmp('schAfVAmp',in_inputMain(1,:))));
ind_afVFreqStr=cell2mat(in_inputMain(2:end,find(strcmp('schAfVFreq',in_inputMain(1,:)))));
ind_afLRefsStr=in_inputMain(2:end,find(strcmp('schAfLRefs',in_inputMain(1,:))));
ind_randSeedsStr=in_inputMain(2:end,find(strcmp('randSeeds',in_inputMain(1,:))));

% Zone volume and gain extraction
ind_geoRef=in_inputMain(2:end,find(strcmp('geoRef',in_inputMain(1,:))));
ind_mixModelStr=in_inputMain(2:end,find(strcmp('mixModel',in_inputMain(1,:))));

for d_perm=1:setup_nPerm
    if (isnumeric(ind_geoRef{d_perm})) % Check if input geometry is referenced
        ind_zoneVol(d_perm,:)=zeros(1,setup_nModelZones);
        ind_zoneVolGain(d_perm,:)=zeros(1,setup_nModelZones);
        for d_zone=1:setup_nModelZones
            d_geoRefZone = strcat(num2str(ind_geoRef{d_perm}),'_',num2str(d_zone));
            ind_zoneVol(d_perm,d_zone) = in_inputGeo{find(strcmp(d_geoRefZone,in_inputGeo(:,1))),find(strcmp('zoneVol',in_inputGeo(1,:)))};
            ind_zoneVolGain(d_perm,d_zone) = 1/ind_zoneVol(d_perm,d_zone);
        end
    end
end

%% Generate full matrices of variables for parallel simulation
for d_perm=1:setup_nPerm
    for d_MC=1:setup_nMC
        d_sim=(d_perm-1)*setup_nMC+d_MC;
        disp(['Assembling variables.... ' num2str(d_sim) '/' num2str(setup_nSim) ' ' datestr(now)]);

        % Set basic non-MC variables
        r_noiseRefsStr(d_sim,:)=ind_noiseRefsStr(d_perm,:);
        r_afType(d_sim,:)=ind_afType(d_perm,:);    
        r_sensorResp(d_sim)=ind_sensorResp(d_perm)/3600;
        
        % Set random number seeding
        d_formatStr='%f_';
        d_temp=sscanf(char(ind_randSeedsStr{d_perm}),d_formatStr);
        for d_i=1:5 % 1=Zone Volume 2=Offset Release Rates 3=Flowrates 4=Sensor Specs 5=Sensor Concentrations
            if d_temp(d_i)==1
                r_randSeeds(d_sim,d_i)=setup_ranSeedOffset+d_i; % Fixed seed
            elseif d_temp(d_i)==2
                r_randSeeds(d_sim,d_i)=setup_ranSeedOffset+d_perm*100000+d_i; % Stepped Perm / Fixed MC
            elseif d_temp(d_i)==3
                r_randSeeds(d_sim,d_i)=setup_ranSeedOffset+d_MC*100+d_i; % Fixed Perm / Stepped MC
            elseif d_temp(d_i)==4
                r_randSeeds(d_sim,d_i)=setup_ranSeedOffset+d_perm*100000+d_MC*100+d_i; % Stepped Perm / Stepped MC
            else
                r_randSeeds(d_sim,d_i)=0; % Random seed
            end
        end
        
        % Sequence Length MC
        if (strcmp(ind_seqLength(d_perm,1),'M'))
            r_seqLength(d_sim)=in_mcSeqLength(ceil(rand(1)*length(in_mcSeqLength)));
        else
            r_seqLength(d_sim)=ind_seqLength(d_perm);
        end

        % Sequence Period MC
        if (strcmp(ind_seqPeriod(d_perm,1),'M'))
            r_seqPeriod(d_sim)=in_mcSeqPeriod(ceil(rand(1)*length(in_mcSeqPeriod)));
        else
            r_seqPeriod(d_sim)=ind_seqPeriod(d_perm);
        end

        % Sequence Multiple
        if (strcmp(ind_seqMultipleStr{d_perm}(1),'F'))
            r_seqMultiple(d_sim)=max(1,floor(r_seqPeriod(d_sim)/r_seqLength(d_sim)/(str2double(ind_seqMultipleStr{d_perm}(3:end))/3600)));
            r_seqMultipleType(d_sim)=1;
        else
            r_seqMultiple(d_sim)=ind_seqMultipleStr{d_perm};
            r_seqMultipleType(d_sim)=0;
        end
        
        % Simulation Step Size
        r_stepSize(d_sim)=r_seqPeriod(d_sim)/(r_seqLength(d_sim)*r_seqMultiple(d_sim))/ceil(r_seqPeriod(d_sim)/(r_seqLength(d_sim)*r_seqMultiple(d_sim))/setup_stepTime); % Round divisor to make timestep a factor of the PRBS dT
        
        % Sequence Averaging
        if (strcmp(ind_nSeqAverage{d_perm}(1),'M'))
            d_temp = sscanf(char(ind_nSeqAverage{d_perm}(3:end)),'%f'); % Highest hours averaged
            d_temp = (d_temp/r_seqPeriod(d_sim)); % Go from hours to number of sequences
            for d_i=0:9
                d_list(d_i+1)=ceil((d_temp)/2^d_i);
            end
            d_list=d_list-1;
            if (max(d_list)>3)
                d_list=[d_list 3];
            end
            d_list(d_list==2)=[];
            d_list(d_list==0)=[];
            r_nSeqAverage{d_sim}=unique(d_list); 
        else
            r_nSeqAverage{d_sim}=sscanf(char(ind_nSeqAverage{d_perm}),'%f')'; % Extracts number from list
        end
            
        % Number of days to simulate
        r_nDays(d_sim)=setup_nDaysStab+max(setup_nDaysRun,ceil(r_seqPeriod(d_sim)*(max(r_nSeqAverage{d_sim})+1)/24));

        % Number of Zones
        if (strcmp(ind_nZones(d_perm,1),'M'))
            if (r_seqLength(d_sim)<31)
                d_zoneMC=ceil(rand*(length(in_mcNZones)-1));
            else
                d_zoneMC=ceil(rand*length(in_mcNZones));
            end
            r_nZones(d_sim)=in_mcNZones(d_zoneMC);
        else
            r_nZones(d_sim)=ind_nZones(d_perm);
            d_zoneMC=find(in_mcNZones==ind_nZones(d_perm));
        end

        % Zone Volume & Gain
        if r_randSeeds(d_sim,1)==0
            rng('shuffle');
        else
            rng(r_randSeeds(d_sim,1));
        end
        
        if (strcmp(ind_geoRef(d_perm,1),'M'))
            for d_zone=1:setup_nModelZones
                if (d_zone<=r_nZones(d_sim))
                    r_zoneVol(d_sim,d_zone)=0.5*in_mcZoneVol(d_zoneMC)+rand*in_mcZoneVol(d_zoneMC);
                    r_zoneVolGain(d_sim,d_zone)=1/r_zoneVol(d_sim,d_zone);
                else
                    r_zoneVol(d_sim,d_zone)=0;
                    r_zoneVolGain(d_sim,d_zone)=0;
                end
            end
        else
            r_zoneVol(d_sim,:)=ind_zoneVol(d_perm,:);
            r_zoneVolGain(d_sim,:)=ind_zoneVolGain(d_perm,:);
        end

        % Tracer release rates
        d_formatStr = '%f_';
        d_releaseT=zeros(1,setup_nModelZones);
        if r_randSeeds(d_sim,2)==0
            rng('shuffle');
        else
            rng(r_randSeeds(d_sim,2));
        end
        
        if (strcmp(ind_releaseRate{d_perm}(1),'M'))
            d_rateSD=sscanf(char(ind_releaseRate{d_perm}(3:end)),'%f %f');
            for d_zone=1:r_nZones(d_sim)
                d_releaseT(d_zone)=in_mcZoneVol(d_zoneMC)*0.5*1000/1000000*2; % Approximately 1000ppm at 0.5ACH
            end
        else
            d_rateSD=sscanf(char(ind_releaseRate{d_perm}(1:9)),'%f %f');
            d_rate=sscanf(char(ind_releaseRate{d_perm}(10:end)),'%f_');
            if (length(d_rate)<r_nZones(d_sim)) % Not all zones supplied, use 1st only
                for d_zone=1:r_nZones(d_sim)
                    d_releaseT(d_zone)=d_rate(1);
                end
            else
                d_releaseT(1:length(d_rate))=d_rate;
            end
        end
        d_releaseBias=d_releaseT*d_rateSD(1)*randn; % Add bias (consistent across zones)
        d_release=d_releaseT+d_releaseBias+d_releaseT*d_rateSD(2).*randn(1,length(d_releaseT)); % Add random uncertainty
        
        r_releaseRate{d_sim} = d_release; % Actual release rate
        r_releaseRateT{d_sim}= d_releaseT; % Theoretical release rate for calculations

        % Sensor Spec MC
        if (strcmp(ind_sensorSpecStr{d_perm}(1),'M'))
            d_temp=sscanf(char(ind_sensorSpecStr{d_perm}(3:end)),'%f');
            r_sensorSpecRefs{d_sim}=d_temp;
            r_sensorSpecType(d_sim)=1;
        else
            d_temp=sscanf(char(ind_sensorSpecStr{d_perm}),'%f,%f,%f_');
            r_sensorSpecRefs{d_sim}=d_temp;
            r_sensorSpecType(d_sim)=2;
        end

        %Mixing Model
        if (strcmp(ind_mixModelStr{d_perm}(1),'Y'))
            r_mixModelRefs(d_sim,:) = str2double(strsplit(ind_mixModelStr{d_perm}(3:end)));
            r_mixModel(d_sim)=1;
            r_tZones(d_sim)=2*r_nZones(d_sim);
            r_mixModelVol(d_sim)=r_mixModelRefs(d_sim,1);
            r_mixModelRate(d_sim)=r_mixModelRefs(d_sim,2);
            for d_zone=1:r_nZones(d_sim)
                r_zoneVol(d_sim,d_zone+r_nZones(d_sim))=r_zoneVol(d_sim,d_zone)*r_mixModelVol(d_sim);
                r_zoneVolGain(d_sim,d_zone+r_nZones(d_sim))=1/r_zoneVol(d_sim,d_zone+r_nZones(d_sim));
                r_zoneVolGain(d_sim,d_zone)=1/(r_zoneVol(d_sim,d_zone)-r_zoneVol(d_sim,d_zone+r_nZones(d_sim)));
            end
        else
            r_mixModel(d_sim)=0;
            r_tZones(d_sim)=r_nZones(d_sim);
            r_mixModelVol(d_sim)=0;
            r_mixModelRate(d_sim)=0;
        end 

        % Airflow MC

        if r_randSeeds(d_sim,3)==0
            rng('shuffle');
        else
            rng(r_randSeeds(d_sim,3));
        end
        
        if (r_afType(d_sim,1)=='M')
            % If airflows are a MC variable
            % Extract permutation airflow references
            d_formatStr='%f';
            d_temp=sscanf(char(ind_afRefsStr{d_perm}),d_formatStr);
            r_afRefs{d_sim}=d_temp;             
            
            d_failFlag=1;
            d_failReselect=-1;
            while d_failFlag>0
                d_ext=random('logn',log(r_afRefs{d_sim}(1)),setup_mFlowExtLogSD,1,1);
                d_int=random('logn',log(r_afRefs{d_sim}(2)),setup_mFlowIntLogSD,1,1);
                
                d_totalVol=sum(r_zoneVol(d_sim,:))*d_ext;
                d_split=sort(rand(r_nZones(d_sim)-1,1));
                d_split=[0; d_split; 1];
                d_split=(d_split(2:end)-d_split(1:end-1))';
%                 d_split=(d_split(2:end)-d_split(1:end-1))'.*r_zoneVol(d_sim,1:r_nZones(d_sim)); Remove volume weighting
                d_split=d_split*d_totalVol/sum(d_split);
                r_flowExt(d_sim,1:r_nZones(d_sim))=d_split./r_zoneVol(d_sim,1:r_nZones(d_sim));

                d_failReselect=d_failReselect+1;
                d_failCount=-1;
                while d_failFlag>0 && d_failCount<50000
                    d_failFlag=0;
                    d_flowTest=[];
                    d_failCount=d_failCount+1;
                    
                    % New weighted assignment of internal flows
                    d_totalVol=sum(r_zoneVol(d_sim,:))*d_int;
                    d_split=sort(rand(r_nZones(d_sim)*(r_nZones(d_sim)-1)-1,1));
                    d_split=[0; d_split; 1];
                    d_zoneWeight=[];
                    for d_zone1=1:r_nZones(d_sim)
                        for d_zone2=1:r_nZones(d_sim)
                            if (d_zone1~=d_zone2)
                                d_zoneWeight(end+1)=r_zoneVol(d_sim,d_zone1)*r_zoneVol(d_sim,d_zone2);
                            end
                        end
                    end
                    d_split=(d_split(2:end)-d_split(1:end-1))';
%                     d_split=(d_split(2:end)-d_split(1:end-1))'.*d_zoneWeight; Remove volume weighting
                    d_split=d_split*d_totalVol/sum(d_split);
                    
                    d_count=1;                              
                    for d_zone1=1:r_nZones(d_sim)
                        for d_zone2=1:r_nZones(d_sim)
                            if (d_zone1==d_zone2)
                                d_flowTest(d_zone1,d_zone2) = r_flowExt(d_sim,d_zone1)*r_zoneVol(d_sim,d_zone1);
                            else
                                d_flowTest(d_zone1,d_zone2) = d_split(d_count);
                                r_flowIntSplit{d_sim}(d_zone1,d_zone2)=d_split(d_count)/r_zoneVol(d_sim,d_zone1);
                                d_count=d_count+1;
                            end
                        end
                    end

                    for d_zone=1:r_nZones(d_sim) % Check internal flows are not driving extra exfiltration
                        if (sum(d_flowTest(d_zone,:),2)<(sum(d_flowTest(:,d_zone),1)-d_flowTest(d_zone,d_zone)))
                            d_failFlag=1;
                        end
                    end
                end
            end
            r_failReselect(d_sim)=d_failReselect;
            d_flow=d_flowTest;            
            
        elseif (r_afType(d_sim,1)=='S')
            % Extract permutation airflow references
            d_formatStr='';
            for d_zone=1:r_nZones(d_sim)
                d_formatStr = strcat(d_formatStr,'%f');
                if d_zone<r_nZones(d_sim)
                    d_formatStr = strcat(d_formatStr,',');
                end
            end
            d_formatStr = strcat(d_formatStr,'_');

            d_temp=sscanf(char(ind_afRefsStr{d_perm}),d_formatStr);
            r_afRefs{d_sim}=d_temp;            
            
            % If airflows are scheduled
            d_formatStr = '%f_%f';
            
            for d_zone1=1:r_nZones(d_sim)
                for d_zone2=1:r_nZones(d_sim)
                    d_i = (d_zone1-1)*r_nZones(d_sim) + d_zone2;
                    d_schLen = size(in_txtSchFlow(:,r_afRefs{d_sim}(d_i)),1);
                    d_temp=[];

                    % Scan values of flow schedule
                    for d_lineMark=1:(d_schLen)
                        d_schStr=in_txtSchFlow{d_lineMark,r_afRefs{d_sim}(d_i)};
                        if (isnan(d_schStr)==0)
                            d_temp(d_lineMark,:)=sscanf(d_schStr,d_formatStr);                    
                        end
                    end
                    
                    d_temp(:,2)=d_temp(:,2)*r_zoneVol(d_sim,d_zone1); % Convert ACH to m3/hr

                    d_x=d_temp(:,1);
                    d_y=d_temp(:,2);
                    d_xx=0:30/3600:r_nDays(d_sim)*24; % Interpolate via spline to 30 second steps
                    d_xx=d_xx';
                    d_yy=pchip(d_x,d_y,d_xx);

                    d_flow(d_zone1,d_zone2,:)=d_yy;
                end
            end
            clear d_x* d_y*;
        elseif (r_afType(d_sim,1)=='F')
            % If airflows are taken from freq spectrum
            
            % Extract permutation airflow references
            d_formatStr='%f';
            d_temp=sscanf(char(ind_afRefsStr{d_perm}),d_formatStr);
            r_afRefs{d_sim}=d_temp; 

            r_freqSpectrum{d_sim}=in_numSchSpect;
            d_freq=linspace(-3,log10(4),100);
            d_freq(end+1)=d_freq(end)+(d_freq(end)-d_freq(end-1));
            r_windBasis{d_sim}(:,1)=10.^(d_freq);
            r_windBasis{d_sim}(:,2)=r_windBasis{d_sim}(:,1)*2*pi;
            r_windBasis{d_sim}(:,3)=pchip(r_freqSpectrum{d_sim}(:,1),r_freqSpectrum{d_sim}(:,2),r_windBasis{d_sim}(:,1));
            r_windBasis{d_sim}(:,4)=r_windBasis{d_sim}(:,3)./r_windBasis{d_sim}(:,2);
            r_windBasis{d_sim}(1:end-1,6)=2/pi*(0.5*(r_windBasis{d_sim}(1:end-1,4)+r_windBasis{d_sim}(2:end,4)).*(r_windBasis{d_sim}(2:end,2)-r_windBasis{d_sim}(1:end-1,2))).^0.5;
            r_windBasis{d_sim}=r_windBasis{d_sim}(1:end-1,:);

            % Generate 30 second dT from 5 minute empirical data            
            d_schLen=r_nDays(d_sim)*24*12+1; % Length of 5 minute empirical dT data required            
            d_x=0:5/60:r_nDays(d_sim)*24;
            d_offset=floor(rand(1)*(length(in_numSchdT)-d_schLen));
            d_y=in_numSchdT(d_offset:d_offset+d_schLen-1,1);
            d_xx=0:30/3600:r_nDays(d_sim)*24; % Interpolate via spline to 30 second steps
            d_xx=d_xx';
            r_stackdT{d_sim}=pchip(d_x,d_y,d_xx);     

            clear d_x* d_y*;

            %% NEW SELECTION METHOD
            % Generate average ACH Values as per MC approach
            d_failFlag=1;
            d_failReselect=-1;
            while d_failFlag>0
                d_ext=r_afRefs{d_sim}(1);
                d_int=r_afRefs{d_sim}(2);
                
                % Split external flow                
                d_totalVol=sum(r_zoneVol(d_sim,:))*d_ext;
                d_split=sort(rand(r_nZones(d_sim)-1,1));
                d_split=[0; d_split; 1];
                d_split=(d_split(2:end)-d_split(1:end-1))'.*r_zoneVol(d_sim,1:r_nZones(d_sim));
                d_split=d_split*d_totalVol/sum(d_split);
                r_flowExt(d_sim,1:r_nZones(d_sim))=d_split./r_zoneVol(d_sim,1:r_nZones(d_sim));

                d_failReselect=d_failReselect+1;
                d_failCount=-1;
                while d_failFlag>0 && d_failCount<50000
                    d_failFlag=0;
                    d_flowTest=[];
                    d_failCount=d_failCount+1;
                    
                    % New weighted assignment of internal flows
                    d_totalVol=sum(r_zoneVol(d_sim,:))*d_int;
                    d_split=sort(rand(r_nZones(d_sim)*(r_nZones(d_sim)-1)-1,1));
                    d_split=[0; d_split; 1];
                    d_zoneWeight=[];
                    for d_zone1=1:r_nZones(d_sim)
                        for d_zone2=1:r_nZones(d_sim)
                            if (d_zone1~=d_zone2)
                                d_zoneWeight(end+1)=r_zoneVol(d_sim,d_zone1)*r_zoneVol(d_sim,d_zone2);
                            end
                        end
                    end
                    d_split=(d_split(2:end)-d_split(1:end-1))'.*d_zoneWeight;
                    d_split=d_split*d_totalVol/sum(d_split);
                    
                    d_count=1;                              
                    for d_zone1=1:r_nZones(d_sim)
                        for d_zone2=1:r_nZones(d_sim)
                            if (d_zone1==d_zone2)
                                d_flowTest(d_zone1,d_zone2) = r_flowExt(d_sim,d_zone1)*r_zoneVol(d_sim,d_zone1);
                            else
                                d_flowTest(d_zone1,d_zone2) = d_split(d_count);
                                r_flowIntSplit{d_sim}(d_zone1,d_zone2)=d_split(d_count)/r_zoneVol(d_sim,d_zone1);
                                d_count=d_count+1;
                            end
                        end
                    end

                    for d_zone=1:r_nZones(d_sim) % Check internal flows are not driving extra exfiltration
                        if (sum(d_flowTest(d_zone,:),2)<(sum(d_flowTest(:,d_zone),1)-d_flowTest(d_zone,d_zone)))
                            d_failFlag=1;
                        end
                    end
                    
                    if (d_failFlag==1)
                        continue;
                    end

                    d_flow=[];
                    d_flowTest=[];
                    for d_zone1=1:r_nZones(d_sim)
                        for d_zone2=1:r_nZones(d_sim)
                            d_i = (d_zone1-1)*r_nZones(d_sim) + d_zone2;
                            d_check=0;
                            while (d_check==0)
                                r_windBasis{d_sim}(:,5)=rand(size(r_windBasis{d_sim},1),1)*2*pi-pi;
                                d_wind=[];
                                parfor d_count=1:round((r_nDays(d_sim)*24*3600/30+1)*2)
                                    d_t=(d_count-1)*30/3600;
                                    d_wind(d_count)=sum(r_windBasis{d_sim}(:,6).*cos(r_windBasis{d_sim}(:,2).*d_t+r_windBasis{d_sim}(:,5)));
                                end
                                d_wind=d_wind*(3.5/max(d_wind));
                                d_wind=d_wind+0.9791; % Offset to mean windspeed
                                d_wind(d_wind<0) = []; % Chop out zeros
                                if (length(d_wind)>r_nDays(d_sim)*24*3600/30+1) % Verify that truncated wind speed is long enough
                                    d_check=1;
                                end
                            end
                            d_wind=d_wind(1:r_nDays(d_sim)*24*3600/30+1)'; % Trim

                            if (d_zone1==d_zone2)
                                d_windAvACH = r_flowExt(d_sim,d_zone1);
                            else
                                d_windAvACH = r_flowIntSplit{d_sim}(d_zone1,d_zone2);
                            end

                            % Use combination of pressures in quadrature to get flowrate
                            d_A=0.9;
                            d_B=1;
                            d_n=0.7;
                            d_B1=-0.333;
                            d_Qs=d_A*r_stackdT{d_sim}.^d_n;
                            d_Qw=d_B*d_wind.^(2*d_n);
                            d_flowVar=(d_Qs.^(1/d_n)+d_Qw.^(1/d_n)+(d_B1*(d_Qs.*d_Qw).^(1/(2*d_n)))).^d_n;
                            d_flowVar=d_flowVar*d_windAvACH/geomean(d_flowVar); % Set average
                            d_flowVar=d_flowVar*r_zoneVol(d_sim,d_zone1); % Convert ACH to m3/hr
                            d_flow(d_zone1,d_zone2,:)=d_flowVar;
                        end
                    end
                    
                    % Correct diagonal terms to total flowrate exiting zone by summing rows
                    d_flowTest=d_flow;
                    for d_zone1=1:r_tZones(d_sim)
                        d_flowTest(d_zone1,d_zone1,:)=-(max(sum(d_flow(d_zone1,:,:),2),(sum(d_flow(:,d_zone1,:),1)-d_flow(d_zone1,d_zone1,:))));
                        d_flowExcess(d_zone1,:)=(sum(d_flow(:,d_zone1,:),1)-d_flow(d_zone1,d_zone1,:))./sum(d_flow(d_zone1,:,:),2)-1;
                        d_flowRaw(:,(d_zone1-1)*(r_tZones(d_sim)+1)+1)=squeeze(-sum(d_flowTest(d_zone1,:,:),2));
                        for d_zone2=1:r_tZones(d_sim)
                            d_flowRaw(:,(d_zone1-1)*(r_tZones(d_sim)+1)+d_zone2+1)=squeeze(d_flow(d_zone1,d_zone2,:));
                        end
                    end
                    
                    d_flowExtTotal=geomean(-squeeze(sum(sum(d_flowTest(:,:,:),2),1)));
                    d_error=(d_flowExtTotal-(sum(r_zoneVol(d_sim,:))*d_ext))/(sum(r_zoneVol(d_sim,:))*d_ext);

                    if(d_error>0.1)
                        disp(['Fail! Count: ' num2str(d_failCount) '/2000....Error ' num2str(d_error*100) '%']);
                        d_failFlag=1;
                    end
                end
            end            
            
        else
            % Extract permutation airflow references
            d_formatStr='';
            for d_zone=1:r_nZones(d_sim)
                d_formatStr = strcat(d_formatStr,'%f');
                if d_zone<r_nZones(d_sim)
                    d_formatStr = strcat(d_formatStr,',');
                end
            end
            d_formatStr = strcat(d_formatStr,'_');

            d_temp=sscanf(char(ind_afRefsStr{d_perm}),d_formatStr);
            r_afRefs{d_sim}=d_temp; 
            
            % Read in ACH values taken directly from afRefs)
            for d_zone1=1:r_nZones(d_sim)
                for d_zone2=1:r_nZones(d_sim)
                    d_i = (d_zone1-1)*r_nZones(d_sim) + d_zone2;
                    d_flow(d_zone1,d_zone2) = r_afRefs{d_sim}(d_i)*r_zoneVol(d_sim,d_zone1);
                end
            end
        end
        
        % Add mixing flows if enabled
        if (r_mixModel(d_sim)==1)
            for d_zone1=1:r_tZones(d_sim)
                for d_zone2=1:r_tZones(d_sim)
                    if ((d_zone1>r_nZones(d_sim)) | (d_zone2>r_nZones(d_sim)))
                        if (d_zone1-d_zone2==r_nZones(d_sim))
                            d_flow(d_zone1,d_zone2,:) = r_mixModelRate(d_sim)*r_zoneVol(d_sim,d_zone1);
                        elseif (d_zone2-d_zone1==r_nZones(d_sim))
                            d_flow(d_zone1,d_zone2,:) = r_mixModelRate(d_sim)*r_zoneVol(d_sim,d_zone2);
                        else
                            d_flow(d_zone1,d_zone2,:) = 0;
                        end
                    end
                end
            end
        end
        
        % Correct diagonal terms to total flowrate exiting zone by summing rows
        for d_zone=1:r_tZones(d_sim)
            d_flow(d_zone,d_zone,:)=-(max(sum(d_flow(d_zone,:,:),2),(sum(d_flow(:,d_zone,:),1)-d_flow(d_zone,d_zone,:))));
        end

        % Linear vector of flowrates for simulation input
        for d_zone1=1:r_tZones(d_sim)
            for d_zone2=1:r_tZones(d_sim)
                d_i = (d_zone1-1)*r_tZones(d_sim) + d_zone2;
                r_flowSim{d_sim}(:,d_i) = squeeze(d_flow(d_zone1,d_zone2,:));
            end
        end
    end
end




%% Save Simulink results
save(strcat(setup_timeStamp, '_setup.mat'), 'setup_*', 'in_*', 'r_*','-v7.3');
cd ..;
cd ../..;

d_procTime=toc
mat_outP1.d_procTime=d_procTime;