% Process Sensor results

%% Read test description
clear;
%#ok<*FNDSB>
tic

d_reqSolve=[1 2 3 4 5];
d_reqImp=[1 2];

mat_testCell=matfile('TestCell.mat');

%% Test Description
sens_results=mat_testCell.T150511;
sens_concZ=[];
sens_ext=sens_results.ext;
sens_concZ{1}=sens_results.z1;
sens_concZ{2}=sens_results.z2;
sens_nZ1=size(sens_concZ{1},2);
sens_nZ2=size(sens_concZ{2},2);

clc_nZones=sens_results.nZones;
clc_seqLength=sens_results.seqLength;
clc_seqPeriod=sens_results.seqPeriod;
clc_seqMultiple=sens_results.seqMultiple;
clc_releaseRate=[1 1];
clc_zoneVol(1:2)=[2.492 5.314];
clc_zoneVolGain=zeros(1,8);
clc_zoneVolGain(1:2)=1./clc_zoneVol;
clc_sensorResp=60;

% Derived variables
clc_cSeqLength=clc_seqLength*clc_seqMultiple;
clc_dSeqLength=floor(clc_seqLength/clc_nZones)*clc_seqMultiple;
clc_nSeq = 24/clc_seqPeriod;
clc_ndt = clc_nSeq * clc_seqLength * clc_seqMultiple;
clc_dt = clc_seqPeriod*60*60/(clc_seqLength*clc_seqMultiple); % Seconds
clc_dth = clc_seqPeriod/(clc_seqLength*clc_seqMultiple); % Hours

clc_nSeqAverage=1;

clc_flow=[];
clc_fit=[];
clc_crossCorr=[];
clc_crossCorrD=[];

%% Generate input sequences used
clc_genPoly=[]; % Initialise loop variables

switch clc_seqLength
    case 7
        clc_genPoly=[3 2 0];
    case 15
        clc_genPoly=[4 3 0];
    case 31
        clc_genPoly=[5 3 0];
    case 63
        clc_genPoly=[6 5 0];
    case 127
        clc_genPoly=[7 6 0];
    case 255
        clc_genPoly=[8 6 5 4 0];
    case 511
        clc_genPoly=[9 5 0];
end

clc_BS = commsrc.pn('GenPoly', clc_genPoly);
set(clc_BS, 'NumBitsOut', clc_seqLength);
clc_PRBS=generate(clc_BS);

% Tracer release day schedules
for d_zone = 1:clc_nZones
    d_count = 1;
    d_temp=[];
    clc_PRBSshift = circshift(clc_PRBS,round((d_zone-1)*clc_seqLength/clc_nZones));
    % Generate 24hr of tracer schedules for calculations, only with entries for changes in state
    for d_j = 0:(clc_seqLength-1)
        for d_k = 0:(clc_seqMultiple-1)
            % Generate full tracer schedule for input file
            sens_input(d_count,d_zone) = clc_PRBSshift(d_j+1);
            d_count=d_count+1;
        end
    end
end

sens_inputFull=[sens_input; sens_input];
   
%% Prepare inputs for cross correlations
% Cross-correlation input variables
clc_crossCorrInput=[];    
clc_a2=[];
clc_seqPosA=[];
clc_seqPosAi=[];

% Static cross correlation variables
for d_i=1:clc_nZones
    clc_crossCorrInput(1:2*clc_cSeqLength,d_i) = sens_inputFull(1:2*clc_cSeqLength,d_i)*clc_releaseRate(d_i)-(clc_releaseRate(d_i)/2);
    clc_a2(d_i)=(clc_releaseRate(d_i)/2)^2;
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
clc_crossCorr=[];
% Calculate cross correlations and flows for all sequence positions required
d_seqA=5;

for d_z1=1:2
    for d_z2=1:2
        sens_conc=[sens_concZ{1}(:,d_z1) sens_concZ{2}(:,d_z2)];
        sens_comb=(d_z1-1)*sens_nZ1+d_z2;

        for d_seqV=1:6
            if (ismember(1,d_reqImp) || ismember(2,d_reqImp))
                for d_relZone=1:clc_nZones
                    for d_concZone=1:clc_nZones
                        clc_crossCorrOutput = (sens_conc(((d_seqV-1)*clc_cSeqLength)+1:(d_seqV+1)*clc_cSeqLength,d_concZone)-mean(sens_ext(((d_seqV-1)*clc_cSeqLength)+1:(d_seqV+1)*clc_cSeqLength)))/1000000;
                        clc_concMean = mean(clc_crossCorrOutput);
                        clc_gain = clc_concMean/(clc_releaseRate(d_relZone)/2);
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


            for d_imp=d_reqImp
                % Calculate flowrates
                cll_crossCorrSum=[];
                cll_D=[];
                cll_E=[];
                cll_X=[];
                cll_Xach=[];
                cll_Xflow=[];
                cll_dt=clc_dth;

                cll_crossCorrCalc=clc_crossCorr{d_imp}(1:(floor(clc_seqLength/clc_nZones)*clc_seqMultiple),:,:,d_seqV);
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
                    clc_flow{sens_comb,1}{d_imp}{1,d_seqA}(d_seqV,(d_i-1)*clc_nZones+1:d_i*clc_nZones)=cll_Xflow(d_i,:);    
                end

                if (ismember(2,d_reqSolve) || ismember(3,d_reqSolve) || ismember(4,d_reqSolve) || ismember(5,d_reqSolve))
                    % Prepare inputs for additional solvers
                    cll_crossCorrCalc=reshape(permute(clc_crossCorr{d_imp}(1:floor(clc_seqLength/clc_nZones)*clc_seqMultiple,:,:,d_seqV),[1 3 2]),floor(clc_seqLength/clc_nZones)*clc_seqMultiple,clc_nZones^2);
                    cll_sumStart = clc_seqMultiple+1+cll_trim;
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
                        clc_flow{sens_comb,2}{d_imp}{1,d_seqA}(d_seqV,:)=reshape(d_flowSort,clc_nZones^2,1);
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
                        clc_flow{sens_comb,4}{d_imp}{1,d_seqA}(d_seqV,:)=reshape(d_flowSort,clc_nZones^2,1);
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
                        clc_flow{sens_comb,3}{d_imp}{1,d_seqA}(d_seqV,:)=d_flowSort;
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
                        clc_flow{sens_comb,5}{d_imp}{1,d_seqA}(d_seqV,:)=d_flowSort;
                    end
                end
            end
        end
    end
end

% %% Direct impulse simulation
% clc_flowTest=clc_flow{1,5}{2}{1,1}(1,:);
% clc_crossCorrTest=[squeeze(clc_crossCorr{2}(:,1,[1 2],1)) squeeze(clc_crossCorr{2}(:,2,[2 1],1))];
% 
% setup_nModelZones=8;
% clc_tZones=2;
% clc_stepSize=1/3600; % 1 second
% 
% clc_Q=[0 clc_flowTest];
% % Fill in blank flowrates into matrix
% for d_i=clc_tZones^2+1:setup_nModelZones^2
%     clc_Q(d_i+1)=0;
% end
% 
%  % Set up initial blank vent gain matrices
% clc_Qgain=cell(1,setup_nModelZones);
% clc_Qgain(1:setup_nModelZones)={zeros(setup_nModelZones,setup_nModelZones^2)};
% 
% for d_i=1:clc_tZones
%     for d_j=1:clc_tZones
%         clc_Qgain{d_j}(d_i,((d_i-1)*clc_tZones)+d_j)=1;
%     end    
% end
% 
% assignin('base','sim_stepSize',clc_stepSize);
% assignin('base','sim_Q',clc_Q);
% assignin('base','sim_Qgain',clc_Qgain);    
% assignin('base','sim_zoneVolGain',clc_zoneVolGain);
% 
% sim_impT=[];
% sim_impX=[];
% sim_impResp=[];
% sim_impTracer=[];        
% clc_impulse=[];
% 
% for d_i=1:clc_tZones
%     d_imp=zeros(1,setup_nModelZones);
%     d_imp(d_i)=clc_zoneVol(d_i)*1000/1000000; % Will start impulse at 1000ppm
%     assignin('base','sim_imp',d_imp);
%     [sim_impT,sim_impX,sim_impResp,sim_impTracer,sim_impFlow]=sim('Impulse_Para',24);
%     clc_impulse(:,((d_i-1)*clc_tZones)+1:(d_i*clc_tZones))=sim_impResp(:,1:clc_tZones);
% end
% 
% clc_shift=[7/15 1 22/15 2 37/15 3 52/15 4 67/15 5 82/15 6; 8/15 1 23/15 2 38/15 3 53/15 4 68/15 5 83/15 6];
% 
% for d_i=1:size(clc_impulse,2)
%     clc_impulseCorr(:,d_i)=clc_impulse(:,d_i);
%     for d_j=1:12
%         clc_impulseCorr(:,d_i)=clc_impulseCorr(:,d_i)+circshift(clc_impulse(:,floor((d_i-1)/2)*2+rem(d_j+(d_i-1),2)+1),clc_shift(rem(d_i,2)+1,d_j)*clc_seqPeriod*3600);
%     end
% end
        
    

%% Save output
d_procTime=toc
