% Process Simulink results

%% Read test description
clear;
%#ok<*FNDSB>
tic

[d_upperPath, d_folderTS, ~] = fileparts(pwd);

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'));
load(strcat(d_folderTS(1:11), '__results.mat'));
cd ..;

clear out_aFlow*;

d_solveResults=[1];
d_reqImp=[1 2 3];
d_reqConc=[1 2 3 4];

for ref_perm = 1:setup_nSim
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
    clc_flowSim=r_flowSim{ref_perm};
    
    clc_cSeqLength=clc_seqLength*clc_seqMultiple;
    clc_nSeq = 24/clc_seqPeriod;
    clc_ndt = clc_nSeq * clc_seqLength * clc_seqMultiple;
    clc_dt = clc_seqPeriod*60*60/(clc_seqLength*clc_seqMultiple); % Seconds
    clc_dth = clc_seqPeriod/(clc_seqLength*clc_seqMultiple); % Hours
    clc_nRunSeq=(clc_nDays-setup_nDaysStab)*24/clc_seqPeriod;
    
    clc_simFlowFull=out_prbsFlow{ref_perm};
    clc_simFlowTimeFull=out_simFlowTimeFull{ref_perm};
    
    clc_simFlow=[];
    clc_simFlowTime=[];
    clc_flowResults=[];
    clc_flowFullRef=[];
    
    clc_simFlow{1}=out_simFlow{ref_perm};
    clc_simFlow{2}=out_simFlow{ref_perm};
    clc_simFlow{3}=out_simFlow{ref_perm};
    clc_simFlowTime{1}=out_simFlowTime{ref_perm};
    clc_simFlowTime{2}=out_simFlowTime{ref_perm};
    clc_simFlowTime{3}=out_simFlowTime{ref_perm};
    
    clc_flow=out_flow{ref_perm};
    
    
    % Total External flow
    for d_solve=d_solveResults
        for d_seqA=clc_nSeqAverage
            if (clc_afType=='S' || clc_afType=='F')
                d_seqVlimS=clc_nRunSeq-d_seqA;
            else
                if (d_seqA==1)
                    d_seqVlimS=max(clc_nSeqAverage);
                else
                    d_seqVlimS=1;
                end
            end
            d_3DFlow=cell(1);
            for d_imp=d_reqImp
                if (d_imp==3)
                    d_seqVlim=1;
                    if (d_seqA>1)
                        continue;
                    end
                else
                    d_seqVlim=d_seqVlimS;
                end
                for d_conc=d_reqConc
                    d_3DFlow{d_imp,d_conc}(:,1)=clc_simFlowTime{d_imp}{d_seqA}(1:d_seqVlim,:);
                    d_3DFlow{d_imp,d_conc}(:,2)=-sum(clc_simFlow{d_imp}{d_seqA}(1:d_seqVlim,:),2);
                    d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{d_seqA},3);
                    for d_noise=1:d_nNoise
                        d_3DFlow{d_imp,d_conc}(:,2+d_noise)=-sum(clc_flow{d_solve}{d_imp,d_conc}{d_seqA}(1:d_seqVlim,:,d_noise),2);
                    end
                end
            end

            clc_flowResults{d_solve,1}{1}{1,find(clc_nSeqAverage==d_seqA)}=d_3DFlow;
            
            d_3DFlow=cell(1);
            for d_imp=d_reqImp
                if (d_imp==3)
                    d_seqVlim=1;
                    if (d_seqA>1)
                        continue;
                    end
                else
                    d_seqVlim=d_seqVlimS;
                end
                for d_conc=d_reqConc
                    d_3DFlow{d_imp,d_conc}(:,1)=clc_simFlowTime{d_imp}{d_seqA}(1:d_seqVlim,:);
                    d_base=clc_flowResults{d_solve,1}{1}{1,1}{d_imp,d_conc}(:,2);
                    d_baseFull=full(spdiags(repmat(d_base',d_seqA,1),0:-1:-length(d_base)+1,length(d_base)+d_seqA-1,d_seqA));
                    d_baseFull(:,2:2:d_seqA)=[];
                    d_baseAve=mean(d_baseFull,2);
                    d_3DFlow{d_imp,d_conc}(:,2)=d_baseAve(d_seqA:d_seqA+d_seqVlim-1);
                    d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{1},3);
                    for d_noise=1:d_nNoise
                        d_base=clc_flowResults{d_solve,1}{1}{1,1}{d_imp,d_conc}(:,2+d_noise);
                        d_baseFull=full(spdiags(repmat(d_base',d_seqA,1),0:-1:-length(d_base)+1,length(d_base)+d_seqA-1,d_seqA));
                        %d_baseFull(:,2:2:d_seqA)=[];
                        d_baseAve=mean(d_baseFull,2);
                        d_3DFlow{d_imp,d_conc}(:,2+d_noise)=d_baseAve(d_seqA:d_seqA+d_seqVlim-1);
                    end
                end
            end

            clc_flowResults{d_solve,1}{1}{2,find(clc_nSeqAverage==d_seqA)}=d_3DFlow;
        end
    end
    
    
    if (clc_afType=='S' || clc_afType=='F')
        clc_flowFullRef{1}{1}(:,1)=clc_simFlowTimeFull;
        clc_flowFullRef{1}{1}(:,2)=-sum(clc_simFlowFull,2);
    end

    % Zone flow - Exfiltration
    for d_zone=1:clc_nZones
        for d_solve=d_solveResults
             for d_seqA=clc_nSeqAverage
                 if (clc_afType=='S' || clc_afType=='F')
                    d_seqVlimS=clc_nRunSeq-d_seqA;
                else
                    d_seqVlimS=1;
                end
                d_3DFlow=[];
                for d_imp=d_reqImp
                    if (d_imp==3)
                        d_seqVlim=1;
                        if (d_seqA>1)
                            continue;
                        end
                    else
                    d_seqVlim=d_seqVlimS;
                end
                    for d_conc=d_reqConc
                        d_3DFlow{d_imp,d_conc}(:,1)=clc_simFlowTime{d_imp}{d_seqA}(1:d_seqVlim,:);
                        d_3DFlow{d_imp,d_conc}(:,2)=-sum(clc_simFlow{d_imp}{d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones),2);
                        d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{d_seqA},3);
                        for d_noise=1:d_nNoise
                            d_3DFlow{d_imp,d_conc}(:,2+d_noise)=-sum(clc_flow{d_solve}{d_imp,d_conc}{d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones,d_noise),2);
                        end
                    end
                end

                clc_flowResults{d_solve,2}{1,d_zone}{find(clc_nSeqAverage==d_seqA)}=d_3DFlow;
             end
        end
        
        if (clc_afType=='S' || clc_afType=='F')
            clc_flowFullRef{2}{1,d_zone}(:,1)=clc_simFlowTimeFull;
            clc_flowFullRef{2}{1,d_zone}(:,2)=-sum(clc_simFlowFull(:,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones),2);
        end
    end
    
    % Zone flow - Infiltration
    for d_zone=1:clc_nZones
        for d_solve=d_solveResults
            for d_seqA=clc_nSeqAverage
                if (clc_afType=='S' || clc_afType=='F')
                    d_seqVlimS=clc_nRunSeq-d_seqA;
                else
                    d_seqVlimS=1;
                end
                d_3DFlow=[];
                for d_imp=d_reqImp
                    if (d_imp==3)
                        d_seqVlim=1;
                        if (d_seqA>1)
                            continue;
                        end
                    else
                        d_seqVlim=d_seqVlimS;
                    end
                    for d_conc=d_reqConc
                        d_3DFlow{d_imp,d_conc}(:,1)=clc_simFlowTime{d_imp}{d_seqA}(1:d_seqVlim,:);
                        d_3DFlow{d_imp,d_conc}(:,2)=-clc_simFlow{d_imp}{d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+d_zone);
                        for d_k=1:clc_nZones
                            if (d_k~=d_zone) 
                                d_3DFlow{d_imp,d_conc}(:,2)=d_3DFlow{d_imp,d_conc}(:,2)-clc_simFlow{d_imp}{d_seqA}(1:d_seqVlim,(d_k-1)*clc_nZones+d_zone);
                            end
                        end
                        d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{d_seqA},3);
                        for d_noise=1:d_nNoise
                            d_3DFlow{d_imp,d_conc}(:,2+d_noise)=-clc_flow{d_solve}{d_imp,d_conc}{d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+d_zone,d_noise);
                            for d_k=1:clc_nZones
                                if (d_k~=d_zone) 
                                    d_3DFlow{d_imp,d_conc}(:,2+d_noise)=d_3DFlow{d_imp,d_conc}(:,2+d_noise)-clc_flow{d_solve}{d_imp,d_conc}{d_seqA}(1:d_seqVlim,(d_k-1)*clc_nZones+d_zone,d_noise);
                                end
                            end
                        end
                    end
                end

                clc_flowResults{d_solve,2}{2,d_zone}{find(clc_nSeqAverage==d_seqA)}=d_3DFlow;
            end
        end
        
        if (clc_afType=='S' || clc_afType=='F')
            clc_flowFullRef{2}{2,d_zone}(:,1)=clc_simFlowTimeFull;
            clc_flowFullRef{2}{2,d_zone}(:,2)=-clc_simFlowFull(:,(d_zone-1)*clc_nZones+d_zone);
            for d_k=1:clc_nZones
                if (d_k~=d_zone) 
                    clc_flowFullRef{2}{2,d_zone}(:,2)=clc_flowFullRef{2}{1,d_zone}(:,2)-clc_simFlowFull(:,(d_k-1)*clc_nZones+d_zone);
                end
            end
        end
    end
    
    % Individual flows
    for d_zone1=1:clc_nZones+1
        for d_zone2=1:clc_nZones
            for d_solve=d_solveResults
                for d_seqA=clc_nSeqAverage
                    if (clc_afType=='S' || clc_afType=='F')
                        d_seqVlimS=clc_nRunSeq-d_seqA;
                    else
                        d_seqVlimS=1;
                    end
                    d_3DFlow=[];
                    for d_imp=d_reqImp
                        if (d_imp==3)
                            d_seqVlim=1;
                            if (d_seqA>1)
                                continue;
                            end
                        else
                            d_seqVlim=d_seqVlimS;
                        end
                        for d_conc=d_reqConc
                            d_3DFlow{d_imp,d_conc}(:,1)=clc_simFlowTime{d_imp}{d_seqA}(1:d_seqVlim,:);
                            if (d_zone1==d_zone2) % Zone exfiltration
                                d_3DFlow{d_imp,d_conc}(:,2)=clc_flowResults{d_solve,2}{1,d_zone2}{find(clc_nSeqAverage==d_seqA)}{d_imp,d_conc}(:,2);
                            elseif (d_zone1==clc_nZones+1) % Zone infiltration
                                d_3DFlow{d_imp,d_conc}(:,2)=clc_flowResults{d_solve,2}{2,d_zone2}{find(clc_nSeqAverage==d_seqA)}{d_imp,d_conc}(:,2);
                            else
                                d_3DFlow{d_imp,d_conc}(:,2)=clc_simFlow{d_imp}{d_seqA}(1:d_seqVlim,(d_zone1-1)*clc_nZones+d_zone2);
                            end
                            d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{d_seqA},3);
                            for d_noise=1:d_nNoise
                                if (d_zone1==d_zone2) % Zone exfiltration
                                    d_3DFlow{d_imp,d_conc}(:,2+d_noise)=clc_flowResults{d_solve,2}{1,d_zone2}{find(clc_nSeqAverage==d_seqA)}{d_imp,d_conc}(:,2+d_noise);
                                elseif (d_zone1==clc_nZones+1) % Zone infiltration
                                    d_3DFlow{d_imp,d_conc}(:,2+d_noise)=clc_flowResults{d_solve,2}{2,d_zone2}{find(clc_nSeqAverage==d_seqA)}{d_imp,d_conc}(:,2+d_noise);
                                else
                                    d_3DFlow{d_imp,d_conc}(:,2+d_noise)=clc_flow{d_solve}{d_imp,d_conc}{d_seqA}(1:d_seqVlim,(d_zone1-1)*clc_nZones+d_zone2,d_noise);
                                end
                            end
                        end
                    end
                    clc_flowResults{d_solve,3}{d_zone1,d_zone2}{find(clc_nSeqAverage==d_seqA)}=d_3DFlow;
                end
            end
            
            if (clc_afType=='S' || clc_afType=='F')
                clc_flowFullRef{3}{d_zone1,d_zone2}(:,1)=clc_simFlowTimeFull;
                if (d_zone1==d_zone2) % Zone exfiltration
                    clc_flowFullRef{3}{d_zone1,d_zone2}(:,2)=clc_flowFullRef{2}{1,d_zone2}(:,2);
                elseif (d_zone1==clc_nZones+1) % Zone infiltration
                    clc_flowFullRef{3}{d_zone1,d_zone2}(:,2)=clc_flowFullRef{2}{2,d_zone2}(:,2);
                else
                    clc_flowFullRef{3}{d_zone1,d_zone2}(:,2)=clc_simFlowFull(:,(d_zone1-1)*clc_nZones+d_zone2);
                end
            end
        end
    end
    
    if ((clc_afType=='S' || clc_afType=='F') && exist('out_pftConc','var'))
        % Calculate PFT flowrate and potential bias
        clc_pftFlow(1)=-mean(sum(out_pftFlow{ref_perm}(setup_nDaysStab*24/clc_stepSize:end,:),2));
        clc_pftFlow(2)=sum(mean(out_pftConc{ref_perm}(setup_nDaysStab*24/clc_stepSize:end,:),1).*clc_zoneVol(1:clc_nZones))/sum(clc_zoneVol(1:clc_nZones));
        clc_pftFlow(3)=sum(out_pftTracer{ref_perm}(1,:));
        clc_pftFlow(4)=clc_pftFlow(3)/(clc_pftFlow(2)/1000000);
        clc_pftFlow(5)=clc_pftFlow(4)/clc_pftFlow(1);
        out_aPftResults(ref_perm,:)=clc_pftFlow;
    end
    
    
%     % Total Exiting flow
%     for d_solve=d_solveResults
%         for d_seqA=clc_nSeqAverage
%             if (clc_afType=='S' || clc_afType=='F')
%                 d_seqVlim=clc_nRunSeq-d_seqA;
%             else
%                 d_seqVlim=1;
%             end
%             d_3DFlow=cell(1);
%             for d_method=d_reqImp
%                 for d_conc=d_reqConc
%                     d_3DFlow{d_method,d_conc}(:,1)=clc_simFlowTime{d_method}{d_seqA}(1:d_seqVlim,:);
%                     d_3DFlow{d_method,d_conc}(:,2)=-sum(diag(reshape(clc_simFlow{d_method}{d_seqA}(1:d_seqVlim,:),clc_nZones,clc_nZones)));
%                     d_nNoise=size(clc_flow{d_solve}{d_method,d_conc}{d_seqA},3);
%                     for d_noise=1:d_nNoise
%                         d_3DFlow{d_method,d_conc}(:,2+d_noise)=-sum(diag(reshape(clc_flow{d_solve}{d_method,d_conc}{d_seqA}(1:d_seqVlim,:,d_noise),clc_nZones,clc_nZones)));
%                     end
%                 end
%             end
% 
%             clc_flowResults{d_solve,4}{1}{d_seqA}=d_3DFlow;
%         end
%     end
%     
%     % Total Mixing flow
%     for d_solve=d_solveResults
%         for d_seqA=clc_nSeqAverage
%             if (clc_afType=='S' || clc_afType=='F')
%                 d_seqVlim=clc_nRunSeq-d_seqA;
%             else
%                 d_seqVlim=1;
%             end
%             d_3DFlow=cell(1);
%             for d_method=d_reqImp
%                 for d_conc=d_reqConc
%                     d_3DFlow{d_method,d_conc}(:,1)=clc_simFlowTime{d_method}{d_seqA}(1:d_seqVlim,:);
%                     d_3DFlow{d_method,d_conc}(:,2)=sum(clc_simFlow{d_method}{d_seqA}(1:d_seqVlim,:),2)-sum(diag(reshape(clc_simFlow{d_method}{d_seqA}(1:d_seqVlim,:),clc_nZones,clc_nZones)));
%                     d_nNoise=size(clc_flow{d_solve}{d_method,d_conc}{d_seqA},3);
%                     for d_noise=1:d_nNoise
%                         d_3DFlow{d_method,d_conc}(:,2+d_noise)=sum(clc_flow{d_solve}{d_method,d_conc}{d_seqA}(1:d_seqVlim,:,d_noise),2)-sum(diag(reshape(clc_flow{d_solve}{d_method,d_conc}{d_seqA}(1:d_seqVlim,:,d_noise),clc_nZones,clc_nZones)));
%                     end
%                 end
%             end
% 
%             clc_flowResults{d_solve,4}{2}{d_seqA}=d_3DFlow;
%         end
%     end
%     
%     % Zone flow - Exiting
%     for d_zone=1:clc_nZones
%         for d_solve=d_solveResults
%              for d_seqA=clc_nSeqAverage
%                  if (clc_afType=='S' || clc_afType=='F')
%                     d_seqVlim=clc_nRunSeq-d_seqA;
%                 else
%                     d_seqVlim=1;
%                 end
%                 d_3DFlow=[];
%                 for d_method=d_reqImp
%                     for d_conc=d_reqConc
%                         d_3DFlow{d_method,d_conc}(:,1)=clc_simFlowTime{d_method}{d_seqA}(1:d_seqVlim,:);
%                         d_3DFlow{d_method,d_conc}(:,2)=-clc_simFlow{d_method}{d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+d_zone);
%                         d_nNoise=size(clc_flow{d_solve}{d_method,d_conc}{d_seqA},3);
%                         for d_noise=1:d_nNoise
%                             d_3DFlow{d_method,d_conc}(:,2+d_noise)=-clc_flow{d_solve}{d_method,d_conc}{d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+d_zone,d_noise);
%                         end
%                     end
%                 end
% 
%                 clc_flowResults{d_solve,5}{1,d_zone}{d_seqA}=d_3DFlow;
%              end
%         end
%     end
        
    out_aFlowResults{ref_perm}=clc_flowResults;
    if (clc_afType=='S' || clc_afType=='F')
        out_aFlowFullRef{ref_perm}=clc_flowFullRef;
    end

end

%% Save output
cd Results;
save(strcat(d_folderTS(1:11), '__results.mat'), 'out_*', '-v7.3');
cd ..;

toc