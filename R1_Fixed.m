% Process Simulink results

%% Read test description
clear;
%#ok<*FNDSB>
tic

[d_upperPath, d_folder, ~] = fileparts(pwd);
if d_folder(2)=='_', d_folderTS=d_folder(5:15);
else d_folderTS=d_folder(1:11); end

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'))
mat_outP3=matfile(strcat(d_folderTS, '__outP3.mat'),'Writable',true);

for d_i=1:length(r_flowSim)
    clc_nZones=r_nZones(d_i);
    d_flow=reshape(r_flowSim{d_i}(1:clc_nZones^2),clc_nZones,clc_nZones)'.*-1;    

    %% 1. PFT Single - No Noise
    out_refTracer{d_i,1}(1:clc_nZones,1)=r_releaseRateT{d_i}(1:clc_nZones)/2;
    out_refConc{d_i,1}(1:clc_nZones,1)=inv(d_flow)*out_refTracer{d_i,1}(1:clc_nZones,1)*1000000;
    d_conc=sum(out_refConc{d_i,1}(1:clc_nZones).*r_zoneVol(d_i,1:clc_nZones)')/sum(r_zoneVol(d_i,1:clc_nZones)); 
    out_refFlow{d_i,1}(1)=sum(out_refTracer{d_i,1}(1:clc_nZones,1))/(d_conc/1000000);
    
    %% 2. PFT Multi - No Noise
    out_refTracer{d_i,2}(1:clc_nZones^2,1)=0;
    for d_zone=1:clc_nZones 
        out_refTracer{d_i,2}((d_zone-1)*clc_nZones+d_zone,1)=r_releaseRateT{d_i}(d_zone)/2;
        out_refConc{d_i,2}((d_zone-1)*clc_nZones+1:d_zone*clc_nZones,1)=d_flow\out_refTracer{d_i,2}((d_zone-1)*clc_nZones+1:d_zone*clc_nZones,1)*1000000;
    end
    d_Q=[];
    for d_zone=1:clc_nZones
        d_Q(d_zone,1:clc_nZones)=(out_refConc{d_i,2}((d_zone-1)*clc_nZones+1:d_zone*clc_nZones,1)/1000000)/out_refTracer{d_i,2}((d_zone-1)*clc_nZones+d_zone,1);
    end
    out_refFlow{d_i,2}(1:clc_nZones^2,1)=-reshape(inv(d_Q),clc_nZones^2,1);
    
    %% 3. Const Conc - No Noise
    out_refConc{d_i,3}(1:clc_nZones,1)=ones(clc_nZones,1)*(1000/1000000);
    out_refTracer{d_i,3}(1:clc_nZones,1)=d_flow*out_refConc{d_i,3}(1:clc_nZones,1);
    out_refFlow{d_i,3}(1:clc_nZones,1)=out_refTracer{d_i,3}(1:clc_nZones,1)./out_refConc{d_i,3}(1:clc_nZones,1);

    for d_noise=1:in_noiseAveNum
        %% 4. PFT Single - Noise
        out_refTracer{d_i,4}(1:clc_nZones,d_noise)=out_refTracer{d_i,1}(1:clc_nZones,1).*(randn(clc_nZones,1)*0.08+1);
        out_refConc{d_i,4}(1:clc_nZones,d_noise)=(d_flow\out_refTracer{d_i,4}(1:clc_nZones,d_noise)*1000000).*(randn(clc_nZones,1)*0.08+1);
        d_conc=sum(out_refConc{d_i,4}(1:clc_nZones,d_noise).*r_zoneVol(d_i,1:clc_nZones)')/sum(r_zoneVol(d_i,1:clc_nZones));
        out_refFlow{d_i,4}(1,d_noise)=sum(out_refTracer{d_i,1}(1:clc_nZones,1))/(d_conc/1000000);
        
        %% 5. PFT Multi - Noise
        out_refTracer{d_i,5}(1:clc_nZones^2,1)=0;
        for d_zone=1:clc_nZones 
            out_refTracer{d_i,5}((d_zone-1)*clc_nZones+d_zone,d_noise)=(r_releaseRateT{d_i}(d_zone)/2)*(randn(1)*0.08+1);
            out_refConc{d_i,5}((d_zone-1)*clc_nZones+1:d_zone*clc_nZones,d_noise)=(d_flow\out_refTracer{d_i,5}((d_zone-1)*clc_nZones+1:d_zone*clc_nZones,d_noise)*1000000).*(randn(clc_nZones,1)*0.08+1);
        end
        d_Q=[];
        for d_zone=1:clc_nZones
            d_Q(d_zone,1:clc_nZones)=(out_refConc{d_i,5}((d_zone-1)*clc_nZones+1:d_zone*clc_nZones,d_noise)/1000000)/out_refTracer{d_i,2}((d_zone-1)*clc_nZones+d_zone,1);
        end
        out_refFlow{d_i,5}(1:clc_nZones^2,d_noise)=-reshape(inv(d_Q)',clc_nZones^2,1);

        %% 6. Const Conc - Noise
        out_refConc{d_i,6}(1:clc_nZones,d_noise)=ones(clc_nZones,1)*(1000/1000000).*(randn(clc_nZones,1)*0.05+1);
        out_refTracer{d_i,6}(1:clc_nZones,d_noise)=(d_flow*(out_refConc{d_i,6}(1:clc_nZones,d_noise))).*(randn(clc_nZones,1)*0.05+1);
        out_refFlow{d_i,6}(1:clc_nZones,d_noise)=out_refTracer{d_i,6}(1:clc_nZones,d_noise)./out_refConc{d_i,3}(1:clc_nZones,1);
    end
end

for ref_bPerm = 1:length(r_flowSim)
    clc_nZones=r_nZones(ref_bPerm);
    clc_flowResults=[];

    clc_simFlow=r_flowSim{ref_bPerm};
    clc_flow=out_refFlow(d_i,:);

    %% External - Total
    for d_refFlow=1:6
        d_3DFlow=[];
        d_3DFlow1=[];
        
        d_3DFlow(:,1)=-sum(clc_simFlow);
        d_3DFlow1{d_refFlow}(:,1)=d_3DFlow;
        d_nNoise=size(clc_flow{1,d_refFlow},2);
        for d_noise=1:d_nNoise
            switch d_refFlow
                case {1,4} ; d_3DFlow1(:,1+d_noise)=clc_flow{1,d_refFlow}(1,d_noise);
                case {2,5}; d_3DFlow1(:,1+d_noise)=-sum(clc_flow{1,d_refFlow}(:,d_noise));
                case {3,6}; d_3DFlow1(:,1+d_noise)=sum(clc_flow{1,d_refFlow}(:,d_noise));
            end
        end

        clc_flowResults{1}{d_refFlow}{1,1}=d_3DFlow1;
    end
    outB_aFlowResults{ref_bPerm}=clc_flowResults;

%     %% External - Exfiltration
%     for d_zone=1:clc_nZones % any(d_refFlow==[2 3 5 6])
%         for d_solve=d_reqSolve
%              for d_seqA=clc_nSeqAverage
%                 if (clc_afType=='S' || clc_afType=='F')
%                     d_seqVlim=clc_nRunSeq-d_seqA;
%                 else
%                     d_seqVlim=1;
%                 end
%                 d_3DFlow=[];
%                 d_3DFlow1=[];
%                 d_3DFlow2=[];
%                 for d_imp=d_reqImp
%                     if (d_imp==3 && d_seqA>1)
%                         continue;
%                     end
%                     for d_conc=d_reqConc
%                         d_3DFlow(:,1)=clc_simFlowTime{d_seqA}(1:d_seqVlim,:);
%                         d_3DFlow(:,2)=-sum(clc_simFlow{d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones),2);
%                         d_3DFlow1{d_imp,d_conc}(:,1:2)=d_3DFlow;
%                         d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA},3);
%                         for d_noise=1:d_nNoise
%                             d_3DFlow1{d_imp,d_conc}(:,2+d_noise)=-sum(clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones,d_noise),2);
%                         end
% 
%                         if (d_seqA>1)
%                             d_3DFlow2{d_imp,d_conc}(:,1:2)=d_3DFlow;
%                             d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA},3);
%                             for d_noise=1:d_nNoise
%                                 d_3DFlow2{d_imp,d_conc}(:,2+d_noise)=-sum(clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+1:d_zone*clc_nZones,d_noise),2);
%                             end
%                         end
% 
%                     end
%                 end
% 
%                 clc_flowResults{d_solve,2}{1,d_zone}{1,find(clc_nSeqAverage==d_seqA)}=d_3DFlow1;
%                 if (d_seqA>1)
%                     clc_flowResults{d_solve,2}{1,d_zone}{2,find(clc_nSeqAverage==d_seqA)}=d_3DFlow2;
%                 end
%             end
%         end
%     end
% 
%     %% External - Infiltration
%     for d_zone=1:clc_nZones
%         for d_solve=d_reqSolve
%              for d_seqA=clc_nSeqAverage
%                 if (clc_afType=='S' || clc_afType=='F')
%                     d_seqVlim=clc_nRunSeq-d_seqA;
%                 else
%                     d_seqVlim=1;
%                 end
%                 d_3DFlow=[];
%                 d_3DFlow1=[];
%                 d_3DFlow2=[];
%                 for d_imp=d_reqImp
%                     if (d_imp==3 && d_seqA>1)
%                         continue;
%                     end
%                     for d_conc=d_reqConc
%                         d_3DFlow(:,1)=clc_simFlowTime{d_seqA}(1:d_seqVlim,:);
%                         d_3DFlow(:,2)=-clc_simFlow{d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+d_zone);
%                         for d_k=1:clc_nZones
%                             if (d_k~=d_zone) 
%                                 d_3DFlow(:,2)=d_3DFlow(:,2)-clc_simFlow{d_seqA}(1:d_seqVlim,(d_k-1)*clc_nZones+d_zone);
%                             end
%                         end
%                         d_3DFlow1{d_imp,d_conc}(:,1:2)=d_3DFlow;
%                         d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA},3);
%                         for d_noise=1:d_nNoise
%                             d_3DFlow1{d_imp,d_conc}(:,2+d_noise)=-clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+d_zone,d_noise);
%                             for d_k=1:clc_nZones
%                                 if (d_k~=d_zone) 
%                                     d_3DFlow1{d_imp,d_conc}(:,2+d_noise)=d_3DFlow1{d_imp,d_conc}(:,2+d_noise)-clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA}(1:d_seqVlim,(d_k-1)*clc_nZones+d_zone,d_noise);
%                                 end
%                             end
%                         end
% 
%                         if (d_seqA>1)
%                             d_3DFlow2{d_imp,d_conc}(:,1:2)=d_3DFlow;
%                             d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA},3);
%                             for d_noise=1:d_nNoise
%                                 d_3DFlow2{d_imp,d_conc}(:,2+d_noise)=-clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA}(1:d_seqVlim,(d_zone-1)*clc_nZones+d_zone,d_noise);
%                                 for d_k=1:clc_nZones
%                                     if (d_k~=d_zone) 
%                                         d_3DFlow2{d_imp,d_conc}(:,2+d_noise)=d_3DFlow2{d_imp,d_conc}(:,2+d_noise)-clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA}(1:d_seqVlim,(d_k-1)*clc_nZones+d_zone,d_noise);
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
% 
%                 clc_flowResults{d_solve,2}{2,d_zone}{1,find(clc_nSeqAverage==d_seqA)}=d_3DFlow1;
%                 if (d_seqA>1)
%                     clc_flowResults{d_solve,2}{2,d_zone}{2,find(clc_nSeqAverage==d_seqA)}=d_3DFlow2;
%                 end
%             end
%         end
%     end
% 
%     %% Internal - Total
%     for d_solve=d_reqSolve
%         for d_seqA=clc_nSeqAverage
%             if (clc_afType=='S' || clc_afType=='F')
%                 d_seqVlim=clc_nRunSeq-d_seqA;
%             else
%                 d_seqVlim=1;
%             end
%             d_3DFlow=[];
%             d_3DFlow1=[];
%             d_3DFlow2=[];
%             for d_imp=d_reqImp
%                 if (d_imp==3 && d_seqA>1)
%                     continue;
%                 end
%                 for d_conc=d_reqConc
%                     d_3DFlow(:,1)=clc_simFlowTime{d_seqA}(1:d_seqVlim,:);
%                     d_3DFlow(:,2)=sum(clc_simFlow{d_seqA}(1:d_seqVlim,:),2);
%                     for d_k=1:clc_nZones 
%                         d_3DFlow(:,2)=d_3DFlow(:,2)-clc_simFlow{d_seqA}(1:d_seqVlim,(d_k-1)*clc_nZones+d_k);
%                     end
%                     d_3DFlow1{d_imp,d_conc}(:,1:2)=d_3DFlow;
%                     d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA},3);
%                     for d_noise=1:d_nNoise
%                         d_3DFlow1{d_imp,d_conc}(:,2+d_noise)=sum(clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA}(1:d_seqVlim,:,d_noise),2);
%                         for d_k=1:clc_nZones 
%                             d_3DFlow1{d_imp,d_conc}(:,2+d_noise)=d_3DFlow1{d_imp,d_conc}(:,2+d_noise)-clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA}(1:d_seqVlim,(d_k-1)*clc_nZones+d_k,d_noise);
%                         end
%                     end
% 
%                     if (d_seqA>1)
%                         d_3DFlow2{d_imp,d_conc}(:,1:2)=d_3DFlow;
%                         d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA},3);
%                         for d_noise=1:d_nNoise
%                             d_3DFlow2{d_imp,d_conc}(:,2+d_noise)=sum(clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA}(1:d_seqVlim,:,d_noise),2);
%                             for d_k=1:clc_nZones 
%                                 d_3DFlow2{d_imp,d_conc}(:,2+d_noise)=d_3DFlow2{d_imp,d_conc}(:,2+d_noise)-clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA}(1:d_seqVlim,(d_k-1)*clc_nZones+d_k,d_noise);
%                             end
%                         end
%                     end
% 
%                 end
%             end
% 
%             clc_flowResults{d_solve,3}{1}{1,find(clc_nSeqAverage==d_seqA)}=d_3DFlow1;
%             if (d_seqA>1)
%                 clc_flowResults{d_solve,3}{1}{2,find(clc_nSeqAverage==d_seqA)}=d_3DFlow2;
%             end
%         end
%     end
% 
%     %% Internal - Individual
%     for d_zone1=1:clc_nZones
%         for d_zone2=1:clc_nZones-1
%             if (d_zone2>=d_zone1)
%                 d_zone2n=d_zone2+1;
%             else
%                 d_zone2n=d_zone2;
%             end
%             for d_solve=d_reqSolve
%                 for d_seqA=clc_nSeqAverage
%                     if (clc_afType=='S' || clc_afType=='F')
%                         d_seqVlim=clc_nRunSeq-d_seqA;
%                     else
%                         d_seqVlim=1;
%                     end
%                     d_3DFlow=[];
%                     d_3DFlow1=[];
%                     d_3DFlow2=[];
%                     for d_imp=d_reqImp
%                         if (d_imp==3 && d_seqA>1)
%                             continue;
%                         end
%                         for d_conc=d_reqConc
%                             d_3DFlow(:,1)=clc_simFlowTime{d_seqA}(1:d_seqVlim,:);
%                             d_3DFlow(:,2)=clc_simFlow{d_seqA}(1:d_seqVlim,(d_zone1-1)*clc_nZones+d_zone2n);
%                             d_3DFlow1{d_imp,d_conc}(:,1:2)=d_3DFlow;
%                             d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA},3);
%                             for d_noise=1:d_nNoise
%                                 d_3DFlow1{d_imp,d_conc}(:,2+d_noise)=clc_flow{d_solve}{d_imp,d_conc}{1,d_seqA}(1:d_seqVlim,(d_zone1-1)*clc_nZones+d_zone2n,d_noise);
%                             end
% 
%                             if (d_seqA>1)
%                                 d_3DFlow2{d_imp,d_conc}(:,1:2)=d_3DFlow;
%                                 d_nNoise=size(clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA},3);
%                                 for d_noise=1:d_nNoise
%                                     d_3DFlow2{d_imp,d_conc}(:,2+d_noise)=clc_flow{d_solve}{d_imp,d_conc}{2,d_seqA}(1:d_seqVlim,(d_zone1-1)*clc_nZones+d_zone2n,d_noise);
%                                 end
%                             end
%                         end
%                     end
%                     clc_flowResults{d_solve,4}{d_zone1,d_zone2}{1,find(clc_nSeqAverage==d_seqA)}=d_3DFlow1;
%                     if (d_seqA>1)
%                         clc_flowResults{d_solve,4}{d_zone1,d_zone2}{2,find(clc_nSeqAverage==d_seqA)}=d_3DFlow2;
%                     end
%                 end
%             end
%         end
%     end

    outB_aFlowResults{ref_bPerm}=clc_flowResults;
end

% try
%     mat_outP3.out_aFlowResults(1,d_batchRun)=outB_aFlowResults;
% catch
%     outB_aFlowResults{d_batchSize}=[];
%     mat_outP3.out_aFlowResults(1,d_batchRun)=outB_aFlowResults;
% end



d_procTime=toc
% mat_outP3.d_procTimePFT=d_procTime;