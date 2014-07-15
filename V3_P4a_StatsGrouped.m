% Process Simulink results

% Read test description
clear;
[d_upperPath, d_folderTS, ~] = fileparts(pwd);

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'));
load(strcat(d_folderTS(1:11), '__results.mat'));
mat_prbsConc=matfile(strcat(d_folderTS(1:11), '__prbsConc.mat'),'Writable',true);
outB_prbsConcDisc=mat_prbsConc.out_prbsConcDisc;

clear out_aFlowError*;
    
d_nSolve=[1];
d_impulseType=[1];
d_concType=[3];

d_grouped=1;

if (d_grouped==1)
    id_1=r_seqLength;
    id_2=r_seqPeriod;
    id_3=r_nZones;

    id_group1=unique(id_1);
    id_group2=unique(id_2);
    id_group3=unique(id_3);

    d_graph=cell(length(id_group1),length(id_group2),length(id_group3));

    for i=1:length(id_1)    
        d_pos1=find(id_group1==id_1(i));
        d_pos2=find(id_group2==id_2(i));
        d_pos3=find(id_group3==id_3(i));

        d_graph{d_pos1,d_pos2,d_pos3}(end+1,2)=out_aFlowResults(i);
        d_graph{d_pos1,d_pos2,d_pos3}{end,1}=i;
    end

    % Individual lines sum across?
    % Results stored in which categories?
    
    d_count1=1;
    
    for d_pos2=1:size(d_graph,2)
        for d_pos1=1:size(d_graph,1)
            d_flowError=[];
            d_errorSummary=[];
            d_count2=ones(1,3);
            for d_pos3=1:size(d_graph,3)
                disp(['Group ' num2str(d_count1) ]);
                d_summary=d_graph{d_pos1,d_pos2,d_pos3};
                if (isempty(d_summary))
                    continue;
                end                
                for d_perm=1:size(d_summary,1)
                    disp(['Processing Test ' num2str(d_perm) '/' num2str(size(d_summary,1))]);
                    for d_solve=d_nSolve
                        for d_conc=d_concType
                            for d_impulse=d_impulseType
                                d_seqAList=1:length(r_nSeqAverage{d_summary{d_perm,1}});
                                for d_seqA=d_seqAList
                                    for d_flowType=1:3 % Separate on each output line
                                        switch d_flowType
                                            case 1
                                                d_flowProcess=d_summary{d_perm,2}{d_solve,d_flowType}(1,1);
                                            case 2
                                                d_flowProcess=d_summary{d_perm,2}{d_solve,d_flowType}(1,:);
                                            case 3
                                                d_flowProcess=d_summary{d_perm,2}{d_solve,d_flowType}(1:end,:);
                                        end
                                        d_nFlow1=size(d_flowProcess,1);
                                        d_nFlow2=size(d_flowProcess,2);
                                        d_nNoise=size(d_flowProcess{1,1}{1,d_seqA}{d_impulse,d_conc},2)-2;
                                        d_ndt=size(d_flowProcess{1,1}{1,d_seqA}{d_impulse,d_conc},1);

                                        for d_noise=1:d_nNoise
                                            for d_dt=1:d_ndt
                                                d_flowTotal=0;
                                                for d_flow1=1:d_nFlow1
                                                    for d_flow2=1:d_nFlow2
                                                        d_flowTotal=d_flowTotal+abs(d_flowProcess{d_flow1,d_flow2}{1,d_seqA}{d_impulse,d_conc}(d_dt,2));
                                                    end
                                                end               

                                                for d_flow1=1:d_nFlow1
                                                    for d_flow2=1:d_nFlow2
                                                        d_flowProcessP=d_flowProcess{d_flow1,d_flow2}{1,d_seqA}{d_impulse,d_conc};
                                                        d_flowError{d_flowType}(d_count2(d_flowType),1)=(d_flowProcessP(d_dt,d_noise+2)-d_flowProcessP(d_dt,2))/d_flowProcessP(d_dt,2);
                                                        d_flowError{d_flowType}(d_count2(d_flowType),2)=abs(d_flowProcessP(d_dt,2))/d_flowTotal/d_nNoise;
                                                        d_count2(d_flowType)=d_count2(d_flowType)+1;
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
            end
            for d_flowType=1:3
                % Weighted statistics
                d_errorSummary(1,(d_flowType-1)*7+1)=size(d_flowError{d_flowType},1);
                d_errorSummary(1,(d_flowType-1)*7+2)=wmean(d_flowError{d_flowType}(:,1),d_flowError{d_flowType}(:,2));
                d_errorSummary(1,(d_flowType-1)*7+3)=wprctile(d_flowError{d_flowType}(:,1),5,d_flowError{d_flowType}(:,2));
                d_errorSummary(1,(d_flowType-1)*7+4)=wprctile(d_flowError{d_flowType}(:,1),25,d_flowError{d_flowType}(:,2));
                d_errorSummary(1,(d_flowType-1)*7+5)=wprctile(d_flowError{d_flowType}(:,1),50,d_flowError{d_flowType}(:,2));
                d_errorSummary(1,(d_flowType-1)*7+6)=wprctile(d_flowError{d_flowType}(:,1),75,d_flowError{d_flowType}(:,2));
                d_errorSummary(1,(d_flowType-1)*7+7)=wprctile(d_flowError{d_flowType}(:,1),95,d_flowError{d_flowType}(:,2));

                % Non-weighted statistics
                d_errorSummary(2,(d_flowType-1)*7+1)=size(d_flowError{d_flowType},1);
                d_errorSummary(2,(d_flowType-1)*7+2)=mean(d_flowError{d_flowType}(:,1));
                d_errorSummary(2,(d_flowType-1)*7+3)=prctile(d_flowError{d_flowType}(:,1),5);
                d_errorSummary(2,(d_flowType-1)*7+4)=prctile(d_flowError{d_flowType}(:,1),25);
                d_errorSummary(2,(d_flowType-1)*7+5)=prctile(d_flowError{d_flowType}(:,1),50);
                d_errorSummary(2,(d_flowType-1)*7+6)=prctile(d_flowError{d_flowType}(:,1),75);
                d_errorSummary(2,(d_flowType-1)*7+7)=prctile(d_flowError{d_flowType}(:,1),95); 
            end
            out_aFlowErrorW(d_count1,1)=d_pos1;
            out_aFlowErrorW(d_count1,2)=d_pos2;
            out_aFlowErrorW(d_count1,3)=d_pos3;
            out_aFlowErrorW(d_count1,4:24)=d_errorSummary(1,:);

            out_aFlowErrorNW(d_count1,1)=d_pos1;
            out_aFlowErrorNW(d_count1,2)=d_pos2;
            out_aFlowErrorNW(d_count1,3)=d_pos3;
            out_aFlowErrorNW(d_count1,4:24)=d_errorSummary(2,:);
            d_count1=d_count1+1;
        end
    end
end

%% Save output
save(strcat(d_folderTS(1:11), '__Results.mat'), 'out_*','-v7.3');
cd ..;

