% Process Simulink results

%% Read test description
clear;
%#ok<*FNDSB>
tic

[d_upperPath, d_folderTS, ~] = fileparts(pwd);

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'), '-regexp', '^(?!r_flowSim)...')
mat_outP3=matfile(strcat(d_folderTS(1:11), '__outP3.mat'),'Writable',true);
mat_outP4=matfile(strcat(d_folderTS(1:11), '__outP4.mat'),'Writable',true);

out_aFlowResults=mat_outP3.out_aFlowResults;
    
d_reqSolve=[1];
d_reqImp=[1 2 3];
d_reqConc=[1 2 3 4];

d_solve=d_reqSolve;
d_impulse=d_reqImp;
d_conc=d_reqConc;
d_seqA=1;

d_grouped=1;

if (d_grouped==1)
    id_1=r_seqLength;
    id_2=r_seqPeriod;

    id_group1=unique(id_1);
    id_group2=unique(id_2);

    d_graph=cell(length(id_group1),length(id_group2));

    for i=1:length(id_1)    
        d_pos1=find(id_group1==id_1(i));
        d_pos2=find(id_group2==id_2(i));

        d_graph{d_pos1,d_pos2}(end+1,2)=out_aFlowResults(i);
        d_graph{d_pos1,d_pos2}{end,1}=i;
    end
    
    for d_impulse=d_reqImp
        for d_conc=d_reqConc
            disp(['Type ' num2str(d_impulse) ':' num2str(d_conc)]);
            d_count1=1;
            for d_pos2=1:size(d_graph,2)
                for d_pos1=1:size(d_graph,1)            
                    d_flowError=[];
                    d_errorSummary=[];
                    d_count2=ones(1,3);
                    disp(['Group ' num2str(d_count1) ]);
                    d_summary=d_graph{d_pos1,d_pos2};
                    if (isempty(d_summary))
                        continue;
                    end                
                    for d_perm=1:size(d_summary,1)
                        disp(['Processing Test ' num2str(d_perm) '/' num2str(size(d_summary,1))]);
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
                    for d_flowType=1:3
                        % Weighted statistics
                        d_errorSummary(1,(d_flowType-1)*11+1)=size(d_flowError{d_flowType},1);
                        d_errorSummary(1,(d_flowType-1)*11+3)=wmean(d_flowError{d_flowType}(:,1),d_flowError{d_flowType}(:,2));
                        d_errorSummary(1,(d_flowType-1)*11+2)=d_errorSummary(1,(d_flowType-1)*11+3)-wstd(d_flowError{d_flowType}(:,1),d_flowError{d_flowType}(:,2));
                        d_errorSummary(1,(d_flowType-1)*11+4)=d_errorSummary(1,(d_flowType-1)*11+3)+wstd(d_flowError{d_flowType}(:,1),d_flowError{d_flowType}(:,2));
                        d_errorSummary(1,(d_flowType-1)*11+5)=wprctile(d_flowError{d_flowType}(:,1),5,d_flowError{d_flowType}(:,2));
                        d_errorSummary(1,(d_flowType-1)*11+6)=wprctile(d_flowError{d_flowType}(:,1),15.87,d_flowError{d_flowType}(:,2));
                        d_errorSummary(1,(d_flowType-1)*11+7)=wprctile(d_flowError{d_flowType}(:,1),25,d_flowError{d_flowType}(:,2));
                        d_errorSummary(1,(d_flowType-1)*11+8)=wprctile(d_flowError{d_flowType}(:,1),50,d_flowError{d_flowType}(:,2));
                        d_errorSummary(1,(d_flowType-1)*11+9)=wprctile(d_flowError{d_flowType}(:,1),75,d_flowError{d_flowType}(:,2));
                        d_errorSummary(1,(d_flowType-1)*11+10)=wprctile(d_flowError{d_flowType}(:,1),84.13,d_flowError{d_flowType}(:,2));
                        d_errorSummary(1,(d_flowType-1)*11+11)=wprctile(d_flowError{d_flowType}(:,1),95,d_flowError{d_flowType}(:,2));
                        [out_aFlowErrorWHist{d_impulse,d_conc}{d_count1,d_flowType}(:,1) out_aFlowErrorWHist{d_impulse,d_conc}{d_count1,d_flowType}(:,2)]=histwc(d_flowError{d_flowType}(:,1),d_flowError{d_flowType}(:,2),[-1:0.05:1]);

                        % Non-weighted statistics
                        d_errorSummary(2,(d_flowType-1)*11+1)=size(d_flowError{d_flowType},1);
                        d_errorSummary(2,(d_flowType-1)*11+3)=mean(d_flowError{d_flowType}(:,1));
                        d_errorSummary(2,(d_flowType-1)*11+2)=d_errorSummary(2,(d_flowType-1)*11+3)-std(d_flowError{d_flowType}(:,1));
                        d_errorSummary(2,(d_flowType-1)*11+4)=d_errorSummary(2,(d_flowType-1)*11+3)+std(d_flowError{d_flowType}(:,1));
                        d_errorSummary(2,(d_flowType-1)*11+5)=prctile(d_flowError{d_flowType}(:,1),5);
                        d_errorSummary(2,(d_flowType-1)*11+6)=prctile(d_flowError{d_flowType}(:,1),15.87);
                        d_errorSummary(2,(d_flowType-1)*11+7)=prctile(d_flowError{d_flowType}(:,1),25);
                        d_errorSummary(2,(d_flowType-1)*11+8)=prctile(d_flowError{d_flowType}(:,1),50);
                        d_errorSummary(2,(d_flowType-1)*11+9)=prctile(d_flowError{d_flowType}(:,1),75);
                        d_errorSummary(2,(d_flowType-1)*11+10)=prctile(d_flowError{d_flowType}(:,1),84.13);
                        d_errorSummary(2,(d_flowType-1)*11+11)=prctile(d_flowError{d_flowType}(:,1),95);
                        [out_aFlowErrorNWHist{d_impulse,d_conc}{d_count1,d_flowType}(:,1) out_aFlowErrorNWHist{d_impulse,d_conc}{d_count1,d_flowType}(:,2)]=hist(d_flowError{d_flowType}(:,1),[-1:0.05:1]);
                    end
                    out_aFlowErrorW{d_impulse,d_conc}(d_count1,1)=id_group1(d_pos1);
                    out_aFlowErrorW{d_impulse,d_conc}(d_count1,2)=id_group2(d_pos2);
                    out_aFlowErrorW{d_impulse,d_conc}(d_count1,3:35)=d_errorSummary(1,:);

                    out_aFlowErrorNW{d_impulse,d_conc}(d_count1,1)=id_group1(d_pos1);
                    out_aFlowErrorNW{d_impulse,d_conc}(d_count1,2)=id_group2(d_pos2);
                    out_aFlowErrorNW{d_impulse,d_conc}(d_count1,3:35)=d_errorSummary(2,:);
                    d_count1=d_count1+1;
                end
            end
        end
    end
end

%% Save output
mat_outP4.out_aFlowErrorW=out_aFlowErrorW;
mat_outP4.out_aFlowErrorNW=out_aFlowErrorNW;
mat_outP4.out_aFlowErrorWHist=out_aFlowErrorWHist;
mat_outP4.out_aFlowErrorNWHist=out_aFlowErrorNWHist;
cd ..;

