d_req_i=req_i;

nGroup(1:5)=1;
group=[];
d_countTotal=1;
for i=1:length(groupDims)
    nGroup(i)=length(d_req_i{groupDims(i)});
    d_countTotal=d_countTotal*nGroup(i);
    group{i,1}=d_req_i{groupDims(i)};
end

d_count=1;
missingData=0;

clear out_group out_summary out_summaryHist;

cd(d_dir{d_figAve});

disp(['Group Count Total: ' num2str(d_countTotal)]);

for g1=1:nGroup(1)
    d_req_i{groupDims(1)}=group{1,1}(g1);
    for g2=1:nGroup(2)
        if (length(groupDims)>=2), d_req_i{groupDims(2)}=group{2,1}(g2); end
        for g3=1:nGroup(3)
            if (length(groupDims)>=3), d_req_i{groupDims(3)}=group{3,1}(g3); end
            for g4=1:nGroup(4)
                if (length(groupDims)>=4), d_req_i{groupDims(4)}=group{4,1}(g4); end
                for g5=1:nGroup(5)
                    if (length(groupDims)>=5), d_req_i{groupDims(5)}=group{5,1}(g5); end
                    d_summary=[];
                    %%% Start looping through permutations
                    for i1=d_req_i{1}
                        for i2=d_req_i{2}
                            for i3=d_req_i{3}
                                for i4=d_req_i{4}
                                    for i5=d_req_i{5}
                                        for i6=d_req_i{6}
                                            for i7=d_req_i{7}
                                                for i8=d_req_i{8}
                                                    for i9=d_req_i{9}
                                                        d_input=[];
                                                        try
                                                            outB_aFlowRef
                                                            d_inputWeight=cat(1,out_results{:,1});
                                                            if (i9==2)
                                                                d_selection=find(d_inputWeight(:,3)>d_inputWeight(:,2));
                                                                d_inputWeight=d_inputWeight(d_selection,6);
                                                                d_inputError=cat(1,out_results{d_selection,2});
                                                            else
                                                                d_inputWeight=d_inputWeight(:,6);
                                                                d_inputError=cat(1,out_results{:,2});
                                                            end
                                                            if (ismember(i8,d_concNoNoise)), d_noise=1;
                                                            else d_noise=100; end
                                                            d_inputWeight=reshape(repmat(d_inputWeight,1,d_noise)',[],1);
                                                            d_input=[d_inputError d_inputWeight];
                                                            d_summary=[d_summary; d_input];
                                                        catch
                                                            missingData=missingData+1;
                                                            continue;
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
                    out_group{g1,g2,g3,g4,g5}=d_summary;
                    d_count=d_count+1;
                end
            end
        end
    end
end

out_summary=[];

for g1=1:size(out_group,1)
    for g2=1:size(out_group,2)
        for g3=1:size(out_group,3)
            for g4=1:size(out_group,4)
                for g5=1:size(out_group,5)
                    d_process=out_group{g1,g2,g3,g4,g5};
                    if isempty(d_process)
                        continue
                    end
                    d_process=d_process(d_process(:,2)>0,:);
                    out_summary(end+1,1)=g1;
                    out_summary(end,2)=g2;
                    out_summary(end,3)=g3;
                    out_summary(end,4)=g4;
                    out_summary(end,5)=g5;
                    
                    out_summary(end,6)=size(d_process,1);
                    out_summary(end,8)=wmean(d_process(:,1),d_process(:,2));
                    out_summary(end,7)=out_summary(end,8)-wstd(d_process(:,1),d_process(:,2));
                    out_summary(end,9)=out_summary(end,8)+wstd(d_process(:,1),d_process(:,2));
                    out_summary(end,10)=wprctile(d_process(:,1),5,d_process(:,2));
                    out_summary(end,11)=wprctile(d_process(:,1),15.87,d_process(:,2));
                    out_summary(end,12)=wprctile(d_process(:,1),25,d_process(:,2));
                    out_summary(end,13)=wprctile(d_process(:,1),50,d_process(:,2));
                    out_summary(end,14)=wprctile(d_process(:,1),75,d_process(:,2));
                    out_summary(end,15)=wprctile(d_process(:,1),84.13,d_process(:,2));
                    out_summary(end,16)=wprctile(d_process(:,1),95,d_process(:,2));
                    out_summary(end,17)=out_summary(end,8);
                    out_summary(end,18)=out_summary(end,13);               
                    out_summary(end,19)=wstd(d_process(:,1),d_process(:,2));
                    out_summary(end,20)=out_summary(end,13)-out_summary(end,11);
                    out_summary(end,21)=out_summary(end,15)-out_summary(end,13);
                    out_summary(end,22)=sum(d_process(:,2)); % Check that weights all sum to 1 for each case.
%                     out_summary(end,22)=kurtosis(d_process(:,1));
%                     [out_summaryHist{g1,g2,g3,g4}(:,1) out_summaryHist{g1,g2,g3,g4}(:,2)]=histwc(d_process(:,1),d_process(:,2),[-1:0.05:1]);
                end
            end
        end
    end
end

disp(['Missing Data: ' num2str(missingData)]);

cd ..
        
        
        
        
        
        





