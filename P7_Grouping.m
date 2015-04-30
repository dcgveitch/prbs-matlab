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

clear out_group out_summary out_summaryHist;

cd(d_dir{d_figAve});

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
                    disp(['Group ' num2str(d_count) '/' num2str(d_countTotal)]);
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
                                                            filename=[num2str(i1) '_' num2str(i2) '_' num2str(i3) '_' num2str(i4) '_' num2str(i5) '_' num2str(i6) '_' num2str(i7) '_' num2str(i8) '_' num2str(i9) '.mat'];
                                                            load(filename);
                                                            d_inputError=cat(1,out_results{:,2});
                                                            d_inputWeight=cat(1,out_results{:,1});
                                                            d_inputWeight=d_inputWeight(:,5);
                                                            if (ismember(i8,d_concNoNoise)), d_noise=1;
                                                            else d_noise=100; end
                                                            d_inputWeight=reshape(repmat(d_inputWeight,1,d_noise)',[],1);
                                                            d_input=[d_inputError d_inputWeight];
                                                            d_summary=[d_summary; d_input];
                                                        catch
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
                    out_summary(end,17)=kurtosis(d_process(:,1));
%                     [out_summaryHist{g1,g2,g3,g4}(:,1) out_summaryHist{g1,g2,g3,g4}(:,2)]=histwc(d_process(:,1),d_process(:,2),[-1:0.05:1]);
                end
            end
        end
    end
end

cd ..
        
        
        
        
        
        





