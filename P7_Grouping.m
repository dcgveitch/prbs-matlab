req_i{1}=d_reqSeqLength;
req_i{2}=d_reqSeqPeriod;
req_i{3}=d_reqNZones;
req_i{4}=d_reqSolve;
req_i{5}=d_reqImp;
req_i{6}=d_reqNSeqA;
req_i{7}=d_reqTSeqA;
req_i{8}=d_reqConc;
req_i{9}=d_reqFlowType;

nGroup(1:5)=1;
group=[];
for i=1:length(groupDims)
    nGroup(i)=length(req_i{groupDims(i)});
    group{i,1}=req_i{groupDims(i)};
end

clear out_group out_summary out_summaryHist;

for g1=1:nGroup(1)
    req_i{groupDims(1)}=group{1,1}(g1);
    for g2=1:nGroup(2)
        if (length(groupDims)>=2), req_i{groupDims(2)}=group{2,1}(g2); end
        for g3=1:nGroup(3)
            if (length(groupDims)>=3), req_i{groupDims(3)}=group{3,1}(g3); end
            for g4=1:nGroup(4)
                if (length(groupDims)>=4), req_i{groupDims(4)}=group{4,1}(g4); end
                for g5=1:nGroup(5)
                    if (length(groupDims)>=5), req_i{groupDims(5)}=group{5,1}(g5); end
                    d_summary=[];
                    %%% Start looping through permutations
                    for i1=req_i{1}
                        for i2=req_i{2}
                            for i3=req_i{3}
                                for i4=req_i{4}
                                    for i5=req_i{5}
                                        for i6=req_i{6}
                                            for i7=req_i{7}
                                                for i8=req_i{8}
                                                    for i9=req_i{9}
                                                        try                                                    
                                                            d_summary=[d_summary; cell2mat(mat_outP6.out_resultsCombined(i1,i2,i3,i4,i5,i6,i7,i8,i9))];
                                                        catch
                                                            try
                                                                d_summary=cell2mat(mat_outP6.out_resultsCombined(i1,i2,i3,i4,i5,i6,i7,i8,i9));
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
                    end
                    out_group{g1,g2,g3,g4,g5}=d_summary;  
                end
            end
        end
    end
end

out_summary=[];

if d_fig==4
    group3Val=[15 30 60];
    group4Val=[1 2 4 8];
    group5Val=[2 4 8 16];
end

for g1=1:size(out_group,1)
    for g2=1:size(out_group,2)
        for g3=1:size(out_group,3)
            for g4=1:size(out_group,4)
                for g5=1:size(out_group,5)
                    d_process=out_group{g1,g2,g3,g4,g5};
                    if isempty(d_process)
%                         out_summary(end,6:16)=0;
                        continue
                    end
                    out_summary(end+1,1)=g1;
                    out_summary(end,2)=g2;
                    out_summary(end,3)=g3;
                    out_summary(end,4)=g4;
                    out_summary(end,5)=g5;
                    
                    out_summary(end,6)=size(d_process,1);
                    out_summary(end,8)=wmean(d_process(:,7),d_process(:,8));
                    out_summary(end,7)=out_summary(end,8)-wstd(d_process(:,7),d_process(:,8));
                    out_summary(end,9)=out_summary(end,8)+wstd(d_process(:,7),d_process(:,8));
                    out_summary(end,10)=wprctile(d_process(:,7),5,d_process(:,8));
                    out_summary(end,11)=wprctile(d_process(:,7),15.87,d_process(:,8));
                    out_summary(end,12)=wprctile(d_process(:,7),25,d_process(:,8));
                    out_summary(end,13)=wprctile(d_process(:,7),50,d_process(:,8));
                    out_summary(end,14)=wprctile(d_process(:,7),75,d_process(:,8));
                    out_summary(end,15)=wprctile(d_process(:,7),84.13,d_process(:,8));
                    out_summary(end,16)=wprctile(d_process(:,7),95,d_process(:,8));
                    if d_fig==4 
                        out_summary(end,17)=group4Val(g4)*group5Val(g5);
                        out_summary(end,18)=group4Val(g4)/group3Val(g3);
                    end
                    [out_summaryHist{g1,g2,g3,g4}(:,1) out_summaryHist{g1,g2,g3,g4}(:,2)]=histwc(d_process(:,7),d_process(:,8),[-1:0.05:1]);
                end
            end
        end
    end
end
        
        
        
        
        
        





