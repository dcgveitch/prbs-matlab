%% Read test description
clear;
%#ok<*FNDSB>
tic

[d_upperPath, d_folder, ~] = fileparts(pwd);
if d_folder(2)=='_', d_folderTS=d_folder(5:15);
else d_folderTS=d_folder(1:11); end

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'), '-regexp', '^(?!r_flowSim)...');
cd P6;
load(strcat(d_folderTS(1:11), '__outP6Single.mat'));
mat_outP6=matfile(strcat(d_folderTS(1:11), '__outP6Single.mat'),'Writable', true);


for d_solve=d_reqSolve
    for d_impulse=d_reqImp
        for d_conc=d_reqConc
            for d_flowType=1:4
                d_process=out_resultsCombined{d_flowType,d_conc}{d_solve,d_impulse};
                for d_i=unique(d_process(:,4))'
                    d_processFilter=d_process(find(d_process(:,4)==d_i),:);
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,1)=size(d_processFilter,1);
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,3)=wmean(d_processFilter(:,2),d_processFilter(:,1));
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,2)=out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,3)-wstd(d_processFilter(:,2),d_processFilter(:,1));
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,4)=out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,3)+wstd(d_processFilter(:,2),d_processFilter(:,1));
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,5)=wprctile(d_processFilter(:,2),5,d_processFilter(:,1));
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,6)=wprctile(d_processFilter(:,2),15.87,d_processFilter(:,1));
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,7)=wprctile(d_processFilter(:,2),25,d_processFilter(:,1));
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,8)=wprctile(d_processFilter(:,2),50,d_processFilter(:,1));
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,9)=wprctile(d_processFilter(:,2),75,d_processFilter(:,1));
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,10)=wprctile(d_processFilter(:,2),84.13,d_processFilter(:,1));
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,11)=wprctile(d_processFilter(:,2),95,d_processFilter(:,1));
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,12)=out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,3);
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,13)=out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,8);               
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,14)=wstd(d_processFilter(:,2),d_processFilter(:,1));
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,15)=out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,8)-out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,6);
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,16)=out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,10)-out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,8);
                    out_resultsCombinedSummary{d_flowType,d_conc}{d_solve,d_impulse}(d_i,17)=sum(d_processFilter(:,1)); % Check that weights all sum to 1 for each case.
                end
            end
        end
    end
end       

mat_outP6.out_resultsCombinedSummary=out_resultsCombinedSummary;

