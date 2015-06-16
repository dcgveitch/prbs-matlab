% Process Sensor results

%% Read test description
clear;
%#ok<*FNDSB>
tic

d_reqSolve=[1 2 3 4 5];
d_reqImp=[1 2];

mat_testCell=matfile('TestCell.mat','Writable',true);
d_testName={'T141207' 'T150511' 'T150522' 'T150528'};
d_simFolder='150611T1602_FixPhysical';
d_folderTS=d_simFolder(1:11);

cd(d_simFolder);
cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'))
load(strcat(d_folderTS(1:11), '__outP1.mat'));
load(strcat(d_folderTS(1:11), '__outP2.mat'));

for d_test=1:length(d_testName)
    %% Test Description
    disp(['Test: ' num2str(d_test)])
    clear sens_* clc_* cl_*;
    
    sens_data=d_testName{d_test};
    eval(['sens_results=mat_testCell.' sens_data ';']);
    
    clc_simRef=sens_results.simRef;
    
    out_concTraces{d_test,1}(:,1)=sens_results.z1(:,1)-sens_results.extSim(:,1);
    out_concTraces{d_test,1}(:,2)=sens_results.z1(:,2)-sens_results.extSim(:,1);
    out_concTraces{d_test,2}(:,1)=sens_results.z2(:,1)-sens_results.extSim(:,2);
    out_concTraces{d_test,2}(:,2)=sens_results.z2(:,2)-sens_results.extSim(:,2);
    
    for d_i=1:7
        for d_zone=1:sens_results.nZones
            out_concTraces{d_test,d_zone}(:,2+d_i)=out_prbsConcDisc{d_test}{d_i}(1:length(out_concTraces{d_test,d_zone}(:,2)),d_zone);
        end
    end
end
    

