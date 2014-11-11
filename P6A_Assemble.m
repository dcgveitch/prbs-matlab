% %% Read test description
clear;
%#ok<*FNDSB>
tic

[d_upperPath, d_folder, ~] = fileparts(pwd);
if d_folder(2)=='_', d_folderTS=d_folder(5:15);
else d_folderTS=d_folder(1:11); end

cd Results;
mat_outP6=matfile(strcat(d_folderTS(1:11), '__outP6.mat'),'Writable',true);
d_nBatch=34;
d_count=1;

% for d_i=1:d_nBatch
%     disp(['Batch ' num2str(d_i)]);
%     d_input=mat_outP6.outM_resultsCombined(1,d_i);
%     if(~isempty(d_input{1}))
%         d_inputConCat{d_count}=vertcat(d_input{1}{:});
%         d_count=d_count+1;
%     end
%     clear d_input;
% end

clear d_inputConCat;

for d_i=1:d_nBatch
    disp(['Batch ' num2str(d_i)]);
    d_input=mat_outP6.outM_resultsCombinedRef(1,d_i);
    if(~isempty(d_input{1}))
        d_inputConCat{d_count}=vertcat(d_input{1}{:});
        d_count=d_count+1;
    end
    clear d_input;
end


