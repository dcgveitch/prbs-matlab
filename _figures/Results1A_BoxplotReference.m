%% Plot boxplots of general grouped results
%% Flow type, Seq Period, Seq Length

% Make boxplots of the non-subsmapled theoretical results
% Interested only in no-noise with differnet combinations of zones

clear
[d_upperPath, d_folder, ~] = fileparts(pwd);
if d_folder(2)=='_', d_folderTS=d_folder(5:15);
else d_folderTS=d_folder(1:11); end

cd Results;
load(strcat(d_folderTS(1:11), '_setup.mat'), '-regexp', '^(?!r_flowSim)...');

nX=4;
nXpad=2;
nY=1;
figPosition=[5, 5, 15, 7];

sizeXMOff=0.05;
sizeYMOff=0.05;
sizeXTotal=0.94;
sizeYTotal=0.85;
sizeXPad=0.01;
sizeYPad=0.01;

sizeXM=(sizeXTotal-((nX-1)*sizeXPad)-(nX/nXpad-1)*sizeXPad)/nX; 
sizeYM=(sizeYTotal-(nY-1)*sizeYPad)/nY;

% What's included in the summary
% If selected as a grouping dimension, the summary is split for each
rPrefix='R';
req_i(1)={1}; % Not Used
req_i(2)={1}; % Not Used
req_i(3)={1:3}; % nZones EXCLUDES 8 ZONES
req_i(4)={1:3}; % Reference Method
req_i(5)={1}; % Not Used
req_i(6)={1}; % Not Used
req_i(7)={1}; % Not Used
req_i(8)={[7]}; % conc
req_i(9)={[1 2 3 4]}; % flowType

% Select grouping dimensions
% 1st=Subplot 2nd=XAxis 3rd=Group
groupDims=[9 3 4];
colours=pmkmp(max(length(req_i{groupDims(3)})+1,3),'CubicL');

d_dir{1}='P6';
d_figAve=1;

fileConc={'Theory' 'LTI' 'Noise' 'Sensor'};
d_concNoNoise=[1 2];
% fileConc={'Theory' '' '' '' 'LTI' 'Sensor' 'Noise' 'Reverse'};
% d_concNoNoise=[1 5];

fileZones={'2' '3' '5' '8' 'All'};

fileDescript=unique(r_nZones);
fileSaveName1=['SImpNoiseOnly.pdf'];
Fig_Boxplot;

% for d_out1=1:3 %nZones
%     for d_out2=[1] %Conc
%         if d_out1<5, req_i(3)={d_out1};
%         else req_i(3)={1:length(unique(r_nZones))};
%         end
%         req_i{8}=d_out2;
%         fileDescript=unique(r_nZones);
%         fileSaveName1=[fileConc{d_out2} '_' fileZones{d_out1} 'Zones.pdf'];
% %         fileSaveName2=['31bit_' fileConc{d_out2} '_' num2str(fileDescript(req_i{3})) 'Zones.pdf'];
%         Fig_Boxplot;
% %         Fig_15bitPerf
%     end
% end