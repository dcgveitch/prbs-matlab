
Ccell=cell(clc_nZones,clc_nZones);
Ccell(:,:)={zeros(length(clc_impAve(:,1)),1)};
Dcell=cell(clc_nZones,1);
Dcell(:,:)={zeros(length(clc_impAve(:,1)),1)};

for i=1:clc_nZones
    Ccell{i,d_zone}=-clc_impAve(:,(i-1)*clc_nZones+d_zone)/clc_zoneVol(d_zone);
    for k=1:clc_nZones
        if (k~=d_zone)
            Ccell{i,k}=clc_impAve(:,(i-1)*clc_nZones+k)/clc_zoneVol(d_zone);
        end
    end
    Dcell{i,1}=clc_impdt(:,(i-1)*clc_nZones+d_zone);
end

C=cell2mat(Ccell);
D=cell2mat(Dcell);            