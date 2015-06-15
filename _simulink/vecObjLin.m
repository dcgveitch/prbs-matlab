
Ccell=cell(clc_nZones^2,clc_nZones^2);
Ccell(:,:)={zeros(length(clc_impAve(:,1)),1)};
Dcell=cell(clc_nZones^2,1);
Dcell(:,:)={zeros(length(clc_impAve(:,1)),1)};

for i=1:clc_nZones
    for j=1:clc_nZones
            Ccell((i-1)*clc_nZones+j,(i-1)*clc_nZones+1:i*clc_nZones)={-clc_impAve(:,(j-1)*clc_nZones+i)/clc_zoneVol(i)};
            Dcell{(i-1)*clc_nZones+j,1}=clc_impdt(:,(j-1)*clc_nZones+i);
            for k=1:clc_nZones
                if (k~=i)
                    Ccell{(i-1)*clc_nZones+j,(k-1)*clc_nZones+i}=clc_impAve(:,(j-1)*clc_nZones+k)/clc_zoneVol(i);
                end
            end
    end
end

C=cell2mat(Ccell);
D=cell2mat(Dcell);

C=[zeros(size(C,1),clc_nZones) C];
D=[D; zeros(clc_nZones,1)];

Cext=zeros(clc_nZones,clc_nZones*(clc_nZones+1));

for i=1:clc_nZones
    Cext(i,i)=-1;
    Cext(i,i*clc_nZones+1:(i+1)*clc_nZones)=1;
    for j=1:clc_nZones
        if (j~=i)
            Cext(i,j*clc_nZones+i)=-1;
        end
    end
end

C=[C; Cext];




            