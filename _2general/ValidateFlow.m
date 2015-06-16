d_pos=1;
for d_perm=1:setup_nSim
    d_temp=reshape(mean(r_flowSim{d_perm},1),r_nZones(d_perm),r_nZones(d_perm))';
    % Total ACH
    out_flowValTotal(d_perm,1)=-sum(sum(d_temp))/sum(r_zoneVol(d_perm,1:r_nZones(d_perm)));
    % Total internal?
    out_flowValTotal(d_perm,2)=(sum(sum(d_temp))-sum(diag(d_temp)))/sum(r_zoneVol(d_perm,1:r_nZones(d_perm)));
    for d_zone=1:r_nZones(d_perm)
        out_flowValZone(d_pos,1)=r_nZones(d_perm);
        out_flowValZone(d_pos,2)=-sum(d_temp(d_zone,:))/r_zoneVol(d_perm,d_zone);
        out_flowValZone(d_pos,3)=-sum(d_temp(:,d_zone))/r_zoneVol(d_perm,d_zone);
        out_flowValZone(d_pos,4)=(sum(d_temp(d_zone,:))-d_temp(d_zone,d_zone))/r_zoneVol(d_perm,d_zone);
        d_pos=d_pos+1;
    end
end

out_intACH=[];

for d_perm=1:setup_nSim
    for d_i=1:size(r_flowIntSplit{d_perm},2)
        out_intACH=[out_intACH; r_flowIntSplit{d_perm}(:,d_i)];
    end
end
            


testflow(:,1)=-sum(r_flowSim{1}(:,1:2),2);
testflow(:,2)=-sum(r_flowSim{1}(:,3:4),2);
testflow(:,3)=-(r_flowSim{1}(:,1)+r_flowSim{1}(:,3));
testflow(:,4)=-(r_flowSim{1}(:,4)+r_flowSim{1}(:,2));
testflow(:,5)=r_flowSim{1}(:,2);
testflow(:,6)=r_flowSim{1}(:,3);





% for d_perm=1:setup_nSim
%     parfor d_i=1:length(r_flowSim{d_perm})
%         d_temp=reshape(r_flowSim{d_perm}(d_i,:),r_nZones(d_perm),r_nZones(d_perm))';
%         dd_fullFlow=[];
%         for d_zone=1:r_nZones(d_perm)
%             dd_fullFlow((d_zone-1)*3+1)=-sum(d_temp(d_zone,:));
%             dd_fullFlow((d_zone-1)*3+2)=-sum(d_temp(:,d_zone));
%             dd_fullFlow((d_zone-1)*3+3)=(sum(d_temp(d_zone,:))-d_temp(d_zone,d_zone));
%         end
%         d_fullFlow(d_i,:)=dd_fullFlow;
%     end
%     out_fullFlow{d_perm}=d_fullFlow;
% end