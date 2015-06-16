function F = vecObj(x,d_zone,nZones,zoneVol,impAve,impdt)
zones=1:nZones;
zones(zones==d_zone)=[];
if nZones==2
    F = [...
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1))-x(d_zone)*impAve(:,d_zone))-impdt(:,d_zone);
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1)+nZones)-x(d_zone)*impAve(:,d_zone+nZones))-impdt(:,d_zone+nZones)];
elseif nZones==3
     F = [...
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1))+x(zones(2))*impAve(:,zones(2))-x(d_zone)*impAve(:,d_zone))-impdt(:,d_zone);
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1)+nZones)+x(zones(2))*impAve(:,zones(2)+nZones)-x(d_zone)*impAve(:,d_zone+nZones))-impdt(:,d_zone+nZones);
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1)+2*nZones)+x(zones(2))*impAve(:,zones(2)+2*nZones)-x(d_zone)*impAve(:,d_zone+2*nZones))-impdt(:,d_zone+2*nZones)];
elseif nZones==5
     F = [...
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1))+x(zones(2))*impAve(:,zones(2))+x(zones(3))*impAve(:,zones(3))+x(zones(4))*impAve(:,zones(4))-x(d_zone)*impAve(:,d_zone))-impdt(:,d_zone);
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1)+1*nZones)+x(zones(2))*impAve(:,zones(2)+1*nZones)+x(zones(3))*impAve(:,zones(3)+1*nZones)+x(zones(4))*impAve(:,zones(4)+1*nZones)-x(d_zone)*impAve(:,d_zone+1*nZones))-impdt(:,d_zone+1*nZones);
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1)+2*nZones)+x(zones(2))*impAve(:,zones(2)+2*nZones)+x(zones(3))*impAve(:,zones(3)+2*nZones)+x(zones(4))*impAve(:,zones(4)+2*nZones)-x(d_zone)*impAve(:,d_zone+2*nZones))-impdt(:,d_zone+2*nZones);
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1)+3*nZones)+x(zones(2))*impAve(:,zones(2)+3*nZones)+x(zones(3))*impAve(:,zones(3)+3*nZones)+x(zones(4))*impAve(:,zones(4)+3*nZones)-x(d_zone)*impAve(:,d_zone+3*nZones))-impdt(:,d_zone+3*nZones);
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1)+4*nZones)+x(zones(2))*impAve(:,zones(2)+4*nZones)+x(zones(3))*impAve(:,zones(3)+4*nZones)+x(zones(4))*impAve(:,zones(4)+4*nZones)-x(d_zone)*impAve(:,d_zone+4*nZones))-impdt(:,d_zone+4*nZones)];
elseif nZones==8
     F = [...
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1))+x(zones(2))*impAve(:,zones(2))+x(zones(3))*impAve(:,zones(3))+x(zones(4))*impAve(:,zones(4))+x(zones(5))*impAve(:,zones(5))+x(zones(6))*impAve(:,zones(6))+x(zones(7))*impAve(:,zones(7))-x(d_zone)*impAve(:,d_zone))-impdt(:,d_zone);
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1)+1*nZones)+x(zones(2))*impAve(:,zones(2)+1*nZones)+x(zones(3))*impAve(:,zones(3)+1*nZones)+x(zones(4))*impAve(:,zones(4)+1*nZones)+x(zones(5))*impAve(:,zones(5)+1*nZones)+x(zones(6))*impAve(:,zones(6)+1*nZones)+x(zones(7))*impAve(:,zones(7)+1*nZones)-x(d_zone)*impAve(:,d_zone+1*nZones))-impdt(:,d_zone+1*nZones);
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1)+2*nZones)+x(zones(2))*impAve(:,zones(2)+2*nZones)+x(zones(3))*impAve(:,zones(3)+2*nZones)+x(zones(4))*impAve(:,zones(4)+2*nZones)+x(zones(5))*impAve(:,zones(5)+2*nZones)+x(zones(6))*impAve(:,zones(6)+2*nZones)+x(zones(7))*impAve(:,zones(7)+2*nZones)-x(d_zone)*impAve(:,d_zone+2*nZones))-impdt(:,d_zone+2*nZones);
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1)+3*nZones)+x(zones(2))*impAve(:,zones(2)+3*nZones)+x(zones(3))*impAve(:,zones(3)+3*nZones)+x(zones(4))*impAve(:,zones(4)+3*nZones)+x(zones(5))*impAve(:,zones(5)+3*nZones)+x(zones(6))*impAve(:,zones(6)+3*nZones)+x(zones(7))*impAve(:,zones(7)+3*nZones)-x(d_zone)*impAve(:,d_zone+3*nZones))-impdt(:,d_zone+3*nZones);
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1)+4*nZones)+x(zones(2))*impAve(:,zones(2)+4*nZones)+x(zones(3))*impAve(:,zones(3)+4*nZones)+x(zones(4))*impAve(:,zones(4)+4*nZones)+x(zones(5))*impAve(:,zones(5)+4*nZones)+x(zones(6))*impAve(:,zones(6)+4*nZones)+x(zones(7))*impAve(:,zones(7)+4*nZones)-x(d_zone)*impAve(:,d_zone+4*nZones))-impdt(:,d_zone+4*nZones);
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1)+5*nZones)+x(zones(2))*impAve(:,zones(2)+5*nZones)+x(zones(3))*impAve(:,zones(3)+5*nZones)+x(zones(4))*impAve(:,zones(4)+5*nZones)+x(zones(5))*impAve(:,zones(5)+5*nZones)+x(zones(6))*impAve(:,zones(6)+5*nZones)+x(zones(7))*impAve(:,zones(7)+5*nZones)-x(d_zone)*impAve(:,d_zone+5*nZones))-impdt(:,d_zone+5*nZones);
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1)+6*nZones)+x(zones(2))*impAve(:,zones(2)+6*nZones)+x(zones(3))*impAve(:,zones(3)+6*nZones)+x(zones(4))*impAve(:,zones(4)+6*nZones)+x(zones(5))*impAve(:,zones(5)+6*nZones)+x(zones(6))*impAve(:,zones(6)+6*nZones)+x(zones(7))*impAve(:,zones(7)+6*nZones)-x(d_zone)*impAve(:,d_zone+6*nZones))-impdt(:,d_zone+6*nZones);
        1/zoneVol(d_zone)*(x(zones(1))*impAve(:,zones(1)+7*nZones)+x(zones(2))*impAve(:,zones(2)+7*nZones)+x(zones(3))*impAve(:,zones(3)+7*nZones)+x(zones(4))*impAve(:,zones(4)+7*nZones)+x(zones(5))*impAve(:,zones(5)+7*nZones)+x(zones(6))*impAve(:,zones(6)+7*nZones)+x(zones(7))*impAve(:,zones(7)+7*nZones)-x(d_zone)*impAve(:,d_zone+7*nZones))-impdt(:,d_zone+7*nZones)];   
end



