ventTC=[1 2 5 10];
ventRate=1./(ventTC.*3600);
sensorResp=[30 60 150 300];
stepSize=[1];
count=1;

for d_sR=sensorResp
    for d_vR=ventRate
        for d_sS=stepSize
            x_in=0:1:1/d_vR*20-1;
            y_in=d_sS*exp(-x_in*d_vR);
            d_a=1/d_sR;
            y_out=filter(d_a, [1 d_a-1], y_in);
            [fit,gf,fitDetails]=F_ExpFix(x_in(10*d_sR:end)*d_vR,y_out(10*d_sR:end));
            temp=coeffvalues(fit);
            out_Results(count,1)=temp(1);
            y_fit=out_Results(count,1)*exp(-x_in*d_vR);
            y_resi=y_fit-y_out;
            out_resi{count}=y_resi;
            out_fitResi{count}=fitDetails.residuals;
            [fit,gf]=F_Exp(x_in,y_resi);
            temp=coeffvalues(fit);
            out_Results(count,2)=temp(1);
            out_Results(count,3)=-temp(2);
            out_Results(count,4)=(1/-temp(2))/d_sR;
            count=count+1;
        end
    end
end