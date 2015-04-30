for test=1:4
    in_summary=in_data{test};
    clear d_*;
    d_state=[7 8];
    d_decay=[1 2;3 4];
    d_support=[5 6 9 10 11 12];
    d_decayCount=zeros(1,length(d_state));

    for d_i=1:size(in_summary,1)-1
        for d_j=1:length(d_state)
            if (in_summary(d_i+1,d_state(d_j))==0 && in_summary(d_i,d_state(d_j))==1)
                d_decayCount(d_j)=d_decayCount(d_j)+1;
            end
            d_decayMap(d_i,d_j)=d_decayCount(d_j);
        end
    end

    d_count1=1;
    for d_i=1:length(d_state)
        for d_j=1:d_decayCount(d_i)
            for d_k=1:size(d_decay,2)
                d_process=in_summary(find(d_decayMap(:,d_i)==d_j),[d_decay(d_i,d_k) d_support(1) d_support(2) d_support(3) d_support(4) d_support(5) d_support(6)]);
                d_process=d_process(1:find(d_process(:,1)==min(d_process(:,1))),:);
                if (~isempty(find(d_process(:,1)==2001))) 
                    d_process=d_process(find(d_process(:,1)==2001,1,'last')+1:end,:);
                end
                d_process=d_process(find(d_process(:,1)>=1000),:);
                d_outC{d_k,d_count1}=d_process;
            end
            d_count1=d_count1+1;
        end
    end

    t=30/3600;
    
    d_count1=1;
    for d_i=1:size(d_outC,1)
        for d_j=1:size(d_outC,2)
            d_process=d_outC{d_i,d_j};
            d_temp=mean(d_process(:,6:7),2);
            Y=d_process(:,1)-mean(d_process(:,2));
            out_decayTrace{test}(1:length(Y),d_count1)=Y;
            d_count1=d_count1+1;
            X=t:t:size(Y,1)*t;
            X=X';

            % Log linear fit
            [fit,gf]=F_Linear(X,log(Y(1)./Y(:)));   

            temp=coeffvalues(fit);
            d_ACH{d_i}(1,d_j)=temp(1);
            GF_ln(1,d_j)=gf.adjrsquare;
            GF_ln(2,d_j)=gf.sse;
            GF_ln(3,d_j)=gf.rmse;
            % Exponential fit
            [fit,gf]=F_Exp(X,Y(:));
%             plot(fit);
%             hold on
%             plot(X,Y);

            temp=coeffvalues(fit);
            d_ACH{d_i}(2,d_j)=-temp(2);
            GF_exp(1,d_j)=gf.adjrsquare;
            GF_exp(2,d_j)=gf.sse;
            GF_exp(3,d_j)=gf.rmse;
            % Exponential fit
            [fit,gf]=F_ExpOff(X,Y(:));

            temp=coeffvalues(fit);
            d_ACH{d_i}(3,d_j)=temp(2);

            GF_expOff(1,d_j)=gf.adjrsquare;
            GF_expOff(2,d_j)=gf.sse;
            GF_expOff(3,d_j)=gf.rmse;
            % Integral
            d_ACH{d_i}(4,d_j)=(Y(1)-Y(end))/((size(Y,1)*t)*mean(Y(:)));
            % 2 Point
            d_ACH{d_i}(5,d_j)=log(Y(1)/Y(end))/(size(Y,1)*t);

            d_ACH{d_i}(6,d_j)=mean(d_temp); % Average temperature in zone
            d_ACH{d_i}(7,d_j)=mean(d_process(:,4)); % Average Fan1 RPM
            d_ACH{d_i}(8,d_j)=mean(d_process(:,5)); % Average Fan2 RPM
        end
    end
    out_ACH{test}=d_ACH;
    out_decay{test}=d_outC;
end

for d_i=1:length(out_decayTrace)
    for d_j=1:size(out_decayTrace{d_i},2)
        d_process=out_decayTrace{d_i}(:,d_j);
        d_process=d_process./max(d_process);
        out_decayTraceNorm{d_i}(1:length(d_process),d_j)=d_process;
    end
end

        
        
    
    
    
    
    

            
        
    