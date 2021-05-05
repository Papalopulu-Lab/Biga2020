function [TH TL]=getKymoPersistence(datamat,time,figopt)
TH=[];
TL=[];
for n=1:size(datamat,2)
    vect0=datamat(:,n);
    Detrend=vect0-mean(vect0);
    % smooth the zero mean data
    Detrend=sgolayfilt(Detrend,1,11);
    Detrend=sgolayfilt(Detrend,1,11);
    idx=GetZeroCross(Detrend);
    if figopt==1
        figure
        plot(time,Detrend,'k');
        hold on
        plot(time(idx),Detrend(idx),'*')
    end
    switch numel(idx)
        case 1
            int1=time(idx)-time(1);
            int2=time(end)-time(idx-1);
            if sum(Detrend(1:idx)>=0) % first int is high
                t_h=int1;
                t_l=int2;
            else
                t_h=int2;
                t_l=int1;
            end
        otherwise
            if idx(1)>15
                idx=[1;idx];
            end
            if numel(vect0)-idx(end)>15
                idx=[idx;numel(vect0)];
            end
            % interval size
            dt=time(idx(2:end))-time(idx(1:end-1));
            % split into high and low
            if sum(Detrend(idx(1):idx(2))>=0) % first int is high
                t_h=dt(1:2:end);
                t_l=dt(2:2:end);
            else
                t_h=dt(2:2:end);
                t_l=dt(1:2:end);
            end
            t_h=mean(t_h);
            t_l=mean(t_l);
    end
    TH=[TH t_h];
    TL=[TL t_l];
end