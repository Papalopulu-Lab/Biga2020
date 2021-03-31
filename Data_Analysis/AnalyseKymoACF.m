function [MeanSacrPkInt,SacrPkInt,LSP_period]=AnalyseKymoACF(KymoG,tif,flag,scaling)
% uncoment this to view the kymo data
% figure,imagesc(KymoG);
% xlabel('Ventral<--Space-->Dorsal');
% ylabel('<--Time');
% title('Kymograph')
%% Find all intensity values for time 0 and plot
MeanIntervalALL = [];
IntervalALL = [];
PgmALL =[];
cALL = [];
MedIntervalALL =[];
SpData_raw=[];
SpData_detrend=[];
Spectral_detrend=[];
LSP_period=NaN*ones(size(KymoG,2),2);
KymoG=double(KymoG);
bl=find(KymoG==0);
pvec=1-0.0001;
if ~isempty(bl)
    disp('Blank pixels excluded');
    KymoG(bl)=NaN;
end
n=1;
f=0:0.0010:0.2;
Pthresh=[];
for i=1:8:size(KymoG,1)-8 % every 2h
    % can be modified to include every time point 
    % for i=1:size(KymoG,1)
    % t0 = double(KymoG(i,:));
    % every x hours
    t0 = double(KymoG(i:i+8,:)); %(rows are time, columns are space in the kymograph)
    t0=mean(t0,1,'omitnan');
    t0_raw=t0;
    totlength = size(t0);
    dist = [1:totlength(1,2)];
    space = dist * scaling ;
    SpData_raw=[SpData_raw t0_raw(:)]; 
    % figure
    % plot(space,t0,'LineWidth',2);
    % title('Intensity', 'FontSize',14 );
    % xlabel('Distance (um)', 'FontSize',14 );
    % ylabel('Venus::HES5 intensity', 'FontSize',14 );
    %% Detrend intensity
    space(isnan(t0))=[];
    t0(isnan(t0))=[];
    p = polyfit(space,t0,4);
    fit = polyval(p,space);
    % figure
    % plot(space,t0,'b');
    % hold on
    % plot(space,fit,'r--');
    % legend('raw','trend');
    Detrend = t0-fit ;
    Detrend_raw=NaN*ones(size(t0_raw));
    Detrend_raw(~isnan(t0_raw))=Detrend;
    totlength = size(t0_raw);
    dist = [1:totlength(1,2)];
    space_raw = dist * scaling ;
    SpData_detrend=[SpData_detrend Detrend_raw(:)]; 
    space=space*scaling;
    [pxx,fx,pth] = plomb(Detrend,space,'Pd',pvec);
    Pthresh=[Pthresh pth];
    pxx1 = plomb(Detrend,1/(space(2)-space(1)),f);
    % export and save LSP output
    Spectral_detrend=[Spectral_detrend bring_to_size(pxx1,[size(KymoG,2)*2,1],NaN)];
    % find dominant 1st and 2nd dominant peaks 
    pxx1=pxx;
    idx1=find(pxx1==max(pxx1),1,'first');
    pxx1(idx1)=-Inf;
    idx2=find(pxx1==max(pxx1),1,'first');
    if pxx(idx1)>pth
        LSP_period(n,1)=1/fx(idx1);
    end
    if pxx(idx2)>pth
        LSP_period(n,2)=1/fx(idx2);
    end
    %Find autocorrelation
    [c,lags] = xcorr(Detrend);
    smacr=c;
    Corrlags = lags * scaling;
    if numel(Detrend)==size(KymoG,2)
        Corrlags_save=Corrlags;
        f_save=f;
    end
    maxlagno=size(KymoG,2)-numel(Detrend);
    cALL = [cALL; [NaN*ones(1,maxlagno), c, NaN*ones(1,maxlagno)]]; % 987
    [pks,locs] = findpeaks(smacr,'minpeakheight',0.01*max(c),'minpeakdistance',10/scaling);
    x_peaks = Corrlags(locs);
    % find peak to peak distances
    PeakInterval = diff(x_peaks);
    MeanInterval = mean(PeakInterval);
    MeanIntervalALL = [MeanIntervalALL; bring_to_size(MeanInterval,[1,size(KymoG,1)],NaN)]; 
    IntervalALL = [IntervalALL; bring_to_size(PeakInterval,[1,size(KymoG,1)],NaN)]; 
    MedInterval=median(PeakInterval);
    MedIntervalALL = [MedIntervalALL; bring_to_size(MedInterval,[1,size(KymoG,1)],NaN)];
    n=n+1;
end
Corrlags=Corrlags_save;
%% analyse individually and save the detected peaks
ACFPeaks_ALL=[];
for i=1:size(cALL,1)
    BootMat=[];
    for k=1:100
        % bootstrap
        kidx=randperm(size(SpData_detrend,1));
        randvec=SpData_detrend(kidx,i);
        tmp=randvec; tmp(isnan(tmp))=[];
        randvec(isnan(randvec))=mean(tmp); % exclude NaNs
        [cboot,lag]=xcorr(randvec);
        BootMat=[BootMat cboot(:)];
    end
    sd=std(BootMat');
    m=-2*sd; m(lag==0)=NaN;
    M=2*sd; M(lag==0)=NaN;
    % end of bootstrap
    Sacr =cALL(i,:);
    Ssmacr = Sacr; 
    [Sacrpks,Sacrlocs] = findpeaks(Ssmacr,'minpeakheight',0.01*max(Sacr),'minpeakdistance',10/scaling); %  minpeakheight 1% of max ACF
    Sacrx_peaks = Corrlags(Sacrlocs);
%     % keep only signficant peaks
    for j=1:numel(Sacrlocs)
        if Sacrpks(j)<M(Sacrlocs(j))
            Sacrx_peaks(j)=NaN;
            Sacrlocs(j)=NaN;
        end
    end
    disp('excluded')
    sum(isnan(Sacrx_peaks))/numel(Sacrx_peaks)
    Sacrx_peaks(isnan(Sacrx_peaks))=[];
    ACFPeaks_ALL=[ACFPeaks_ALL;bring_to_size(Sacrlocs,[1,50],NaN)];
    % find peak to peak distances
    SacrPeakInterval = diff(Sacrx_peaks);
    SacrPkInt(i) = mean(SacrPeakInterval) % mean per ACF
end
SacrPkInt(isnan(SacrPkInt))=[];
MeanSacrPkInt=mean(SacrPkInt); % mean per slice
disp ('Saving results...')
save(strtok(tif,'.'));