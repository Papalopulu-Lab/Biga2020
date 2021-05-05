function Smat=AnalyseBandsFixedSize(tif,L,scaling,figopt,sz)
sz=round(sz/scaling);
% analysis has not been done
Kymo=imread(tif,'tiff');
KymoG=rgb2gray(Kymo);
KymoG=imgaussfilt(KymoG,round(2/scaling));
% reorder as space by time
KymoG=flipud(KymoG'); 
% read the first 4 timepoints and get bands using Hilbert
vect0=mean(KymoG(:,1:16)');
dist = [1:numel(vect0)];
space = dist * scaling;
p = polyfit(space,vect0,4);
fit = polyval(p,space);
if figopt==1
    h=figure,imshow(KymoG*2),colormap(jet),title(L.filename,'Interpreter','None');
    xlabel('Time');
    ylabel('Space');
    hold on
    figure,subplot(2,1,1),plot(space,vect0)
    ylabel('Venus::HES5 Intensity');
    xlabel('Position on DV axis (um)');
    hold on
    plot(space,fit,'r--');
    legend('raw','trend');
    suptitle(L.filename);%,'Interpreter','None') 
end
Detrend = vect0-fit ;
Detrend=Detrend-mean(Detrend);
if figopt==1
    subplot(2,1,2),plot(space,Detrend);
    hold on
end
% smooth the data
Detrend=sgolayfilt(Detrend,1,11);
Detrend=sgolayfilt(Detrend,1,11);
% divide into equal sized intervals
idx=1:sz:size(KymoG,1);
if figopt==1
    plot(space,Detrend,'k');
    % plot areas over kymo and export data
    figure(h)
    for j=1:numel(idx)
        plot([1,size(KymoG,2)],[idx(j) idx(j)],'LineWidth',2,'Color',[0 0 0]);
    end
end
% export the average signal over time
Smat=[];
if figopt==1
   h=get(gcf);
   figno=h.Number+2;
   figure(figno);
end
% %% check this output is correct
for j=1:numel(idx)-1
    chunk=KymoG(idx(j):idx(j+1),:);
    sig=mean(chunk);
    if (figopt==1)&& (j<=9)
        subplot(3,3,j),plot(sig);
    end
    Smat=[Smat; sig];
end
if figopt==1
    suptitle(L.filename);%,'Interpreter','None')
end