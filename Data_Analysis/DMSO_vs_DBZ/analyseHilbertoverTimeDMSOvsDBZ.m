%% Importing kymograph tiffs
close all, clear all, clc 
preanalysed=1; % 1 means has been analysed
dmso_exptlist={'p1_DMSO_kymo','p3_DMSO_kymo','supersand_DMSO_kymo'};
dmso_scaling = 0.594; 
KOmat=[];
maxT=10*60/15;
for n=1:numel(dmso_exptlist)
    M=[]; S=[]; K=[]; L=[];
    clear cell head
    thisdir=dmso_exptlist{n}; % each has LHS and RHS and an ouput folder
    %% analyse content of LHS
    l1=MultiSortFileNames([thisdir,'\LHS']); 
    %% analyse content of RHS
    l2=MultiSortFileNames([thisdir,'\RHS']);
    %merge the two lists
    L=[l1 l2];
    for i=1:numel(L)
        tmp=L(i).filename; tmp=strtok(tmp,'.');
        Kymo=imread([L(i).dirname '\' tmp],'tiff');
        KymoG=rgb2gray(Kymo);
        % detrend the kymo signal
        SpData_detrend=[];
        for j=1:size(KymoG,1)
            t0 = double(KymoG(j,:));
            totlength = size(t0);
            dist = [1:totlength(1,2)];
            space = dist * dmso_scaling;
            p = polyfit(space,t0,4); % poly order 4
            fit = polyval(p,space);
            Detrend = t0-fit ;
            SpData_detrend=[SpData_detrend Detrend(:)];
        end
        s=size(SpData_detrend);
        % get hilbert phase reconstruction
        PhaseMap=[];
        for k=1:maxT
            dvect1=SpData_detrend(:,k);
            dvect=dvect1;
            for l=1:50
                dvect=sgolayfilt(dvect,3,5); 
            end
            smooth=dvect;
            %
            hb=hilbert(smooth);
            amp=abs(hb);
            ph=angle(hb);
            % store ph into a matrix
            PhaseMap=[PhaseMap; ph'];
        end
        % calculate phase synchronisation in the same region
        KO=abs(sum(exp(1i*PhaseMap))/s(1));
        kom=mean(KO);
        KOmat=[KOmat, kom];
    end
end
%% additional exclusions
KOmat([10,11])=[];
figure,subplot(1,2,1),boxplot(KOmat)
set(gca,'XTickLabel','DMSO');
set(gca,'YLim',[0,0.15]);
ylabel('Phase synchronisation');
disp('DMSO Mean')
mean(KOmat)
%% analyse DBZ
dbz_exptlist={'14022018_p3_2umDBZ_kymo','20112017_2umDBZ_kymo','30012018_p1_2umDBZ_kymo','30012018_p3_2umDBZ_kymo'};   
dbz_scaling= 0.59; 
% timepoints excluded based on ACF
censoring
KOmat=[];
for n=1:numel(dbz_exptlist)
    M=[]; S=[]; K=[]; L=[];
    clear cell head
    thisdir=dbz_exptlist{n}; % each has LHS and RHS and an ouput folder
    %% analyse content of LHS
    l1=MultiSortFileNames([thisdir,'\LHS']); 
    %% analyse content of RHS
    l2=MultiSortFileNames([thisdir,'\RHS']);
    %merge the two lists
    L=[l1 l2];
    KOmPerExpt=[];
    for i=1:numel(L)
        tmp=L(i).filename; tmp=strtok(tmp,'.');
        Kymo=imread([L(i).dirname '\' tmp],'tiff');
        KymoG=rgb2gray(Kymo);
        % detrend the kymo signal
        SpData_detrend=[];
        for j=1:size(KymoG,1)
            t0 = double(KymoG(j,:));
            totlength = size(t0);
            dist = [1:totlength(1,2)];
            space = dist * dbz_scaling;
            p = polyfit(space,t0,4); % poly order 4
            fit = polyval(p,space);
            Detrend = t0-fit ;
            SpData_detrend=[SpData_detrend Detrend(:)];
        end
        s=size(SpData_detrend);
        % get hilbert phase reconstruction
        PhaseMap=[];
        if ~isnan(censor(n,i))
            maxT=censor(n,i);
        else
            maxT=s(2);
        end
        for k=1:maxT
            dvect1=SpData_detrend(:,k);
            dvect=dvect1;
            for l=1:50
                dvect=sgolayfilt(dvect,3,5); 
            end
            smooth=dvect;
            %
            hb=hilbert(smooth);
            amp=abs(hb);
            ph=angle(hb);
            % store ph into a matrix
            PhaseMap=[PhaseMap; ph'];
        end
        % calculate phase synchronisation in the same region
        KO=abs(sum(exp(1i*PhaseMap))/s(1));
        kom=mean(KO);
        KOmat=[KOmat, kom];
    end
end
% further exclusions based on ACF
KOmat([2,10,19])=[];
subplot(1,2,2),boxplot(KOmat)
set(gca,'XTickLabel','DBZ');
set(gca,'YLim',[0,0.15]);
ylabel('Phase synchronisation');
disp('DBZ mean')
mean(KOmat)
    
    
