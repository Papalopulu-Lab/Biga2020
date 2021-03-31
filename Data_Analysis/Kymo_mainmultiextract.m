%% Importing kymograph tiffs
close all, clear all, clc 
preanalysed=0; % 1 means has been analysed
exptlist={'040417_p1_kymo','161216_p5_kymo','04082016_p2_kymo','26072016_p4_kymo','28032017_p1_kymo','28032017_p6_kymo'}; 
scaling{1} = 0.581; % change this to the pixel to um scaling factor
scaling{2}=0.488;
scaling{3}=0.488;
scaling{4}=0.254;
scaling{5}=0.488;
scaling{6}=0.488; 
for n=1:numel(exptlist)
    M=[]; S=[]; K=[]; L=[];
    clear cell head
    thisdir=exptlist{n}; % each has LHS and RHS and an ouput folder
    %% analyse content of LHS
    l1=MultiSortFileNames([thisdir,'\LHS']); 
    %% analyse content of RHS
    l2=MultiSortFileNames([thisdir,'\RHS']);
    %merge the two lists
    L=[l1 l2];
    %% analyse each kymo file and store results
    for i=1:numel(L)
        i
        tif=[L(i).dirname,'\',L(i).filename];
        if preanalysed==0
            % analysis has not been done
            Kymo=imread(tif,'tiff');
            KymoG=rgb2gray(Kymo);
            [M(i),M_ACF,LSP]=AnalyseKymoACF(KymoG,tif,flag,scaling{n});
            K(i).kymo=KymoG;
            LSP_period=LSP(1:size(M_ACF,2),:); 
            M_LSP(i)=mean(LSP_period(:),'omitnan');
        else
            % pre-analysed load the files
            tmp=L(i).filename; tmp=strtok(tmp,'.');
            tif=[L(i).dirname,'\',tmp,'.mat'];
            load(tif,'KymoG','MeanSacrPkInt','spacing','SpData_raw',...
            'SpData_detrend','cALL','Corrlags','ACFPeaks_ALL','space',...
            'LSP_period','Spectral_detrend','f_save');
            time=[0:size(SpData_detrend,2)-1]*0.25;
            % export these variables
            M(i)=MeanSacrPkInt;
            K(i).kymo=KymoG;
            LSP_period=LSP_period(1:size(cALL,1),:);
            M_LSP(i)=mean(LSP_period(:),'omitnan');
            % uncomment to save each file info
            xlswrite([L(i).dirname,'\',tmp],SpData_raw,'IntRaw');
            xlswrite([L(i).dirname,'\',tmp],SpData_detrend,'IntDetrend');
            xlswrite([L(i).dirname,'\',tmp],space(:),'Dist_X');
            xlswrite([L(i).dirname,'\',tmp],cALL','ACF');
            xlswrite([L(i).dirname,'\',tmp],ACFPeaks_ALL','ACF_XPeaks');
            xlswrite([L(i).dirname,'\',tmp],Corrlags(:),'ACF_lag');
            xlswrite([L(i).dirname,'\',tmp],f_save(:),'Frequency');
            xlswrite([L(i).dirname,'\',tmp],Spectral_detrend,'LSP');
            xlswrite([L(i).dirname,'\',tmp],M_LSP,'Period_LSP');
            % extra info -band size per tiff
            INTall=[];
            for j=1:size(ACFPeaks_ALL,1)
                idx=ACFPeaks_ALL(j,:);
                idx(isnan(idx))=[];
                bob=Corrlags(idx);
                int=bob(2:end)-bob(1:end-1);
                int=int(numel(int)/2+1:end);
                int=bring_to_size(int,[1,10],NaN);
                INTall=[INTall int(:)];
            end
            xlswrite([L(i).dirname,'\',tmp],INTall,'Band_size');
        end  
        % split by Apical Medium Basal
        apical_acf(n)=mean(M(1:3:end));
        basal_acf(n)=mean(M(2:3:end));
        medium_acf(n)=mean(M(3:3:end));
        % 
        apical_lsp(n)=mean(M_LSP(1:3:end));
        basal_lsp(n)=mean(M_LSP(2:3:end));
        medium_lsp(n)=mean(M_LSP(3:3:end));
    end
end
%% plot auto-correlation averages per experiment
figure, boxplot([apical_acf(:), medium_acf(:), basal_acf(:)]);
title('Auto-correlation analysis');
ylabel('Peak:Peak distance per experiment');
set(gca,'XTickLabel',{'Apical','Interm.','Basal'});
%% plot lomb-scargle averages per experiment
figure, boxplot([apical_lsp(:), medium_lsp(:), basal_lsp(:)]);
title('Lomb-Scargle periodogram analysis');
ylabel('Spatial periodicity per experiment');
set(gca,'XTickLabel',{'Apical','Interm.','Basal'});

    
    
