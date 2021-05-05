close all hidden, clear all, clc 
exptlist={'040417_p1_kymo','161216_p5_kymo','04082016_p2_kymo','26072016_p4_kymo','28032017_p1_kymo'}; 
scaling{1} = 0.581; % change this to the pixel to um scaling factor
scaling{2}=0.488;
scaling{3}=0.488;
scaling{4}=0.254;
scaling{5}=0.488;
figopt=0;
sz{1}=20;
sz{2}=20
sz{3}=20;
for n=1:numel(exptlist)
    M=[]; S=[]; K=[]; L=[];
    %clear cell head
    thisdir=exptlist{n}; % each has LHS and RHS and an ouput folder
    %% analyse content of LHS
    l1=MultiSortFileNames([thisdir,'\LHS']); 
    %% analyse content of RHS
    l2=MultiSortFileNames([thisdir,'\RHS']);
    %merge the two lists
    L=[l1 l2];
    %% analyse each kymo file and store results
    SMatApical=[];
    SMatMedium=[];
    SMatBasal=[];
    HeaderApical=[];
    HeaderMedium=[];
    HeaderBasal=[];
    for i=1:numel(L) 
        tif=[L(i).dirname,'\',L(i).filename];
        if strcmp(L(i).width,'15um') | strcmp(L(i).width,'15umwide')
            side=L(i).dirname(end-2:end);
            % store into correct matrix
            switch L(i).position
                case 'apical'
                     % analyse the kymo
                    mat=AnalyseBandsFixedSize(tif,L(i),scaling{n},figopt,sz{1});
                    % include position in DV
                    for j=1:size(mat,1)
                        header{j}=[side,'_', L(i).z,'_' num2str(j)];
                    end
                    SMatApical=[SMatApical; mat];
                    HeaderApical=[HeaderApical header];
                case 'medium'
                    mat=AnalyseBandsFixedSize(tif,L(i),scaling{n},figopt,sz{2});
                    % include position in DV
                    for j=1:size(mat,1)
                        header{j}=[side,'_', L(i).z,'_' num2str(j)];
                    end
                    SMatMedium=[SMatMedium; mat];
                    HeaderMedium=[HeaderMedium header];
                case 'basal'
                    mat=AnalyseBandsFixedSize(tif,L(i),scaling{n},figopt,sz{3});
                    % include position in DV
                    for j=1:size(mat,1)
                        header{j}=[side,'_', L(i).z,'_' num2str(j)];
                    end
                    SMatBasal=[SMatBasal; mat];
                    HeaderBasal=[HeaderBasal header];
            end
        end
        clear header
    end
    % export the data to file
    warning off
    dr=strtok(exptlist{n},'kymo');
    xlswrite([exptlist{n}, '\Output\' dr 'FixedBandInt.xls'],SMatApical','Apical')
    xlswrite([exptlist{n}, '\Output\' dr 'FixedBandInt.xls'],HeaderApical,'HeadApical')
    xlswrite([exptlist{n}, '\Output\' dr 'FixedBandInt.xls'],SMatMedium','Medium')
    xlswrite([exptlist{n}, '\Output\' dr 'FixedBandInt.xls'],HeaderMedium,'HeadMedium')
    xlswrite([exptlist{n}, '\Output\' dr 'FixedBandInt.xls'],SMatBasal','Basal')
    xlswrite([exptlist{n}, '\Output\' dr 'FixedBandInt.xls'],HeaderBasal,'HeadBasal')
end
%%
   
    
