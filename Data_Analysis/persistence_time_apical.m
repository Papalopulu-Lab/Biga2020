clear all, close all hidden, clc
% please run Kymo_trackbandsfixedsize.m before this to generate 20um regs
exptlist={'040417_p1_kymo','161216_p5_kymo','04082016_p2_kymo','26072016_p4_kymo','28032017_p1_kymo'}; 
TH=[]; TL=[];
for n=1:numel(exptlist)
    header=strtok(exptlist{n},'kymo');
    fname=strcat([exptlist{n},'\Output\',header,'FixedBandInt.xls']);
    num1=xlsread(fname,'Apical');
    [~,head1]=xlsread(fname,'HeadApical');
    time=([1:size(num1,2)]-1)/4; % fps is 15min
    data=num1';
    [th tl]=getKymoPersistence(data,time,0);
    TH=[TH; bring_to_size(th,[1,100],NaN)];
    TL=[TL; bring_to_size(tl,[1,100],NaN)];    
end
figure,subplot(1,3,1), hist(TH(:)), title('Time in high state');
xlabel('Time');ylabel('Counts');
subplot(1,3,2),hist(TL(:)),title('Time in low state');
xlabel('Time');ylabel('Counts');
subplot(1,3,3),hist(TH(:)./TL(:)),title('Ratio high vs low');
xlabel('Time');ylabel('Counts');
suptitle('Apical');

