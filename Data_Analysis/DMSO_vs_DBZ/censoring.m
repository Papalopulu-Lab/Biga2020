fr=60/15;
censor=NaN*ones(4,6);
% DBZ1 14022018p0
censor(1,1)=8*fr;%8h
censor(1,2)=2*fr;
censor(1,3)=10*fr;
censor(1,4)=10*fr;
censor(1,5)=10*fr;
censor(1,6)=10*fr;
% DBZ220112017p0
censor(2,1)=6*fr;
censor(2,2)=NaN;
censor(2,3)=NaN;%none
censor(2,4)=10*fr;
censor(2,5)=8*fr;
censor(2,6)=6*fr;
% DBZ30012018p1
censor(3,1)=NaN;
censor(3,2)=NaN;
censor(3,3)=8*fr;
censor(3,4)=NaN;
% DBZ30012018p3
censor(4,1)=12*fr;
censor(4,2)=NaN;
censor(4,3)=6*fr;
censor(4,4)=12*fr;
save acf_nonsignificant