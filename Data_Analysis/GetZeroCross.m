function [idx]=GetZeroCross(data,thresh)
err=sign(data(1:end-1))+sign(data(2:end));
idx=find(err==0);
% find points that are close together and keep the last one
list=[1:numel(idx)];
label=zeros(1,numel(idx));
while numel(list)>1
    keep=list(1);
    cx=idx(list(1));
    % remove this from list
    list(1)=[];
    dist=sqrt((cx-idx(list)).^2);
    tmp=find(dist<20);
    label(keep)=1;
    list(tmp)=[];
    %pause
end
label(list)=1;
idx=idx(logical(label));
