function [struct]=MultiSortFileNames(dirname)
l=ls(dirname);
% get image names inside folder
not_tifs=[];
for j=1:size(l,1)
    try
        tmp=imfinfo([dirname,'\',l(j,:)]);
    catch
        not_tifs=[not_tifs,j];
    end
end
l(not_tifs,:)=[];
% identify z, width, and position 
for j=1:size(l,1)
    strname=l(j,:);
    struct(j).filename=strname;
    struct(j).dirname=dirname;
    [tok,rem]=strtok(strname,'z');
    struct(j).header=tok;
    [tok,rem]=strtok(rem,'_');
    struct(j).z=tok;
    [tok,rem]=strtok(rem,'_');
    struct(j).width=tok;
    [tok,rem]=strtok(rem,'_');
    [tok,tmp]=strtok(tok,'.');
    struct(j).position=tok;
  
end