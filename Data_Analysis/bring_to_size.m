function vect=bring_to_size(vect,rs,val)
%% resizes a vector by padding the missing elements with val
%% rs=desired size
%% val= value used for padding ex: -1, 0,NaN
s=size(vect);
if ((s(1)~=rs(1))|(s(2)~=rs(2)))
        if (s(1)>rs(1)) % oversized
            vect=vect(1:rs(1),:);
        else
            vect=[vect;val*ones(rs(1)-s(1),s(2))];
        end
        
        ss=size(vect);
        if (s(2)>rs(2)) %undersized
            vect=vect(:,1:rs(2));
        else
            vect=[vect,val*ones(ss(1),rs(2)-s(2))];
        end
       
end
    