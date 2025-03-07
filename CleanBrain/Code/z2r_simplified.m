
%Copyright © 2022 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 
%This function back transforms correlation data using tanh in the context of a fisher ztransform. 
%It brings the cell data of the single subjects into a common matrix depending on the dimension of data and renames data.

function [mat, dat] = z2r_simplified(biggerdata)

l=length(biggerdata);

bigdata=cell2mat(biggerdata);
namie=fieldnames(bigdata);
newname = strrep(namie, '_trans', '');
corstruc=strcmp(newname,namie);
ln=length(namie);


%for the number of vars in biggerdata
for i = 1 : ln
    
    temp=[bigdata.(namie{i})];
    s=size(temp);
    
    % if the matrix is an connectivity image
    if s(2)> l
       
        %replace nan with 0
        index=isnan(temp)
        temp(index)=0;
        
        % reshape the 2d format into 3d
        temp_3d=reshape(temp,s(1),s(1),l);
        
        if corstruc(i)==0
        % calc mean of temp_3d along the z axis = 3 dimension
        ztransie=tanh(nanmean(temp_3d,3));
        else
        
        ztransie=(nanmean(temp_3d,3));
            
        end
        
        temp=[];
        
        for k = 1 : (s(1)-1)
            temp=[(ztransie(k+1:end,k));temp];
        end
        
        
        dat.(newname{i})=temp;
        mat.(newname{i})=ztransie;
        
        
    else
        
       if corstruc(i)==0
           
         dat.(newname{i})=tanh(temp);
           
           
       else
           
         dat.(newname{i})= (temp);
           
       end
       
    end
    
end

