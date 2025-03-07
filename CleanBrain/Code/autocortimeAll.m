%Copyright © 2022 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

% Input: time course
% Output: auto correlation lags 1:4


function [ auto ] = autocortimeAll(time)

s=size(time);
unit=((s(2)+1));
grasp=[1:unit:(s(2)*s(2))];
[r_ww,~] = xcorr(time,4,'coeff');
auto=r_ww(1:4,grasp);

end

