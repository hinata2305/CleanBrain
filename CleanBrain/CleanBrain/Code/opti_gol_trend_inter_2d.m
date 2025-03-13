%Copyright © 2022 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

% This function finds optimal detrend SG filter for a given timecourse. 

% Input: this function needs 1 denoised time course set (Pred) 1 grey matter timecourse set (time) 
% 1 noise time course set (noise). Upper and lower boundaries of the SG window search space under study (windowa windowb). 
% Upper and lower boundaries for the SG polynomial order search space under study (ordera orderb). 
% 2 files that describe the beginning of a psychological event (expera experb) and file that describes the timecourse section of cognitive interest (part).  

% Output: Correlation between predictor and observed timecourse for a given SG filter.


function [ opti ] = opti_gol_trend_inter_2d(time,noise,pred,windowa,windowb,ordera,orderb,experA,experB,part)

% Create predictor timecourse which is a concatenation of event related averages
pred_time=extract_section_single(pred,experB,part);

underwindow=windowa;
upperwindow=round(((windowb-windowa)./2)-1);

opti = zeros(upperwindow,orderb-ordera+1,size(time,2));

for ii = 1:upperwindow
    
    span = underwindow+((ii-1)*2);
    
    for jj = ordera:min(orderb,ii-1)        
     
     order=jj;
     % SG filtered timecourse.
     trend=detrend_gol(time,span,order,1);  
     % Remove noise and trend from timecourse.
     cleantrend=remove_noise(time,noise,trend);
     % Remove jitter from time course.
     cleantrend_time=extract_section_single(cleantrend,experA,part);
     % Correlate predicted and observed timecourse.
     opti(ii,jj,:)=diag(corr([cleantrend_time.section],[pred_time.pred]));
    
     
    end
    
end

end

