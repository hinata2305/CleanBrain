%Copyright © 2022 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

% Input: This function needs a timecourse, a file that defines when the beginning of an event occurs (encode 1), 
% and a file that defines the number of TRs that should be grepped after the beginning of an event (part).

% Output: this function extracts the relevant part of the time course and makes event related averages, predictor time courses, 
% and observed courses without jitter.

function [toi] = extract_section_single(timecourse1, encode1, part)

    % Length of relevant trial starts expressed in Tr.
    l=length(encode1);
        % Timecourse section of interest. 
        part1=part(1);
        part2=part(2);
            
            for i=1:l
                
                % Start of trial.
                a1=encode1(i)+part1;
                % End of trial.
                e1=encode1(i)+part2;
                
                temptime1(i,:,:)=timecourse1((a1:e1),:);
                 
            end
            
            stemp = size(temptime1);
            
            for i=1:stemp(3)
                % 
                tempie1 = squeeze(temptime1(:,:,i));
                % Event related average.
                toi(i).event_av=(mean(tempie1));
                % Empirical predictor function.
                toi(i).pred=(repmat(toi(i).event_av,1,l))';
                stemp1=size(tempie1);
                % Timecourse without jitter.
                toi(i).section=reshape(tempie1',(stemp1(1)*stemp1(2)),1);
               
            end          
end
 