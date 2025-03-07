%Copyright © 2022 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

% Function:
% This function applies a Savitzky-Golay filter to a given time course.
% It has two options: type 1 applies an inverted window at the ends of the
% time course, while type 2 applies no such window.
% 
% Input: 
% timec = time course matrix, 
% span = window size, 
% order = polynomial order.
% type = 1 applies an inverted window at the ends of the time course; 
% type = 2 applies no such window.
% 
% Output: 
% detrend_time = a filtered signal.

function [detrend_time] = detrend_gol(timec, span, order, type)

% Size of the time course matrix.
s = size(timec);
sizeroi = s(2);
% Half the window size used for ends.
m = (span - 1) / 2;

if type == 1
    
    if size(span) == 1
    
        % Loop over the columns of the matrix.
        parfor it = 1:sizeroi
            
            % Half the window size at the beginning of the signal.
            begint = flipud(timec(1:m, it));
            % Half the window size at the end of the signal.
            endt = flipud(timec(end - m + 1:end, it));
            % Combine the time course pieces with inverted ends.
            timecWrand = [begint; timec(:, it); endt];
            % Apply the Savitzky-Golay filter.
            trend = smooth(timecWrand, span, 'sgolay', order);
            % Extract the relevant part of the filtered data.
            detrend_time(:, it) = trend(m + 1:(m + s(1)));
            
        end
        
    else
        
        parfor it = 1:sizeroi
            
            % Half the window size at the beginning of the signal.
            begint = flipud(timec(1:m, it));
            % Half the window size at the end of the signal.
            endt = flipud(timec(end - m + 1:end, it));
            % Combine the time course pieces with inverted ends.
            timecWrand = [begint; timec(:, it); endt];
            % Apply the Savitzky-Golay filter with varying span and order.
            trend = smooth(timecWrand, span(it), 'sgolay', order(it));
            % Extract the relevant part of the filtered data.
            detrend_time(:, it) = trend(m + 1:(m + s(1)));
            
        end

    end
     
else
   
    parfor it = 1:sizeroi
        
        % Apply the Savitzky-Golay filter without inverted ends.
        detrend_time(:, it) = smooth(timec(:, it), span, 'sgolay', order);
        
    end
    
end
    