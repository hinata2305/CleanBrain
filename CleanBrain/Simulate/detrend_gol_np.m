% Warranty this code has not undergone peer code review yet. We do not take responsibility for its 
% correctness. This code may not be used for medical applications. �Koten Sch�ppen 22/11/2017
% 
%  Function:
%  This function applies a Savitzky-Golay-Filter to a given timecourse
%  It has two options; type 1 = applies an inverted window at the rand of the
%  timecourse. type = 2 applies no such window.
% 
%  Input: 
%  timec=timecourse matrix, span = window size, order=polynomial order
%  type 1 = applies an inverted window at the rand of the
%  timecourse. type = 2 applies no such window.
% 
%  Output: A filtered signal


function [detrend_time] = detrend_gol_np(timec,span,order,type);

% Size of the time course matrix.
s=size(timec);
sizeroi=s(2);
% Half the window size used for rand.
m=(span-1)./2;

if type==1
    
    if size(span)==1
    
    % Loop over matrix.
    for it=1:sizeroi
        
        % Half the window size at the beginning of the signal.
        begint=flipud(timec(1:m,it));
        % Half the window size at the end of the signal.
        endt=flipud(timec(s(1)-m:s(1),it));
        % Put the timecourse pieces together.
        timecWrand=[begint; (timec(:,it)); endt];
        % Apply the filter
        trend = smooth(timecWrand,span,'sgolay',order);
        % Take the relevant part of the filtered data
        detrend_time(:,it)=trend(m+1:(m+s(1)));
        
    end
    
    
    else
        
    for it=1:sizeroi
        
        % Half the window size at the beginning of the signal.
        begint=flipud(timec(1:m,it));
        % Half the window size at the end of the signal.
        endt=flipud(timec(s(1)-m:s(1),it));
        % Put the timecourse pieces together.
        timecWrand=[begint; (timec(:,it)); endt];
        %  Apply the filter.
        trend = smooth(timecWrand,span(it),'sgolay',order(it));
        % Take the relevant part of the filtered data.
        detrend_time(:,it)=trend(m+1:(m+s(1)));
        
    end
     
    end
     
else
   
    for it=1:sizeroi
        
        detrend_time(:,it) = smooth(timec(:,it),span,'sgolay',order);
        
        
    end
    
end

    
    