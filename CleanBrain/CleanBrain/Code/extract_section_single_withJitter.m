%Copyright © 2022 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

% Input:
% - timecourse: A matrix representing the time course data to be analyzed.
% - encode1: A vector indicating the time points (TRs) where events start (encoded as 1).
% - part: A 2-element vector specifying the time duration to be included before 
%         and after the event onset (defined in TR units).
% - ende1: A vector indicating the end points of the analysis for each event.

% Output:
% - predictorWithJitter: A matrix containing the modified time course with 
%   event-related averages and jitter applied.

function [predictorWithJitter] = extract_section_single_withJitter(timecourse, encode1, part, ende1)

% Extract the relevant section of the time course based on event timing
timecourse1 = timecourse(encode1(1):ende1(1), :);

% Get the number of events and define the trial length
l = length(encode1); 

% Extract time intervals for data selection
part1 = part(1);
part2 = part(2);

% Adjust encoding to start from the first event
encode1 = encode1 - (encode1(1) - 1);

% Calculate last sample of each event and create a shift vector
last = encode1 + part(2);
shift = encode1;
shift(1) = [];
shift = shift - 1;
last(end) = [];

% Initialize a temporary storage for the selected time segments
for i = 1:l
    % Define the start and end points for the current trial
    a1 = encode1(i) + part1;
    e1 = encode1(i) + part2;
    
    % Extract the relevant time course segment for the current trial
    temptime1(i, :, :) = timecourse1((a1:e1), :);
end

% Determine the dimensions of the extracted segments
stemp = size(temptime1);

% Calculate event-related averages for each time course segment
for i = 1:stemp(3)
    % Squeeze the current trial data to one dimension
    tempie1 = squeeze(temptime1(:, :, i));
    % Compute the average across trials for the current event
    event_av(:, i) = mean(tempie1, 1);  % Average across the first dimension 
end    

% Update the time course with event-related averages
for i = 1:stemp(1)
    a1 = encode1(i) + part1;
    e1 = encode1(i) + part2;
    timecourse1(a1:e1, :) = event_av;  % Replace the original segment with averages
end

% Introduce jitter into the time course based on the shift values
for i = 1:stemp(1) - 1
    temp = last(i);
    tempval = timecourse1(temp, :);        % Retrieve the latest average value
    m = (shift(i) - last(i));              % Calculate the number of rows to be filled
    tempmat = tempval .* (ones(m, 34));    % Expand the values across the required rows
    anfang = temp + 1;                     % Start point for insertion
    ende = shift(i);                       % End point for insertion
    timecourse1(anfang:ende, :) = tempmat; % Insert the expanded values into the time course
    clear tempmat temp;                     % Clear temporary variables to free memory
end

% Return the modified time course with jitter
predictorWithJitter = timecourse1;

end
