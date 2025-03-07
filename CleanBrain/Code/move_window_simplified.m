function [windowmean, max_move] = move_window_simplified(act, mov, window)
% move_window_simplified computes the mean reliability of head motion sub-samples using a sliding window approach.
% 
% This function processes a time series of reliability values and corresponding head motion parameters 
% to compute the mean reliability and maximum head motion within specified sub-samples (windows).
%
% Inputs:
%   act    - Vector of reliability values (e.g., measures of signal quality or activation).
%   mov    - Vector of head motion parameters over time.
%   window  - Size of the sliding window used to create sub-samples.
%
% Outputs:
%   windowmean - Mean reliability values calculated for each sub-sample.
%   max_move   - Maximum head motion observed within each sub-sample.
%
% Note:
% The function starts by sorting the head motion parameters and then iteratively computes
% the mean reliability and maximum head motion for each moving window until the end of the input vector.
%
% Algorithm Steps:
% 1. Sort the head motion parameters and keep track of the original indices.
% 2. For each valid starting position of the window:
%    - Define the current sub-sample indices based on the specified window size.
%    - Extract the range of head motion values corresponding to the sub-sample.
%    - Compute the maximum head motion within this sub-sample.
%    - Transform the reliability values using atanh to obtain activation sections.
%    - Calculate the mean reliability for the current sub-sample.
%
% This sliding window approach allows for analyzing how reliability changes in relation to head motion over time.

% Sorting the head motion parameters alongside their original indices
[C,I] = sort(mov);

% Iterating through the entire sample using moving windows of the specified sub-sample size
for j = 1 :(length(I) - (window - 1))

    % Determine the starting index of the current sub-sample
    under = j;
    % Determine the ending index of the current sub-sample
    upper = j + (window - 1);
    
    % Extract the indices of the current sub-sample
    samp = I([under:upper]);
    
    % Compute the largest head motion within the current sub-sample
    max_move(j) = max(C([under:upper]));

    % Extract motion parameters for the current sub-sample
    move_section = mov(samp);

    % Transform the reliability values for the current sub-sample
    activation_section = atanh(act(samp));
    
    % Calculate the mean reliability for the current sub-sample
    windowmean(:,j) = tanh(mean(activation_section));
end
end

