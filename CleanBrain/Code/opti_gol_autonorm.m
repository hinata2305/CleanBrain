%Copyright © 2022 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

% Input:
% This function requires two cleaned grey matter timecourse sets (time, pred) and two
% files that define the timing of psychological events (experA and experB).
% It also needs the 'part' variable which specifies the duration of the timecourse section
% to analyze after the onset of a psychological event.

% Output:
% The function computes the optimal low-pass Savitzky-Golay (SG) filter for a denoised 
% and detrended timecourse by estimating the correlation between the observed timecourse 
% and the predictor. It generates three outputs:
% 1. opti: Correlation values for the optimal filter settings.
% 2. autofull: Autocorrelation of the complete timecourse.
% 3. autosection: Autocorrelation of timecourses excluding jitter segments.

function [ opti, autofull, autosection] = opti_gol_autonorm(time,pred,experA,experB,part)

% Create a predictor function based on the relevant timecourse section for cognitive interest.
% experB indicates the onset of the cognitive event, while part defines the window length 
% for cognitive interest in TR (time resolution).

% Generate the predictor timecourse by extracting the specified section from the prediction data.
pred_time = extract_section_single(pred, experB, part);
predictor = [pred_time.pred];  % Extract the predictor values.

% Initialize matrices to store results for different filter configurations.
opti = zeros(242, 50, size(time, 2));          % Optimal correlation results for different configurations.
autofull = zeros(242, 50, 4, size(time, 2));  % Autocorrelation results for the entire timecourse.
autosection = zeros(242, 50, 4, size(time, 2)); % Autocorrelation results for specific sections.

% Define the start and end of the timecourse based on experA.
start_timecourse = experA(1);       % Start of the timecourse in TR.
end_timecourse = experA(24) + 17;   % End of the timecourse in TR.

% Iterate through a range of window spans to find the optimal settings for the Savitzky-Golay filter.
for ii = 1:242
    % Define the window span for the SG filter.
    span = (ii * 2 + 1);
    parfor jj = 1:min((span - 1), 50)
        % Set the order of the Savitzky-Golay filter.
        order = jj;

        % Apply the detrending filter to the timecourse data.
        cleanfilter = detrend_gol(time, span, order, 1);

        % Calculate the autocorrelation for the entire timecourse.
        autofull(ii, jj, :, :) = autocortimeAll(cleanfilter(start_timecourse:end_timecourse, :));

        % Extract the segment of the timecourse that aligns with the defined cognitive events.
        cleantrend_time = extract_section_single(cleanfilter, experA, part);

        % Compute the autocorrelation for the section without jitter.
        autosection(ii, jj, :, :) = autocortimeAll([cleantrend_time.section]);

        % Calculate the correlation between the observed and predicted timecourses.
        opti(ii, jj, :) = diag(corr([cleantrend_time.section], predictor));
    end
end

end
