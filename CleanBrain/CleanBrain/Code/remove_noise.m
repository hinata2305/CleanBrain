%Copyright © 2022 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

% Description: This function removes noise from a time series data. It utilizes MATLAB's regress function
% to perform linear regression and effectively eliminate noise components. 

function [clean] = remove_noise(timeraw, noise, detrenddata, highfreqnoise)

    % Get the size of the raw time series data
    s = size(timeraw);
    
    switch nargin  % Determine the number of input arguments to handle different cases
        
        case 4  % If all four arguments are provided
            % Use parallel processing to speed up the computation for each time series channel
            parfor j = 1:(s(2))
                % Construct the design matrix X by combining a column of ones (intercept), 
                % noise, detrended data, and high frequency noise
                X = ([ones((s(1)), 1), noise, detrenddata(:,j), highfreqnoise(:,j)]);
                y = (timeraw(:,j));  % Extract the corresponding raw time series data

                % Suppress warning for rank-deficient design matrix
                warning('off', 'stats:regress:RankDefDesignMat');
                b = regress(y, X);  % Perform linear regression to get coefficients
                WeightTimec = X * diag(b);  % Calculate the weighted time course
                clean(:,j) = y - (sum(WeightTimec, 2));  % Cleaned signal by removing noise components
            end
            
        case 3  % If three arguments are provided
            parfor j = 1:(s(2))
                % Construct the design matrix X with intercept, noise, and detrended data only
                X = ([ones((s(1)), 1), noise, detrenddata(:,j)]);
                y = (timeraw(:,j));

                warning('off', 'stats:regress:RankDefDesignMat');
                b = regress(y, X);
                WeightTimec = X * diag(b);
                clean(:,j) = y - (sum(WeightTimec, 2));
            end
            
        case 2  % If only two arguments (timeraw and noise) are provided
            parfor j = 1:(s(2))
                % Construct the design matrix X with intercept and noise only
                X = ([ones((s(1)), 1), noise]);
                y = (timeraw(:,j));

                warning('off', 'stats:regress:RankDefDesignMat');
                b = regress(y, X);
                WeightTimec = X * diag(b);
                clean(:,j) = y - (sum(WeightTimec, 2));
            end
            
    end

end
