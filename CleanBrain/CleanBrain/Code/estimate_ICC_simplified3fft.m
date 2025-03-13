%Copyright © 2022 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

% Function: This function creates first-level statistics for fMRI reliability and connectivity studies.

% Input: Requires two time courses and the beginning (anfang) and end (ende) indices of a time course.

% Output: This function estimates the following:
% - FFT (Fast Fourier Transform)
% - TSNR (Temporal Signal-to-Noise Ratio)
% - Autocorrelation of time courses
% - Connectivity metrics
% - True connectivity estimation
% - Reliability (ICC) and confidence intervals of ICC
% - Detectable connectivity given reliability
% - Absolute overestimation of connectivity paths
% - Number of corrupt paths and number of overestimated paths.

function [sub] = estimate_ICC_simplified3fft(time1, time2, anfang, ende)

if nargin == 4
    % Extract the time course data for the specified indices
    timesec1 = (time1(anfang(1):ende(1), :));
    timesec2 = (time2(anfang(2):ende(2), :));
else
    timesec1 = time1;
    timesec2 = time2; 
end

s = size(timesec1);
L = s(1);
C = s(2);

% Estimate the Fourier Transform of the time course (MATLAB TOOLBOX)
NFFT = 2^nextpow2(L);
temp = mean(2 * abs(fft(zscore(timesec1), NFFT) / L), 2);
fftA = temp(1:((NFFT / 2) + 1));

L = length(timesec2);
NFFT = 2^nextpow2(L);
temp = mean(2 * abs(fft(zscore(timesec2), NFFT) / L), 2);
fftB = temp(1:((NFFT / 2) + 1));

% Estimate TSNR from untransformed data (meaningful in raw data)
sub.TSNRa = mean(mean(timesec1) ./ std(timesec1));
sub.TSNRb = mean(mean(timesec2) ./ std(timesec2));

% Estimate standard deviation from untransformed data
sub.STDa = mean(std(timesec1));
sub.STDb = mean(std(timesec2));

% Z-transform data to make different sessions comparable
timesec1 = zscore(timesec1);
timesec2 = zscore(timesec2);

% Estimate within ROI reliability using the ICC of choice (See new_icc for options)
% In this case, we use simple Pearson correlation
rel = arrayfun(@(ff) corr(timesec1(:, ff), timesec2(:, ff)), 1:C)';

reltrans = atanh(rel);

% Estimate connectivity and apply Fisher Z-transform
connectiona = atanh(corr(timesec1));
connectionb = atanh(corr(timesec2));

% Initialize matrices for analysis
s = size(timesec1);
reldif = zeros(s(2));
relmax = zeros(s(2));
contrue = zeros(s(2));
sigover_min = nan(s(2));
corrupt = zeros(s(2));
tres = 0.05 ./ (((s(2)) * (s(2) - 1)) / 2);

% Estimate the autocorrelations of the test time courses
unit = (s(2) + 1);
grasp = 1:unit:(s(2) * s(2));
[r_ww, ~] = xcorr(timesec1, 4, 'coeff');
auto = r_ww(1:4, grasp);
sub.auto1_trans = mean(atanh(auto'), 2);

% Estimate the autocorrelations of the retest
[r_ww, ~] = xcorr(timesec2, 4, 'coeff');
auto = r_ww(1:4, grasp);
sub.auto2_trans = mean(atanh(auto'), 2);

% Loop over paths a and b to estimate the detectable upper bound of connectivity
counter = 0;

for k = 1:s(2)
    for j = (k + 1):s(2)
        counter = counter + 1;
        
        % Estimate the detectable upper bound of connectivity
        cona = connectiona(j, k);
        conb = connectionb(j, k);
        meancon = (cona + conb) / 2;
        signie = sign(meancon);
        both = abs(meancon);
        conab(counter, 1) = cona;
        conab(counter, 2) = conb;

        % Only estimate paths with positive reliability
        if (reltrans(j) > 0) && (reltrans(k) > 0)
            % Estimate upper bound
            relmaxtemp = sqrt(reltrans(j) * reltrans(k));
            relmax(j, k) = relmaxtemp;
            % Estimate absolute overestimation
            reldif(j, k) = relmaxtemp - both;
            
            TOI = [timesec1(:, k) timesec1(:, j) timesec2(:, k) timesec2(:, j)];
            p1 = r_test_paired_abs(TOI(:, 1), TOI(:, 2), TOI(:, 3), 1, 1);
            p2 = r_test_paired_abs(TOI(:, 2), TOI(:, 1), TOI(:, 4), 1, 1);
            p3 = r_test_paired_abs(TOI(:, 3), TOI(:, 4), TOI(:, 1), 1, 1);
            p4 = r_test_paired_abs(TOI(:, 4), TOI(:, 3), TOI(:, 2), 1, 1);
            
            sigover_min(j, k) = min([p1 p2 p3 p4]);
    
            % Replace true connectivity with upper bound if greater than upper bound
            if both > relmaxtemp
                contrue(j, k) = relmaxtemp * signie;
            else
                contrue(j, k) = meancon; % Connectivity remains if below upper bound
            end
        else
            reldif(j, k) = 0 - both;
            relmax(j, k) = 0;
            contrue(j, k) = 0;
            corrupt(j, k) = 1;
        end
    end
end

[whole_ICC, whole_LB, whole_UB, ~, ~, ~, ~] = new_ICC(conab, 'A-1', 0.05, 0);

cortemp = corrcoef(conab);
sub.rel_trans = reltrans;
sub.connectiona_trans = connectiona;
sub.connectionb_trans = connectionb;
sub.connection_trans = (connectiona + connectionb) / 2;
sub.relmax_trans = relmax;
sub.reldif_trans = reldif;
sub.contrue_trans = contrue;
sub.corrupt = sum(sum(corrupt));
sub.overestimated = sum(sum(reldif < 0));
sub.pearwholebrain_trans = cortemp(2, 1);
sub.iccwholebrain_trans = atanh(whole_ICC);
sub.iccwholebrainLB_trans = atanh(whole_LB);
sub.iccwholebrainUB_trans = atanh(whole_UB);
sub.relover = (reldif ./ contrue) * 100;
sub.sigover_min = sum(sum(sigover_min < tres));
sub.fftA = fftA;
sub.fftB = fftB;

clear relmaxtemp maxcon cona conb sigover

end
