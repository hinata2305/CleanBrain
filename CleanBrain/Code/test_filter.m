%Copyright © 2022 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

% This function processes time course data for a given set of noise variables and filter parameters. 
 
% Input: Requires two grey matter time course datasets (test and retest) 
%        and two noise time course datasets (test and retest), filter parameters,
%        and presentation output files to define the time course of interest.  

% Output: Applies classic SPM and Savitzky-Golay (SG) filters on data. 
%         Estimates statistics including auto-correlation, connectivity, true connectivity, 
%         reliability (ICC), and confidence intervals of ICC, among others.

clear all

% Specify the installation path for the PLOS_DATA_CODE folder.

CurrentScript=mfilename("fullpath");
[parentDir,~,~]=fileparts(CurrentScript);
addpath(parentDir, fullfile(parentDir,'..' ,'Logfiles'),fullfile(parentDir,'..','Timecourses'))

parpool('local')  % Initialize parallel computing pool

% Cognitive section of time course
part = [0 17];

% SG filter parameters
span = 311;   % Optimal filter detrend span
order = 40;   % Optimal filter detrend order
span2 = 3;    % Low pass filter optimal span
order2 = 1;   % Low pass filter optimal order
span3 = 69;   % Sub-optimal filter detrend span
order3 = 6;   % Sub-optimal filter detrend order
span4 = 15;   % Sub-optimal low pass filter span
order4 = 8;   % Sub-optimal low pass filter order

% SPM filter parameters
K{1}.RT = 1.24;      % Repetition time
K{1}.LChoice = 'none';   % No low pass filter applied
K{1}.HChoice = 'specify'; % High pass filter specified
K{1}.HParam = 128;  % Standard cutoff frequency

H{1}.RT = 1.24;      % Repetition time
H{1}.LChoice = 'hrf';  % HRF low pass filter
H{1}.HChoice = 'none';  % No high pass filter

G{1}.RT = 1.24;      % Repetition time
G{1}.LChoice = 'Gaussian'; % Gaussian low pass filter
G{1}.LParam = 2.48;  % Standard cutoff frequency
G{1}.HChoice = 'none';  % No high pass filter

% Load subject names from the provided file
load namelist_short67.mat

fol = length(namelist);  % Number of subjects

% Subject processing loop
for k = 1:fol
    % Load the data for the current subject
    subject = namelist{k};
    load([subject '.mat']);
    time1 = datafile.time1; % Test time series
    time2 = datafile.time2; % Retest time series
    
    % Extract reaction time data from presentation file (not used in analysis)
    [testa, testb, encode1a, encode2a, encode1b, encode2b, reactST, reactSR, reactMT, reactMR, respST, respSR, respMT, respMR] = extract_rt_rel_response_filter(subject);
    
    % Store reaction times in respective arrays
    stroopTest(:, k) = reactST';
    stroopRetest(:, k) = reactSR';
    memoryTest(:, k) = reactMT';
    memoryRetest(:, k) = reactMR';
    
    % Check if lengths of time courses for test and retest runs are correct (1 = correct)
    if (testa + testb) == 2
        % Define the actual start and endpoint for the time course analysis
        anfang1 = [encode1a(1) encode1b(1)];
        ende1 = [(encode1a(1) + 487) (encode1b(1) + 487)];
   
        % Obtain head motion data for Euclidean distance
        mov1 = datafile.mov1(:, [4 5 6 9]);
        mov2 = datafile.mov2(:, [4 5 6 9]);
        
        % Estimate head motion
        movepar1(k) = sum_move(mov1, mov2, anfang1, ende1);

        % First-level statistics for raw time courses
        raw_trans{k} = estimate_ICC_simplified3fft(time1, time2, anfang1, ende1);
        
        % trend of grey matter signal with optimal SG filter
        detrend_lib_1 = detrend_gol(time1, span, order, 1);
        detrend_lib_2 = detrend_gol(time2, span, order, 1);
        
        % trend grey matter signal with sub-optimal SG filter
        detrend_con_1 = detrend_gol(time1, span3, order3, 1);
        detrend_con_2 = detrend_gol(time2, span3, order3, 1);
        
        % Read noise signals obtained from freesurfer PCA (WM and ventricular and head motion components)
        noise1 = [datafile.noise_pca1(:, 1:10) datafile.movpca1(:, 1:2)];
        noise2 = [datafile.noise_pca2(:, 1:10) datafile.movpca2(:, 1:2)];
                    
        % Remove noise from time course without detrending
        denoise_1 = remove_noise(time1, noise1);
        denoise_2 = remove_noise(time2, noise2);
       
        % First-level statistics for the denoised time course
        denoise_trans{k} = estimate_ICC_simplified3fft(denoise_1, denoise_2, anfang1, ende1);

        % Denoise and detrend time course with optimal filter
        cleanLib_1 = remove_noise(time1, noise1, detrend_lib_1);
        cleanLib_2 = remove_noise(time2, noise2, detrend_lib_2);
        
        % First-level statistics for denoised and detrended time courses
        clean_lib_trans{k} = estimate_ICC_simplified3fft(cleanLib_1, cleanLib_2, anfang1, ende1);
   
        % Apply low pass filter using SG filter
        clean_Libfilter_1 = detrend_gol(cleanLib_1, span2, order2, 1);
        clean_Libfilter_2 = detrend_gol(cleanLib_2, span2, order2, 1);
        
        % First-level statistics for denoised, detrended, and low pass filtered time course
        clean_filter_lib_trans{k} = estimate_ICC_simplified3fft(clean_Libfilter_1, clean_Libfilter_2, anfang1, ende1);
        
        % Estimate pure noise from the difference of cleaned time courses
        noiseLib1 = cleanLib_1 - clean_Libfilter_1;
        noiseLib2 = cleanLib_2 - clean_Libfilter_2;
        
        % First-level statistics for pure noise time course
        noise_filter_lib_trans{k} = estimate_ICC_simplified3fft(noiseLib1, noiseLib2, anfang1, ende1);
           
        % Clean and detrend time course using suboptimal filter
        cleanCon_1 = remove_noise(time1, noise1, detrend_con_1);
        cleanCon_2 = remove_noise(time2, noise2, detrend_con_2);
        
        % First-level statistics for the cleaned, detrended time course
        clean_con_trans{k} = estimate_ICC_simplified3fft(cleanCon_1, cleanCon_2, anfang1, ende1);
        
        % Apply low pass filter using suboptimal filter
        clean_Confilter_1 = detrend_gol(cleanCon_1, span4, order4, 1);
        clean_Confilter_2 = detrend_gol(cleanCon_2, span4, order4, 1);
        
        % First-level statistics for the cleaned and low pass filtered time course
        clean_filter_con_trans{k} = estimate_ICC_simplified3fft(clean_Confilter_1, clean_Confilter_2, anfang1, ende1);
     
        % Estimate pure noise time courses using suboptimal filter
        noise1c = cleanCon_1 - clean_Confilter_1;
        noise2c = cleanCon_2 - clean_Confilter_2;
        
        % First-level statistics for pure noise time course with suboptimal filter
        noise_filter_con_trans{k} = estimate_ICC_simplified3fft(noise1c, noise2c, anfang1, ende1);
        
        % Apply classic low pass filter of SPM
        K{1}.row = 1:length(time1); 
        nK = spm_filter('set', K);
        spmfil = spm_filter('apply', nK, time1);  
        detrendHP1 = time1 - spmfil;    

        K{1}.row = 1:length(time2);    
        nK = spm_filter('set', K);
        spmfil = spm_filter('apply', nK, time2);  
        detrendHP2 = time2 - spmfil;    

        % Denoise and remove trend
        cleantempSPM1 = remove_noise(time1, noise1, detrendHP1);
        cleantempSPM2 = remove_noise(time2, noise2, detrendHP2);
        
        % First-level statistics for denoised SPM time course
        clean_SPM_trans{k} = estimate_ICC_simplified3fft(cleantempSPM1, cleantempSPM2, anfang1, ende1);
   
        % Apply HRF filter from SPM
        H{1}.row = 1:length(time1);  
        hK = spm_filter('set', H);
        clean_filtertempHRF1 = spm_filter('apply', hK, cleantempSPM1);
        
        H{1}.row = 1:length(time2); 
        hK = spm_filter('set', H);
        clean_filtertempHRF2 = spm_filter('apply', hK, cleantempSPM2);
       
        % First-level statistics for cleaned and HRF filtered time course
        clean_HRF_trans{k} = estimate_ICC_simplified3fft(clean_filtertempHRF1, clean_filtertempHRF2, anfang1, ende1);
        
        % Estimate noise from HRF filtered time courses
        noiseHRF1 = cleantempSPM1 - clean_filtertempHRF1;
        noiseHRF2 = cleantempSPM2 - clean_filtertempHRF2;
        
        % First-level statistics for noise time course after HRF filtering
        noise_HRF_trans{k} = estimate_ICC_simplified3fft(noiseHRF1, noiseHRF2, anfang1, ende1);
        
        % Apply Gaussian low pass filter
        G{1}.row = 1:length(time1); 
        gK = spm_filter('set', G);
        clean_filtertempGauss1 = spm_filter('apply', gK, cleantempSPM1);
   
        G{1}.row = 1:length(time2); 
        gK = spm_filter('set', G);
        clean_filtertempGauss2 = spm_filter('apply', gK, cleantempSPM2);
    
        % First-level statistics for Gaussian filtered time course
        clean_Gauss_trans{k} = estimate_ICC_simplified3fft(clean_filtertempGauss1, clean_filtertempGauss2, anfang1, ende1);
        
        % Estimate noise from Gaussian filtered time courses
        noiseGauss1 = cleantempSPM1 - clean_filtertempGauss1;
        noiseGauss2 = cleantempSPM2 - clean_filtertempGauss2;
        
        % First-level statistics for noise time course after Gaussian filtering
        noise_Gauss_trans{k} = estimate_ICC_simplified3fft(noiseGauss1, noiseGauss2, anfang1, ende1);
        
        % Apply a bad FFT filter
        clean_fft020_1 = filterFFT(cleanLib_1, 0, 0.20);
        clean_fft020_2 = filterFFT(cleanLib_2, 0, 0.20);
        
        % Estimate Fast Fourier transform for non-optimized filtering
        clean_fft_020_trans{k} = estimate_ICC_simplified3fft(clean_fft020_1, clean_fft020_2, anfang1, ende1);
        
        % Estimate noise from the clean time course using the FFT filter
        noise_grey20_p_1 = cleanLib_1 - clean_fft020_1;
        noise_grey20_p_2 = cleanLib_2 - clean_fft020_2;
        
        % First-level statistics for noise time course after FFT filtering
        noise_fft_020_trans{k} = estimate_ICC_simplified3fft(noise_grey20_p_1, noise_grey20_p_2, anfang1, ende1);
        
        % Extract sections of the cleaned time course with jitter
        bold1CleanLib{k} = extract_section_single_withJitter(cleanLib_1, encode1a, part, ende1(1));
        bold2CleanLib{k} = extract_section_single_withJitter(cleanLib_2, encode1b, part, ende1(2));
       
    else
        % Log subjects with incorrect length of time courses
        wronglength{k} = subject;
    end
end

subnum = length(clean_filter_con_trans);
expnum = length(raw_trans{1}.connectiona_trans);  % Number of experiments

% Combine statistics: raw data (not Fisher z back transformed)
raw_m = cell2mat(raw_trans);
% Statistics for Fisher z back transformation
[raw_bild, raw_dat] = z2r_simplified(raw_trans);

denoise_m = cell2mat(denoise_trans);
[denoise_bild, denoised_dat] = z2r_simplified(denoise_trans);

cleanLib_m = cell2mat(clean_lib_trans);
[clean_lib_bild, clean_lib_dat] = z2r_simplified(clean_lib_trans);

cleanCon_m = cell2mat(clean_con_trans);
[clean_con_bild, clean_con_dat] = z2r_simplified(clean_con_trans);

clean_filter_lib_m = cell2mat(clean_filter_lib_trans);
[clean_filter_lib_bild, clean_filter_lib_dat] = z2r_simplified(clean_filter_lib_trans);

clean_filter_con_m = cell2mat(clean_filter_con_trans);
[clean_filter_con_bild, clean_filter_con_dat] = z2r_simplified(clean_filter_con_trans);

cleanSPM_m = cell2mat(clean_SPM_trans);
[cleanSPM_bild, clean_SPM_dat] = z2r_simplified(clean_SPM_trans);

clean_filter_HRF_m = cell2mat(clean_HRF_trans);
[clean_filter_HRF_bild, clean_filter_HRF_dat] = z2r_simplified(clean_HRF_trans);

clean_filter_Gauss_m = cell2mat(clean_Gauss_trans);
[clean_filter_Gauss_bild, clean_filter_Gauss_dat] = z2r_simplified(clean_Gauss_trans);

% Prepare noise data from different filtering methods
noiseLib = cell2mat(noise_filter_lib_trans);
noiseCon = cell2mat(noise_filter_con_trans);
noiseGauss = cell2mat(noise_Gauss_trans);
noiseHRF = cell2mat(noise_HRF_trans);
noiseFFT = cell2mat(noise_fft_020_trans);
bold_pred1 = (cell2mat(bold1CleanLib));  % Cleaned time course for condition 1
bold_pred2 = (cell2mat(bold2CleanLib));  % Cleaned time course for condition 2

GrandBold = mean(bold_pred1, 2) + mean(bold_pred2, 2); % Average across conditions

% Estimating power spectral density using FFT
L = length(GrandBold);
NFFT = 2^nextpow2(L); % Next power of 2 from the length of GrandBold
temp = mean(2 * abs(fft(zscore(GrandBold), NFFT) / L), 2); % FFT computation
fftBold = temp(1:((NFFT / 2) + 1))';  % Retain only the positive frequencies

% Calculate mean FFT results from clean filters
fftlibfilt = (mean([clean_filter_lib_m.fftA]') + mean([clean_filter_lib_m.fftB]')) ./ 2;
fftconfilt = (mean([clean_filter_con_m.fftA]') + mean([clean_filter_con_m.fftB]')) ./ 2;

% Calculate mean FFT results for noise data
fftLibNoise = (mean([noiseLib.fftA]') + mean([noiseLib.fftB]')) ./ 2; 
fftConNoise = (mean([noiseCon.fftA]') + mean([noiseCon.fftB]')) ./ 2;
fftGNoise = (mean([noiseGauss.fftA]') + mean([noiseGauss.fftB]')) ./ 2;
fftHNoise = (mean([noiseHRF.fftA]') + mean([noiseHRF.fftB]')) ./ 2;
fftFNoise = (mean([noiseFFT.fftA]') + mean([noiseFFT.fftB]')) ./ 2;

% Remove first entry to avoid NaN results in the plot
fftBold(1) = nan;
fftLibNoise(1) = nan;
fftConNoise(1) = nan;
fftGNoise(1) = nan;
fftHNoise(1) = nan;
fftFNoise(1) = nan;

% Prepare frequency vector for plotting
num = (0:0.4/256:0.4);

% Plot connectivity and noise statistics
plot(num, [fftLibNoise; fftConNoise; fftGNoise; fftHNoise; fftFNoise; fftBold]');
plot(num, [fftLibNoise; fftConNoise; fftBold]');


% Compute absolute and Euclidean head motion estimates
moveabs = (([movepar1.abs1] + [movepar1.abs2])') ./ 2; 
mov = (([movepar1.cart1] + [movepar1.cart2])') ./ 2;

% Compute average autocorrelation for the complete time course (not reported in paper!)
auto1 = flipud([mean(tanh([raw_m.auto1_trans])), mean(tanh([denoise_m.auto1_trans])), mean(tanh([cleanLib_m.auto1_trans])), mean(tanh([cleanCon_m.auto1_trans])), mean(tanh([clean_filter_lib_m.auto1_trans])), mean(tanh([clean_filter_con_m.auto1_trans])), mean(tanh([cleanSPM_m.auto1_trans])), mean(tanh([clean_filter_HRF_m.auto1_trans])), mean(tanh([clean_filter_Gauss_m.auto1_trans]))]);
auto2 = flipud([mean(tanh([raw_m.auto2_trans])), mean(tanh([denoise_m.auto2_trans])), mean(tanh([cleanLib_m.auto2_trans])), mean(tanh([cleanCon_m.auto2_trans])), mean(tanh([clean_filter_lib_m.auto2_trans])), mean(tanh([clean_filter_con_m.auto2_trans])), mean(tanh([cleanSPM_m.auto2_trans])), mean(tanh([clean_filter_HRF_m.auto2_trans])), mean(tanh([clean_filter_Gauss_m.auto2_trans]))]);   

% Compile statistics for within-subject and between-path analysis
connectivity = [raw_dat.connection, denoised_dat.connection, clean_lib_dat.connection, clean_con_dat.connection, clean_filter_lib_dat.connection, clean_filter_con_dat.connection, clean_SPM_dat.connection, clean_filter_HRF_dat.connection, clean_filter_Gauss_dat.connection];
reliability = [(mean([raw_m.rel_trans], 2)), (mean([denoise_m.rel_trans], 2)), (mean([cleanLib_m.rel_trans], 2)), (mean([cleanCon_m.rel_trans], 2)), (mean([clean_filter_lib_m.rel_trans], 2)), (mean([clean_filter_con_m.rel_trans], 2)), (mean([cleanSPM_m.rel_trans], 2)), (mean([clean_filter_HRF_m.rel_trans], 2)), (mean([clean_filter_Gauss_m.rel_trans], 2))];

% Compile various reliability metrics
relmax = [raw_dat.relmax, denoised_dat.relmax, clean_lib_dat.relmax, clean_con_dat.relmax, clean_filter_lib_dat.relmax, clean_filter_con_dat.relmax, clean_SPM_dat.relmax, clean_filter_HRF_dat.relmax, clean_filter_Gauss_dat.relmax];
reldif = [raw_dat.reldif, denoised_dat.reldif, clean_lib_dat.reldif, clean_con_dat.reldif, clean_filter_lib_dat.reldif, clean_filter_con_dat.reldif, clean_SPM_dat.reldif, clean_filter_HRF_dat.reldif, clean_filter_Gauss_dat.reldif];
true_connectivity = [raw_dat.contrue, denoised_dat.contrue, clean_lib_dat.contrue, clean_con_dat.contrue, clean_filter_lib_dat.contrue, clean_filter_con_dat.contrue, clean_SPM_dat.contrue, clean_filter_HRF_dat.contrue, clean_filter_Gauss_dat.contrue];
rel_over = [raw_dat.relover, denoised_dat.relover, clean_lib_dat.relover, clean_con_dat.relover, clean_filter_lib_dat.relover, clean_filter_con_dat.relover, clean_SPM_dat.relover, clean_filter_HRF_dat.relover, clean_filter_Gauss_dat.relover];
corrupt = [(mean([raw_m.corrupt])), (mean([denoise_m.corrupt])), (mean([cleanLib_m.corrupt])), (mean([cleanCon_m.corrupt])), (mean([clean_filter_lib_m.corrupt])), (mean([clean_filter_con_m.corrupt])), (mean([cleanSPM_m.corrupt])), (mean([clean_filter_HRF_m.corrupt])), (mean([clean_filter_Gauss_m.corrupt]))];
NumberOfOverEst = [mean([raw_m.overestimated]), mean([denoise_m.overestimated]), mean([cleanLib_m.overestimated]), mean([cleanCon_m.overestimated]), mean([clean_filter_lib_m.overestimated]), mean([clean_filter_con_m.overestimated]), mean([cleanSPM_m.overestimated]), mean([clean_filter_HRF_m.overestimated]), mean([clean_filter_Gauss_m.overestimated])];
NumberOfOverEst_min = [mean([raw_m.sigover_min]), mean([denoise_m.sigover_min]), mean([cleanLib_m.sigover_min]), mean([cleanCon_m.sigover_min]), mean([clean_filter_lib_m.sigover_min]), mean([clean_filter_con_m.sigover_min]), mean([cleanSPM_m.sigover_min]), mean([clean_filter_HRF_m.sigover_min]), mean([clean_filter_Gauss_m.sigover_min])];

% Report averages for the within-subject, within-path analysis
crulmtdo = [tanh(mean(atanh(connectivity))); 
             tanh(mean(reliability)); 
             tanh(mean(atanh(relmax))); 
             tanh(mean(atanh(true_connectivity))); 
             tanh(mean(atanh(reldif))); 
             (tanh(mean(atanh(reldif))) ./ tanh(mean(atanh(true_connectivity)))) * -100; 
             (corrupt) ./ 561 * 100; 
             (NumberOfOverEst) ./ 561 * 100; 
             (NumberOfOverEst_min) ./ 561 * 100];

% Calculate number of ROIs detected with different reliability standards
tres = 0.4;  
fair = [mean((sum((tanh([raw_m.rel_trans])) > tres))); 
        mean((sum((tanh([denoise_m.rel_trans])) > tres)));  
        mean((sum((tanh([cleanLib_m.rel_trans])) > tres))); 
        mean((sum((tanh([cleanCon_m.rel_trans])) > tres)));
        mean((sum((tanh([clean_filter_lib_m.rel_trans])) > tres))); 
        mean((sum((tanh([clean_filter_con_m.rel_trans])) > tres))); 
        mean((sum((tanh([cleanSPM_m.rel_trans])) > tres))); 
        mean((sum((tanh([clean_filter_HRF_m.rel_trans])) > tres))); 
        mean((sum((tanh([clean_filter_Gauss_m.rel_trans])) > tres)))];

tres = 0.6;
good = [mean((sum((tanh([raw_m.rel_trans])) > tres))); 
        mean((sum((tanh([denoise_m.rel_trans])) > tres)));  
        mean((sum((tanh([cleanLib_m.rel_trans])) > tres))); 
        mean((sum((tanh([cleanCon_m.rel_trans])) > tres))); 
        mean((sum((tanh([clean_filter_lib_m.rel_trans])) > tres))); 
        mean((sum((tanh([clean_filter_con_m.rel_trans])) > tres))); 
        mean((sum((tanh([cleanSPM_m.rel_trans])) > tres))); 
        mean((sum((tanh([clean_filter_HRF_m.rel_trans])) > tres))); 
        mean((sum((tanh([clean_filter_Gauss_m.rel_trans])) > tres)))];

tres = 0.75;
excellent = [mean((sum((tanh([raw_m.rel_trans])) > tres)));
          mean((sum((tanh([denoise_m.rel_trans])) > tres)));
          mean((sum((tanh([cleanLib_m.rel_trans])) > tres)));  
          mean((sum((tanh([cleanCon_m.rel_trans])) > tres))); 
          mean((sum((tanh([clean_filter_lib_m.rel_trans])) > tres))); 
          mean((sum((tanh([clean_filter_con_m.rel_trans])) > tres))); 
          mean((sum((tanh([cleanSPM_m.rel_trans])) > tres))); 
          mean((sum((tanh([clean_filter_HRF_m.rel_trans])) > tres))); 
          mean((sum((tanh([clean_filter_Gauss_m.rel_trans])) > tres)))];

% Compile the fair, good, and excellent detection table
siz = (([fair; good; excellent]) ./ expnum) * 100;

% Probability analysis for within-subject, within-path
allrel = [mean(tanh([raw_m.rel_trans])); 
          mean(tanh([denoise_m.rel_trans])); 
          mean(tanh([cleanLib_m.rel_trans])); 
          mean(tanh([cleanCon_m.rel_trans])); 
          mean(tanh([clean_filter_lib_m.rel_trans])); 
          mean(tanh([clean_filter_con_m.rel_trans])); 
          mean(tanh([cleanSPM_m.rel_trans])); 
          mean(tanh([clean_filter_HRF_m.rel_trans])); 
          mean(tanh([clean_filter_Gauss_m.rel_trans]))];

allrelt = [mean(tanh([raw_m.rel_trans]')); 
           mean(tanh([denoise_m.rel_trans]')); 
           mean(tanh([cleanLib_m.rel_trans]')); 
           mean(tanh([cleanCon_m.rel_trans]')); 
           mean(tanh([clean_filter_lib_m.rel_trans]')); 
           mean(tanh([clean_filter_con_m.rel_trans]')); 
           mean(tanh([cleanSPM_m.rel_trans]')); 
           mean(tanh([clean_filter_HRF_m.rel_trans]')); 
           mean(tanh([clean_filter_Gauss_m.rel_trans]'))];


% Count the number of paths within subject with fair, good, and excellent reliability
fge_nopath_wthinsub = reshape((([fair; good; excellent]) ./ expnum) * 100,9,3)';
% Count the number of subjects showing fair, good, and excellent reliability in average connectivity
numsub = (reshape((sum([(allrel') > 0.4, (allrel') > 0.6, (allrel') > 0.75])) / subnum * 100, 9, 3))';
% Count the number of paths showing fair, good, and excellentreliability for average path
num_average_path = (reshape((sum([(allrelt') > 0.4, (allrelt') > 0.6, (allrelt') > 0.75])) / expnum * 100, 9, 3))';

% Compile total, fair, good, and excellent tables
fge = [fge_nopath_wthinsub; num_average_path; numsub];

% Create histograms for within-subject, within-ROI statistics
h(1) = figure;
hist(connectivity)
xlabel({'Connectivity Coefficient'}, 'FontSize', 20);
ylabel({'Number of Paths'}, 'FontSize', 20);

h(2) = figure;
hist(reliability)
xlabel({'Reliability Coefficient'}, 'FontSize', 20);
ylabel({'Number of Regions'}, 'FontSize', 20);

h(3) = figure;
hist(true_connectivity)
xlabel({'True Connectivity Coefficient'}, 'FontSize', 20);
ylabel({'Number of Paths'}, 'FontSize', 20);

h(4) = figure;
hist((atanh(reldif)) ./ (atanh(true_connectivity)) * -100)
xlabel({'Relative Overestimation Percentage'}, 'FontSize', 20);
ylabel({'Number of Paths'}, 'FontSize', 20);

window = 6;

% Estimate the effects of head motion on reliability using a sliding window approach
% The function calculates the mean reliability for subsamples with minimal head motion
[windowmean_r, move_r] = move_window_simplified(tanh(mean(atanh((raw_dat.rel)))'), mov, window);
[windowmean_hc] = move_window_simplified(tanh(mean(atanh((denoised_dat.rel)))'), mov, window);
[windowmean_c] = move_window_simplified(tanh(mean(atanh((clean_lib_dat.rel)))'), mov, window);
[windowmean_cfc] = move_window_simplified(tanh(mean(atanh((clean_filter_lib_dat.rel)))'), mov, window);

% Plot mean reliability against head motion
nn(1) = figure;
plot(move_r, [windowmean_r; windowmean_hc; windowmean_c; windowmean_cfc]')
xlabel({'Head Motion'}, 'FontSize', 18);
ylabel({'Within Subject Reliability Mean'}, 'FontSize', 18);

delete(gcp) % Close the parallel computing pool
