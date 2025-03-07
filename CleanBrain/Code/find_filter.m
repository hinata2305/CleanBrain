%Copyright © 2022 Koten and Schüppen All rights reserved
%Important Notice: This code is not intended for medical applications 
%and does not have legal approval for such use. We strongly recommend 
%using FDA-approved software for any medical purposes. 

% Input: This function requires two noise time-course sets (test and retest) 
% for white matter, ventricles, head motion, etc., as well as two grey matter 
% time-course sets (test and retest), and two presentation output files utilized 
% for the presentation of fMRI stimuli and the recording of behavioral events.

% Output: The code identifies an optimal Savitzky-Golay (SG) filter. 
% Two versions are available: 
% Version 1 = Finds the optimal low-pass filter given a specific detrend parameter. 
% Optimal detrending parameters must be obtained from Version 2 for input into Version 1. 
% Version 2 = Determines the optimal detrend filter for a specific time course.

% Optimal filter findings for 67 subjects: 
% - For a search space of 245-6: optimal filter is 69/6 
% - For a search space of 483-18: optimal filter is 157/18 
% - For a search space of 483-40: optimal filter is 311/40 

% Specify the installation path for the PLOS_DATA_CODE folder.
% for the plos paper we used exp. 2 for the age paper exp. 1

CurrentScript=mfilename("fullpath");
[parentDir,~,~]=fileparts(CurrentScript);
addpath(parentDir, fullfile(parentDir,'..' ,'Logfiles'),fullfile(parentDir,'..','Timecourses'))

parpool('local');

version = 1; 

% Define the optimal detrend parameters for filtering.
span = 69;
order = 6;

% Specify the time range (begin-end) of the event of interest in TRs (Time Repetitions).
part = [0 17];

% Define the search space for the SG filter; windowa and windowb are the minimum 
% and maximum window sizes, ordera and orderb are the polynomial orders for filtering.
windowa = 3;
windowb = 485;
ordera = 1;
orderb = 484;

% Two experiments can be run: exp1 = spatial working memory, exp2 = verbal working memory.
exp = '2';


% Load the names of the subjects from the specified file.
load namelist_short67.mat;

fol = length(namelist);

for j = 1:fol
    subject = namelist{j};
    load([subject '.mat']);
  
    % Extract reaction time data, time course beginnings and endings for analysis.
    [testa, testb, encode1a, encode2a, encode1b, encode2b, reactST, reactSR, reactMT, reactMR, respST, respSR, respMT, respMR] = extract_rt_rel_response_filter(subject);
        
    if (testa + testb) == 2
       % Pre-process data for Version 1 of the analysis.
       if version == 1
            % Define the beginning of trials for a total of 24 trials.
            % Identify trial beginnings for test (A) and retest (B) runs.
            sub{j}.expA = eval(['encode' exp 'a']);
            sub{j}.expB = eval(['encode' exp 'b']);
                
            % Read the raw time course data.
            time1 = datafile.time1;
            time2 = datafile.time2;
                
            % Apply the optimal detrending parameter to the data.
            detrend1 = detrend_gol(time1, span, order, 1);
            detrend2 = detrend_gol(time2, span, order, 1);
                
            % Read the noise signal obtained from FreeSurfer PCA (1:5 are WM components, 
            % 6:10 are ventricle components, and 2 head motion components).
            noise1 = [datafile.noise_pca1(:, 1:10) datafile.movpca1(:, 1:2)];
            noise2 = [datafile.noise_pca2(:, 1:10) datafile.movpca2(:, 1:2)];
                
            % Remove noise and trends from grey matter timecourse data.
            cleantemp1 = remove_noise(time1, noise1, detrend1);
            cleantemp2 = remove_noise(time2, noise2, detrend2);
                 
            sub{j}.clean1 = cleantemp1;
            sub{j}.clean2 = cleantemp2;  
             
       else
            % Pre-process data for Version 2 of the analysis.
            sub{j}.expA = eval(['encode' exp 'a']);
            sub{j}.expB = eval(['encode' exp 'b']);
                
            time1 = datafile.time1;
            time2 = datafile.time2;
                
            sub{j}.Time1 = time1;
            sub{j}.Time2 = time2;
                
            % Read noise signals from FreeSurfer PCA.
            noise1 = [datafile.noise_pca1(:, 1:10) datafile.movpca1(:, 1:2)];
            noise2 = [datafile.noise_pca2(:, 1:10) datafile.movpca2(:, 1:2)];
                
            sub{j}.Noise1 = noise1;
            sub{j}.Noise2 = noise2;
                
            % Remove noise from grey matter timecourse using the regression function.
            sub{j}.clean1 = remove_noise(time1, noise1);
            sub{j}.clean2 = remove_noise(time2, noise2);
                
       end
     
    else
        % Record any subjects with inconsistent data lengths.
        wronglength{j} = subject;
    end
    
    clear datafile;
end

NOs = length(sub);

% Determine the optimal low-pass filter for the detrended data set and estimate the serial correlations of the predictor function.
if version == 1
    for k = 1:NOs
        tic
        [opti1, auto1, autosection1] = opti_gol_autonorm(sub{k}.clean1, sub{k}.clean2, sub{k}.expA, sub{k}.expB, part);
        opticor1(k, :,:,:,:) = opti1;
        autocor1(k, :,:,:,:) = auto1;
        autosec1(k, :,:,:,:) = autosection1;

        [opti2, auto2, autosection2] = opti_gol_autonorm(sub{k}.clean2, sub{k}.clean1, sub{k}.expB, sub{k}.expA, part);
        opticor2(k, :,:,:,:) = opti2;
        autocor2(k, :,:,:,:) = auto2;
        autosec2(k, :,:,:,:) = autosection2;
        
        toc
    end
               
    % Save the optimal low-pass filter results to a MAT file.
    save([parentDir 'optimal_lowpass.mat'], 'opticor1', 'opticor2', 'autocor1', 'autocor2', 'autosec1', 'autosec2');
    delete(gcp);  
else
    % Determine the optimal detrend filter for given time courses.
    parfor k = 1:NOs
        opti1 = opti_gol_trend_inter_2d(sub{k}.Time1, sub{k}.Noise1, sub{k}.clean2, windowa, windowb, ordera, orderb, sub{k}.expA, sub{k}.expB, part);
        opticor1(k, :,:,:,:) = ([opti1]);

        opti2 = opti_gol_trend_inter_2d(sub{k}.Time2, sub{k}.Noise2, sub{k}.clean1, windowa, windowb, ordera, orderb, sub{k}.expB, sub{k}.expA, part);
        opticor2(k, :,:,:,:) = ([opti2]);
    end
        
    % Save the optimal detrend filter results to a MAT file.
    save([parentDir 'optimal_detrend.mat'], 'opticor1', 'opticor2');
    delete(gcp);
end
