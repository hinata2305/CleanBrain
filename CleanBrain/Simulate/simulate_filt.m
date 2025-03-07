
% Specify the installation path for the PLOS_DATA_CODE folder.
% Windows convention: current = 'my path\PLOS_DATA_CODE\'
% Linux convention: current= 'my path/PLOS_DATA_CODE/'
current = '/data/backup/Graz/FWF/SG_FILTER/PLOS_DATA_CODE/';
addpath([current 'Code'], [current 'Simulate'])

% length of simulated timecourse
timeL=487;

% paramteres of Gaussian filter
G{1}.RT = 1.24;                    % TR time
G{1}.LChoice = 'Gaussian';         % SPM HRF low pass filter
G{1}.LParam = 2.48;                % Standard cut of frequency
G{1}.HChoice = 'none';   

G{1}.row = 1:(timeL); 

gK = spm_filter('set',G);

% parameters of SG filter
span4=15;
order4=8;

% number of simulations per iteration
samp=1000000;

% iteration range
sim=100;

% pre define matrix for speed
rel_Auto=nan(sim,samp);
rel_Auto_filt=nan(sim,samp);
obsAuto1=nan(sim,samp);
obsAuto2=nan(sim,samp);
obsAutoFilt1=nan(sim,samp);
obsAutoFilt2=nan(sim,samp);

obsAutoFiltGauss1=nan(sim,samp);
obsAutoFiltGauss2=nan(sim,samp);

rel_Auto_filtGauss=nan(sim,samp);


for k = 1:samp
    
    parfor g = 1:sim
        
        % target auto correlation
        bb=0.01*g;

        % simualte autocorrelated timecourse and observed autocorrelation
        % of simulated timecourse
        [timecour1Auto, obsAuto1(g,k)]=createAutoCor(0,bb,1,timeL);
        [timecour2Auto, obsAuto2(g,k)]=createAutoCor(0,bb,1,timeL);
        
        % reliability of autocorrelated timecourse
        rel_Auto(g,k)=corr(timecour1Auto,timecour2Auto);
        
        % filter timecourse
        timecour1AutoFilt=detrend_gol_np(timecour1Auto,span4,order4,1);
        timecour2AutoFilt=detrend_gol_np(timecour2Auto,span4,order4,1);
        
        % auto correlation of filtered timecourse
        obsAutoFilt1(g,k)=corr(timecour1AutoFilt(1:end-1),timecour1AutoFilt(2:end));
        obsAutoFilt2(g,k)=corr(timecour2AutoFilt(1:end-1),timecour2AutoFilt(2:end));
        
        % reliability of filtered time course
        rel_Auto_filt(g,k)=corr(timecour1AutoFilt,timecour2AutoFilt);
        
        % apply guassian filter
        Clean_SPM_Gauss_p_1 = (spm_filter('apply',gK, timecour1Auto));
        Clean_SPM_Gauss_p_2 = (spm_filter('apply',gK, timecour2Auto));
       
        % auto correlation of filtered timecourse
        obsAutoFiltGauss1(g,k)=corr(Clean_SPM_Gauss_p_1(1:end-1), Clean_SPM_Gauss_p_1(2:end));
        obsAutoFiltGuass2(g,k)=corr(Clean_SPM_Gauss_p_2(1:end-1), Clean_SPM_Gauss_p_2(2:end));
        
        % reliability of filtered time course
        rel_Auto_filtGauss(g,k)=corr(Clean_SPM_Gauss_p_1,Clean_SPM_Gauss_p_2);
      


    end
    
end

save autosim.mat obsAuto1 obsAuto2 obsAutoFilt1 obsAutoFilt2 obsAutoFiltGauss1 obsAutoFiltGauss2 rel_Auto rel_Auto_filt  rel_Auto_filtGauss

% clear

Ra=rel_Auto';
Raf=rel_Auto_filt';
Rag=rel_Auto_filtGauss';

% mean reliabilty
 MeanRa=tanh(mean(atanh(Ra)));
 MeanRaf=tanh(mean(atanh(Raf)));
 MeanRag=tanh(mean(atanh(Rag)));

% standard deviation of relaibilty distribution
 StdRa=tanh(std(atanh(Ra)));
 StdRaf=tanh(std(atanh(Raf)));
 StdRag=tanh(std(atanh(Rag)));

 % standard deviation of relaibilty difference
 stdF_diff=tanh((std(atanh(Raf)))-(std(atanh(Ra))));
 stdG_diff=tanh((std(atanh(Rag)))-(std(atanh(Ra))));

% mean autocor of the two distribution
Oa=tanh    (mean((atanh(obsAuto1')+atanh(obsAuto2'))./2)); 
Oaf=tanh   (mean((atanh(obsAutoFilt1')+atanh(obsAutoFilt2'))./2)); 
Oag=tanh   (mean((atanh(obsAutoFiltGauss1')+atanh(obsAutoFiltGuass2'))./2));

plot([ StdRa;  StdRaf;  StdRag]')
plot([ MeanRa;  MeanRaf;  MeanRag]')
plot([ Oa;  Oaf;  Oag]')

hist([Ra(:,51) Raf(:,51) Rag(:,51)],100)




