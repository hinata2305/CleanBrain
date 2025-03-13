function XFFTFilt = filterFFT(time, freqa, freqb)
% filterFFT - Applies a bandpass filter to the input time-series data using the FFT method.
%
% Syntax: XFFTFilt = filterFFT(time, freqa, freqb)
%
% Inputs:
%   time  - Vector containing the time-series data to be filtered.
%   freqa - Lower cutoff frequency for the bandpass filter. To implement a 
%           lowpass filter, set freqa = 0.
%   freqb - Upper cutoff frequency for the bandpass filter. To implement a 
%           highpass filter, set freqb = 0.
%
% Outputs:
%   XFFTFilt - Filtered time-series data as a real-valued vector.
%
% Notes:
%   The function uses Fast Fourier Transform (FFT) to perform filtering and 
%   is set up to handle lowpass and highpass filtering based on input parameters.
%   The function zeroes out the FFT components outside the specified frequency 
%   range and then uses the Inverse FFT (IFFT) to reconstruct the time-domain 
%   signal.

Fs = 1/1.240; % Sampling frequency based on the time resolution of the input data
s = size(time); 
L = s(1); % Length of the input time series
NFFT = 2^nextpow2(L); % Next power of 2 based on the length of the input signal for efficient FFT

% Compute the Fourier Transform of the zero-centered time series
Y = fft(zscore(time), NFFT);

% Define the frequency unit based on the sampling frequency and NFFT size
unit = Fs / NFFT;

% Assign default values for the cutoff frequencies if unspecified
if freqa == 0
   freqa = unit; % Set lower cutoff to the lowest frequency unit if lowpass
end
if freqb == 0
   freqb = NFFT * unit; % Set upper cutoff to the highest frequency if highpass
end

% Create a zero-padded array for the filtered FFT coefficients
YFilt = zeros(size(Y) + [1, 0]);
Y(NFFT + 1, :) = 0 + 0 * 1i; % Ensure the Nyquist frequency component is zero
YFilt = Y; % Initialize YFilt with the original FFT coefficients

% Zero out frequencies below the lower cutoff and above the upper cutoff
YFilt(1:ceil(freqa/unit), :) = 0; % Remove low frequencies
YFilt(ceil(freqb/unit):NFFT + 2 - round(freqb/unit), :) = 0; % Remove high frequencies
YFilt(NFFT + 2 - round(freqa/unit):NFFT + 1, :) = 0; % Remove frequencies above the upper cutoff within the Nyquist region

% Perform inverse FFT to obtain the filtered time series in the time domain
XFFTFilt = real(ifft(YFilt, NFFT));
% Ensure the output vector is truncated to the original signal length
XFFTFilt(L + 1:end, :) = []; % Remove any excess samples generated during filtering

end

