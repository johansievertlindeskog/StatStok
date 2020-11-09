%% Stationary stochastic processes, computer exercise 2
% close all
% clear all
% clc

% 2. Periodogram
load unknowndata.mat
m = mean(data);
disp('2. Periodogram')
disp(['mean: ' num2str(m)])
% Subtract mean
x = data-mean(data);
% Peiodogram
X = fft(x);
n = length(x);
Rhat = (X .* conj(X)) / n;
f = (0:n-1)/n;
% Plot with freq. 0 - 1
figure(1)
subplot(3,1,1)
plot(f,Rhat)
title('2. Periodogram')
% Plot with freq. -0.5 - 0.5
subplot(3,1,2)
plot(f-0.5,fftshift(Rhat))
title('2. Periodogram, fftshift')
% Periodogram with zero-padding
NFFT = 512;
X = fft(x,NFFT);
Rhat = (X .* conj(X)) / n;
f = (0:NFFT-1)/NFFT;
subplot(3,1,3)
plot(f,Rhat)
xlabel('Freq.')
title('2. Periodogram, zero-padding')
% Covariance function
rhat = ifft(Rhat);
figure(2)
subplot(2,1,1)
stem(0:20,rhat(1:21),'fill')
hold on
plot(0:16,(1 - (0:16)/16)*0.53,'r')
hold off
title('2. Estimated covariance function')
R = zeros(NFFT,1);
R(104) = 135; R(end-102) = 135;
r = ifft(R);
subplot(2,1,2)
stem(0:20,r(1:21),'fill')
title('2. True covariance function')
xlabel('Lags')
% NOTE: The true frequency is f_0 = 0.2, compare estimates with and without
% zero-padding
% NOTE: The data is only 16 samples long, the estimated covariance function
% is therefore not good, remember no lags > 15 can be measured and
% the factor 1 - |tau|/n


% 3. Investigation of different spectrum estimation techniques
% Filtered noise process
e = randn(500,1);
model.A = [1 -2.39 3.35 -2.34 0.96];
model.C = [1 0 1];
x = filter(model.C,model.A,e);
figure(3)
subplot(2,1,1)
plot(x)
title('3. Filtered noise realisation')
xlabel('Time')
% Spectral density
[H,w] = freqz(model.C,model.A,2048);
R = abs(H).^2;
subplot(2,1,2)
plot(w/2/pi,10*log10(R))
title('3. True spectral density')
xlabel('Freq.')
ylabel('10log_{10}(R)')

% 3. Periodogram vs modified periodogram
% Periodogram
NFFT = 4096;
figure(4)
subplot(2,1,1)
periodogram(x,[],NFFT);
xlabel('')
axis([0 1 -50 50])
title('3. Periodogram')
% Modified periodogram using a Hanning window
subplot(2,1,2)
periodogram(x,hanning(500),NFFT);
axis([0 1 -50 50])
title('3. Modified periodogram, Hanning window')
% NOTE: Compare width of peaks and if the zero (dip) is visible
% NOTE: The normal periodogram can be seen as having a rectangle window,
% while the Hanning window is tapered (goes to 0 at the edges)

% 3. Linear vs logarithmic scale
[Rhat,f] = periodogram(x,[],NFFT,1);
figure(5)
subplot(2,1,1)
plot(f,Rhat)
title('3. Periodogram, linear scale')
subplot(2,1,2)
semilogy(f,Rhat)
title('3. Periodogram, logarithmic scale')
xlabel('Freq.')

% 3. Welch method (window length L = 90)
figure(6)
subplot(2,1,1)
periodogram(x,[],NFFT);
xlabel('')
axis([0 1 -50 50])
title('3. Periodogram')
subplot(2,1,2)
pwelch(x,hanning(90),[],NFFT);
axis([0 1 -50 50])
title('3. Welch method, Hanning window')
% NOTE: Compare width of main lobes, smoothness, visibility of the zero and
% variance between realisations

% 3. Variance reduction
Rhate = periodogram(e,[],NFFT);
Rhatew = pwelch(e,hanning(90),[],NFFT);
disp('3. Variance reduction')
disp('Ratio of variance, periodogram and Welch method')
disp(var(Rhate)./var(Rhatew))
% NOTE: K = 10, so we expect around 10

% 4. Spectral analysis of EEG signal
% Light flicker at 12 Hz
load eegdata12.mat
fs = 1/data1.t(2);
L = floor(2*length(data1.x) / (3+1));
figure(7)
subplot(3,2,1)
periodogram(data1.x,[],NFFT,fs);
axis([0 30 -20 20])
title('4. data1, 12 Hz')
xlabel('')
subplot(3,2,2)
pwelch(data1.x,hanning(L),[],NFFT,fs);
axis([0 30 -20 20])
title('4. data1, 12 Hz')
xlabel('')
subplot(3,2,3)
periodogram(data2.x,[],NFFT,fs);
axis([0 30 -20 20])
title('4. data2, 12 Hz')
xlabel('')
subplot(3,2,4)
pwelch(data2.x,hanning(L),[],NFFT,fs);
axis([0 30 -20 20])
title('4. data2, 12 Hz')
xlabel('')
subplot(3,2,5)
periodogram(data3.x,[],NFFT,fs);
axis([0 30 -20 20])
title('4. data3, 12 Hz')
subplot(3,2,6)
pwelch(data3.x,hanning(L),[],NFFT,fs);
axis([0 30 -20 20])
title('4. data3, 12 Hz')
% Unknown light flicker
load eegdatax.mat
fs = 1/data1.t(2);
L = floor(2*length(data1.x) / (3+1));
figure(8)
subplot(3,2,1)
periodogram(data1.x,[],NFFT,fs);
axis([0 30 -20 20])
title('4. data1, x Hz')
xlabel('')
subplot(3,2,2)
pwelch(data1.x,hanning(L),[],NFFT,fs);
axis([0 30 -20 20])
title('4. data1, x Hz')
xlabel('')
subplot(3,2,3)
periodogram(data2.x,[],NFFT,fs);
axis([0 30 -20 20])
title('4. data2, x Hz')
xlabel('')
subplot(3,2,4)
pwelch(data2.x,hanning(L),[],NFFT,fs);
axis([0 30 -20 20])
title('4. data2, x Hz')
xlabel('')
subplot(3,2,5)
periodogram(data3.x,[],NFFT,fs);
axis([0 30 -20 20])
title('4. data3, x Hz')
subplot(3,2,6)
pwelch(data3.x,hanning(L),[],NFFT,fs);
axis([0 30 -20 20])
title('4. data3, x Hz')
% NOTE: Students should do this analysis in spekgui, but these are the
% relevant plots
% NOTE: Only 3 windows are used for the Welch method, otherwise we get too
% much smoothing of the spectral estimate


% 5. Identificaton of a stationary Gaussian process
load threeprocessdata.mat
% Covariance functions
Lag = 50;
[r1, lags] = xcov(y1,Lag+1);
[r2, lags] = xcov(y2,Lag+1);
[r3, lags] = xcov(y3,Lag+1);
figure(9)
subplot(3,1,1)
stem(lags,r1,'fill')
title('5.1 Covariance funtion, y1')
subplot(3,1,2)
stem(lags,r2,'fill')
title('5.1 Covariance funtion, y2')
subplot(3,1,3)
stem(lags,r3,'fill')
title('5.1 Covariance funtion, y3')
xlabel('Lags')
% Periodogram
[Rhat1,f] = periodogram(y1,[],NFFT,1);
Rhat2 = periodogram(y2,[],NFFT);
Rhat3 = periodogram(y3,[],NFFT);
% Welch method, K = 8
L = floor(2*length(y1) / (8+1));
Rhat1w = pwelch(y1,hanning(L),[],NFFT);
Rhat2w = pwelch(y2,hanning(L),[],NFFT);
Rhat3w = pwelch(y3,hanning(L),[],NFFT);
figure(10)
subplot(3,2,1)
semilogy(f,Rhat1)
title('5.1 Periodogram, y1')
subplot(3,2,2)
semilogy(f,Rhat1w)
title('5.1 Welch method, y1')
subplot(3,2,3)
semilogy(f,Rhat2)
title('5.1 Periodogram, y2')
subplot(3,2,4)
semilogy(f,Rhat2w)
title('5.1 Welch method, y2')
subplot(3,2,5)
semilogy(f,Rhat3)
title('5.1 Periodogram, y3')
xlabel('Freq.')
subplot(3,2,6)
semilogy(f,Rhat3w)
title('5.1 Welch method, y3')
xlabel('Freq.')
% NOTE: It can be benifitial to lower the number of windows for the Welch
% method from 10 slightly
% NOTE: The frequencies of all processes are very similar, therefore the
% covariance functions are very similar

% 5.2 Coherence spectrum
[Cx1y1,f] = mscohere(x1,y1,hanning(90),[],NFFT,1);
Cx3y1 = mscohere(x3,y1,hanning(90),[],NFFT);
figure(11)
subplot(1,2,1)
plot(f,Cx1y1)
axis([0 0.5 0 1])
title('5.2 Coherence spectrum, x1 and y1')
xlabel('Freq.')
subplot(1,2,2)
plot(f,Cx3y1)
axis([0 0.5 0 1])
title('5.2 Coherence spectrum, x3 and y1')
xlabel('Freq.')
% NOTE: Remember that the coherence spectrum = 1 for linear filters with no
% added noise
% NOTE: We are looking at finite realisations of random processes,
% therefore the results will not be exactly like the theory