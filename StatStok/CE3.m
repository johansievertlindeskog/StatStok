%% Stationary stochastic processes, computer exercise 3
close all
clear all
clc
addpath('all_files')

% 2.1 AR(1)-process
C = 1;
A1 = [1 0.5];
A2 = [1 -0.5];
% Spectrum
[H1,w] = freqz(C,A1);
R1 = abs(H1).^2;
H2 = freqz(C,A2);
R2 = abs(H2).^2;
% Covariance function
H1 = freqz(C,A1,512,'whole');
Rd1 = abs(H1).^2;
r1 = ifft(Rd1);
H2 = freqz(C,A2,512,'whole');
Rd2 = abs(H2).^2;
r2 = ifft(Rd2);
% Realisation
x1 = filter(C,A1,normrnd(0,1,1,400));
x2 = filter(C,A2,normrnd(0,1,1,400));
% Figure, a1 = 0.5
figure(1)
subplot(3,1,1)
plot(w/2/pi,R1)
title('2.1 Spectral density, a_1 = 0.5')
xlabel('Frequency')
subplot(3,1,2)
stem(0:10,r1(1:11))
title('2.1 Covariance function, a_1 = 0.5')
xlabel('Time lag')
subplot(3,1,3)
plot(x1)
title('2.1 Realisation, a_1 = 0.5')
xlabel('Sample')
% Figure, a1 = -0.5
figure(2)
subplot(3,1,1)
plot(w/2/pi,R2)
title('2.1 Spectral density, a_1 = -0.5')
xlabel('Frequency')
subplot(3,1,2)
stem(0:10,r2(1:11))
title('2.1 Covariance function, a_1 = -0.5')
xlabel('Time lag')
subplot(3,1,3)
plot(x2)
title('2.1 Realisation, a_1 = -0.5')
xlabel('Sample')

% 2.2 ARMA(2)-process
% Poles and zeros
P = roots([1 -1 0.5]);
Z = roots([1 1 0.5]);
% Spectrum
[H1,w] = freqz([1 1 0.5],[1 -1 0.5]);
R1 = abs(H1).^2;
H2 = freqz([1 -1 0.5],[1 1 0.5]);
R2 = abs(H2).^2;
figure(3)
subplot(2,2,1)
zplane(Z,P)
title('2.2 Pole-zero plot')
subplot(2,2,2)
plot(w/2/pi,R1)
axis([0 0.5 0 40])
title('2.2 Spectral density')
subplot(2,2,3)
zplane(P,Z)
subplot(2,2,4)
plot(w/2/pi,R2)
axis([0 0.5 0 40])
xlabel('Frequency')

% 3.1 AR(2) processes
% Students work with the graphical interface, armagui, provided in the
% course code package, they should realise the following:
% A pole at a frequency, increases the amount of that frequency - the
% spectrum will have a peak at that frequency
% The angle between the pole and the positive real axis in the pole-zero
% plane is the angular frequency, omega = 2pif
% The distance from the origin of the unit circle to the pole corresponds
% to how much of that frequency the process has - the closer to the unit
% circle, the larger amount and higher peak in spectrum
% Moving the poles close to the origin gives smaller values of a1 and a2 in
% the AR polynomial, thus giving almost Xt = et, thus a white noise process
% Poles have to be inside the unit circle (and NOT on the unit circle) for
% the process to be stable, if poles are on or outside the unit circle the
% amplitude of the realisation will increase with each period

% 3.2 MA(q) processes
% Students work with the graphical interface, armagui, provided in the
% course code package, they should realise the following:
% The covariance function of an MA(q) is zero for time lags larger than q
% A zero at a frequency, lowers the amount of that frequency - the spectrum
% will have a dip in that frequency
% Moving a zero close to the unit circle more effectively cancels out a
% frequency, a zero close to the origin of the unit circle barely lowers
% the amount of the corresponding frequency

% 4 Speech modelling













close all force
clear all
clc
% Audio file
[x,Fs] = audioread('fa.wav');       % Fs = sampling frequency.
%sound(x,Fs)
n = length(x);                 
t = (0:n-1)/Fs;
figure
plot(t,x)
title('Speech signal')              
xlabel('Time (s)')
dt = 20*10^(-3);                    % Set length of time-windows (20 ms)
dn = Fs*dt;                         % Sample length of each time-window
                                    % Fs = 8000 Hz = 8000 samples/sec 
                                    % => 160 samples/20*10^(-3) s, so dn=160.
N_sec = floor(n/dn);                % Number of sections with length 20ms
M = 20;                             % Set AR model order
ar_param = zeros(N_sec,M+1);        % a_0 = 1, a_1,...,a_M have to be estimated, 
                                    % so we have a total of M+1 parameters
                                    % in N_sec sections.
for i = 1:N_sec                           % For each time section
    x_temp = x((i-1)*dn+1:i*dn);          % Pick out the i:th  20 ms section
    % 1st iteration: x(1): x(160); 2nd iteration: x(161):x(320);...
    % Student code
    U = zeros(size(x_temp,1)-20,20);      % U is a (160-20)x20 matrix
    for k = 1:20
        U(:,k) = -x_temp(21-k:end-k);
    end
    ar_temp = (U'*U) \ U'*x_temp(21:end);
    Q = (x_temp(21:end) - U*ar_temp)' * (x_temp(21:end) - U*ar_temp);
    e_temp = Q / size(U,1);
    ar_param(i,:) = [ar_temp ; e_temp];     % and save in ar_param. Why do we not need to save the first one?
end
whos x ar_param




% 4.1 Estimation of AR(20) parameters
[x,Fs] = audioread('fa.wav');
sound(x,Fs)
pause(3)
n = length(x);                 
t = (0:n-1)/Fs;
dt = 20*10^(-3);                    % Set length of time-windows (20 ms)
dn = Fs*dt;                         % Sample length of each time-window
N_sec = floor(n/dn);                % Number of sections with length 20ms

M = 20;                             % Set AR model order
ar_param = zeros(N_sec,M+1);

for i = 1:N_sec                           % For each time frame:
    x_temp = x((i-1)*dn+1:i*dn);          % Pick out the i:th  20 ms section
    
    % STUDENT CODE
    U = zeros(size(x_temp,1)-M,M);
    for k = 1:M
        U(:,k) = -x_temp(M+1-k:end-k);
    end
    ar_temp = (U'*U) \ U'*x_temp(M+1:end);
    Q = (x_temp(M+1:end) - U*ar_temp)' * (x_temp(M+1:end) - U*ar_temp); 
    e_temp = Q / size(U,1);     % felet är x_temp(21:end) - U*ar_temp)
    % END STUDENT CODE          % och sen är estimated variance av felet
                                % e_temp = Q / size(U,1)
    % ALTERNATIVE STUDENT CODE
    [ar_temp, e_temp] = arcov(x_temp,M);
    % First parameter of ar_temp is a0 = 1, it should be omitted
    % Either the transpose of ar_temp is needed or we need to save
    % [ar_temp , e_temp]
    ar_temp = ar_temp(2:end)';
    
    ar_param(i,:) = [ar_temp ; e_temp];
    
end

disp('4.1 Memory saved')
whos x ar_param

% 4.2 Calculate spectrum for some 20 ms
i = floor(N_sec/2);                % Pick the middle section.
x_temp = x((i-1)*dn+1:i*dn);
x_rec_temp = filter(1,[1 ar_param(i,1:end-1)],sqrt(ar_param(i,end))*randn(dn,1)); 
N_fft = 1024; % vi väljer att förklara felet från regressionen med white noise
              % som har variansen från regressionsfelet. 

f = (0:N_fft-1)/N_fft;
Px = abs(fft(x_temp,N_fft)).^2/N_sec ;

w  = exp(2i*pi*f);
Pa = (ar_param(i,end)) ./abs( polyval([1 ar_param(i,1:end-1)],w).' ).^2;

figure(4)
semilogy(f,Px);
hold on
semilogy(f,Pa,'r');
hold off
legend('Speech','AR Reconstruction');
title('4.2 Spectrum for 20ms of speech (and reconstruction)')

% 4.3 Recreate the audio in each time frame using the AR parameters
x_rec = zeros(n,1);
for jj = 1:N_sec
   x_rec((jj-1)*dn+1:jj*dn) = filter(1,[1 ar_param(jj,1:end-1)],sqrt(ar_param(jj,end))*randn(dn,1)); 
end   % och kom ihåg att e_temp är precis det som vi plockar ut i ar_param(jj,end) 

figure(5)
subplot(211);
plot(t,x);
title('4.3 Original sound')
subplot(212);
plot(t,x_rec);
title('4.3 Reconstructed sound')
xlabel('Time (s)')

% If stable, play sound
sound(x_rec,Fs)