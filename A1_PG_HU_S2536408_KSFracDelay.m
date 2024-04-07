%-------------------------------------------------
% PBMMI Matlab Assignment 1 - Part 2
%
% Coding the tuning-corrected Karplus-Strong algorithm
%
% Yiming HU 
%-------------------------------------------------


% Parameters/initial values
% -------------------------------------------------------------------------
Fs = 44.1e3;                 % Sampling rate in Hz
dur = 2;                     % Duration of simulation in seconds
f0 = 882;                    % Fundamental frequency of the string in Hz
rho = 0.998;                 % Loss parameter
R = 0.95;                    % Dynamics parameter

M = round(Fs * dur);         % Calculate the duration of simulation in samples
Nexact = Fs / f0 - 0.5;      % Exact delay line length calculation
N = floor(Nexact);           % Delay line length, truncated to remove fractional part
P = Nexact - N;              % Fractional delay
C = (1 - P) / (1 + P);       % All-pass filter coefficient for fractional delay compensation

%rng(0)                      % Uncomment to set the random number generator seed for reproducibility
v = 2 * rand(1, N + 1) - 1;  % Generate a vector of white noise of length N + 1
y = zeros(1, M);             % Initialize the output of the KS algorithm with zeros

% Dynamics Filter
% -------------------------------------------------------------------------
x1 = 0;                      % Initialize state variable, representing the sample x[n - 1]
for n = 0:N
    x0 = (1 - R) * v(n + 1) + R * x1; % Calculate x0 (x[n]) using dynamics parameter
    y(n + 1) = x0;                    % Write x0 into the output vector
    x1 = x0;                          % Update the state variable for the next iteration
end

% Main Karplus-Strong algorithm
% -------------------------------------------------------------------------
yp1 = 0;                     % Initialize the previous output sample for the loop
for n = N + 1:M - 1
    yp0 = C*y(n-N+1) + y(n-N) - C*yp1; % Calculate the new sample using the all-pass filter
    y(n + 1) = 0.5*rho * (yp0 + yp1);  % Apply the loss factor and averaging
    yp1 = yp0;                         % Update the previous output sample for the next iteration
end

% Play the sound 
soundsc(y, Fs);

% Plot the output waveform and its spectrum
% -------------------------------------------------------------------------
t = (0:M - 1) / Fs;                        % Generate a time vector for the x-axis
f = linspace(0, Fs / 2, floor(M / 2) + 1); % Generate a frequency vector for the x-axis of the spectrum plot
Y = abs(fft(y))/max(abs(fft(y)));          % Calculate the normalized FFT of the output signal

% Waveform plot
subplot(2, 1, 1);
plot(t, y);
xlabel('Time (s)');
ylabel('Amplitude');
title('Output Waveform');
xlim([0 dur]);
legend('Output signal');

% Spectrum plot
subplot(2, 1, 2);
p = plot(f, abs(Y(1:floor(M / 2) + 1)));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum');
xlim([0 Fs / 2]);
legend('Output signal');
vline = line([f0 f0], ylim, 'Color', 'red', 'LineStyle', '--'); % Mark the fundamental frequency
legend([p, vline], 'Output signal', 'f0');
