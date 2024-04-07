%-------------------------------------------------
% PBMMI Matlab Assignment 1 - Part 2
%
% Coding the basic Karplus-Strong algorithm
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

M = round(Fs * dur);         % Duration of simulation in samples
N = round(Fs / f0 - 0.5);    % Delay line length, truncated


%rng(0)                      % Uncomment to set the random number generator seed for reproducibility
v = 2 * rand(1, N + 1) - 1;  % Vector of white noise of length N + 1
y = zeros(1, M);             % Output of the KS algorithm, initialized with zeros


%Dynamics Filter
% -------------------------------------------------------------------------
x1 = 0;                      % Initialize state variable, representing the sample x[n - 1]
% Loop over n = 0:N
for n = 0:N
    x0 = (1 - R) * v(n+1) + R * x1; % Calculate x0 (x[n]) x1 = v1 + x0
    y(n+1) = x0;                    % Write x0 into the output vector
    x1 = x0;                        % Update the state variable
end

% Main Karplus-Strong algorithm
% -------------------------------------------------------------------------
for n = N + 1:M - 1
    % Implement the difference equation
    y(n+1) = rho * (y(n - N) + y(n - N + 1)) / 2;
end

% Play the sound
soundsc(y, Fs);

% Plot the output waveform and its spectrum
t = (0:M - 1) / Fs; % Time vector
f = linspace(0, Fs / 2, floor(M / 2) + 1); % Frequency vector
Y = abs(fft(y)) / max(abs(fft(y)));        % Calculate the normalized FFT of the output signal

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
xlim([0 Fs/2]);
vline = line([f0 f0], ylim, 'Color', 'red', 'LineStyle', '--');
legend([p, vline], 'Output signal', 'f0');
