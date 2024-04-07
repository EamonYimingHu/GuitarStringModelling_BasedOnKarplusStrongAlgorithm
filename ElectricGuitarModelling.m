%--------------------------------------------------------------------------
% PBMMI Matlab Assignment 1 - Beyond Basics
%
% This code implements a basic electric guitar model.
% Distortion (soft clipping) and feedback are added so that the output of the
% string is used as input to the delay line.
%
% Yiming HU
%--------------------------------------------------------------------------


% Customize from 0-100 to regulate the degree of distortion
distortion_gain = 70; 

% Guitar strings frequencies in Hz (E2, A2, D3, G3, B3)
string_freqs = [82.41, 110, 146.83, 196, 246.94];

% Parameters/initial values
% -------------------------------------------------------------------------
% Simulate each string and sum the sounds
Fs = 44100;                 % Sampling frequency in Hz
dur = 1.4;                  % Duration in seconds
rho = 0.275;                % Loss parameter
R = 0.95;                   % Dynamics parameter
distortion_volume = 0.81;   % Post-gain volume adjustment
feedback_gain = 0.5;        % Feedback gain
feedback_amount = 0.24;     % Feedback amount
y_total = [];               % Initialize an empty array to concatenate sounds
y_chord = zeros(1, round(Fs*dur));
distortion_gain = distortion_gain/1000+3.9;
for f0 = string_freqs
    y_string = karplus_strong_guitar_string(f0, Fs, dur, rho, R, distortion_gain, distortion_volume, feedback_gain, feedback_amount);
    y_total = [y_total, y_string];
    y_chord = y_chord + y_string;
end

soundsc([y_total, y_chord], Fs); % Play the sound

function y = karplus_strong_guitar_string(f0, Fs, dur, rho, R, gain, volume, feedback_gain, feedback_amount)
    % Calculate lengths and delays
    M = round(Fs * dur);      % Duration of simulation in samples
    Nexact = Fs / f0 - 0.5;
    N = floor(Nexact);        % Delay line length, truncated
    P = Nexact - N;
    C = (1 - P) / (1 + P);

    % Initialize signal vectors
    rng(0);
    v = 2 * rand(1, N + 1) - 1; % Vector of white noise of length N + 1
    y = zeros(1, M);            % Output of the KS algorithm, initialized with zeros
    feedback_signal = zeros(1, M); % Initialize the feedback signal buffer

    % Dynamics filter to simulate the plucking of the string
    x1 = 0; % Initialize state variable, representing the sample x[n - 1]
    for n = 0:N
        x0 = (1 - R) * v(n + 1) + R * x1; % Calculate x0 (x[n]) using dynamics parameter R
        y(n + 1) = x0;                   % Write x0 into the output vector
        x1 = x0;                         % Update the state variable for the next iteration
    end

    % Main Karplus-Strong algorithm with distortion and feedback
    % -------------------------------------------------------------------------
    yp1 = 0; % Initialize previous output sample for the loop
    for n = N + 1:M - 1
        % Karplus-Strong difference equation
        yp0 = C * y(n - N + 1) + y(n - N) - C * yp1;
        
        % Apply distortion to the signal
        y_distorted = soft_clip(0.5 * rho * (yp0 + yp1), gain, volume);
        
        % Mix in the feedback signal
        y_with_feedback = y_distorted + feedback_amount * feedback_signal(max(n - N, 1));
        
        % Store the output in the signal buffer
        y(n + 1) = y_with_feedback;
        
        % Update the feedback signal
        feedback_signal(n + 1) = feedback_gain * y_with_feedback;
        
        % Update the previous output sample
        yp1 = yp0;
    end
end

% Soft clipping function for distortion
% -------------------------------------------------------------------------
function y = soft_clip(x, gain, volume)
    % Apply pre-gain
    x = gain * x;
    % Soft clipping
    if x > 1
        y = 2 / 3;
    elseif x < -1
        y = -2 / 3;
    else
        y = x - (2 / 3) * x^3;
    end
    % Apply volume adjustment
    y = volume * y;
end
