%--------------------------------------------------------------------------
% PBMMI Matlab Assignment 1 - Beyond Basics
%
% This code implements a basic acoustic guitar model. 
% The signal is convolved with a pluck which is then convolved with the IR of a guitar body.
% According to Jaffe's paper, automated the calculation of the rho parameter based on a user-selected T60 value
% Tried building two sets of my own convolution functions, but they don't run as efficiently 
% as matlab's self-contained functions
% Tried create chords and short musical sequences by mixing and normalising the sum of the output vectors
%
% Yiming HU
%--------------------------------------------------------------------------

% Parameters/initial values
% -------------------------------------------------------------------------
Fs = 44.1e3;               % Sampling rate in Hz
Ts = 1/Fs;                 % Sampling period in seconds
t60 = 3;                   % Duration of simulation in seconds
f0Array = [164.81 196.00 146.83 196.00 130.81 164.81 146.83 196.00]; % Frequencies of plucks
pluckInterval = 0.2;       % Time interval between the starts of subsequent plucks in seconds

R = 0.95;                  % Dynamics parameter
M = round(Fs * dur);       % Calculate the duration of simulation in samples
pluck = audioread('pluck.wav');     % Recorded pluck
GBIR = audioread('Guitar_IR.wav');  % Guitar body IR

% Adjust the total duration calculation for plucks
totalDur = dur + (length(f0Array) - 1) * pluckInterval; % Total duration including overlaps
Mtotal = round(Fs * totalDur); % Total samples needed
yTotal = zeros(1, Mtotal);     % Initialize the total output vector

% Loop through each frequency in f0Array to simulate plucks
for i = 1:length(f0Array)
    f0 = f0Array(i);
    rho = exp(-6.91/(t60*f0))/cos(pi*f0*Ts);               % Loss parameter
    Nexact = Fs / f0 - 0.5;
    N = floor(Nexact);
    P = Nexact - N;
    C = (1 - P) / (1 + P);
    v = 2 * rand(1, N + 1) - 1; % Generate white noise
    y = zeros(1, M);            % Initialize the KS algorithm output


    % Dynamics Filter
    % ---------------------------------------------------------------------
    x1 = 0;
    for n = 0:N
        x0 = (1 - R) * v(n + 1) + R * x1;
        y(n + 1) = x0;
        x1 = x0;
    end

    % Main Karplus-Strong algorithm
    % ---------------------------------------------------------------------
    yp1 = 0;
    for n = N + 1:M - 1
        yp0 = C*y(n-N+1) + y(n-N) - C*yp1;
        y(n + 1) = 0.5*rho * (yp0 + yp1);
        yp1 = yp0;
    end
    % Convolution with pluck and GuitarBody IR
    % ---------------------------------------------------------------------
    y = conv(y,pluck);
    y = conv(y,GBIR);
    y = y/max(abs(y));

    % Calculate the start index for the current pluck
    startIndex = round((i-1) * Fs * pluckInterval) + 1;
    endIndex = startIndex + M - 1;

    % Ensure the array bounds are not exceeded
    endIndex = min(endIndex, length(yTotal));

    % Mix into the total output vector
    yTotal(startIndex:endIndex) = yTotal(startIndex:endIndex) + y(1:(endIndex-startIndex+1));
end

% Normalize and apply effects
yTotal = yTotal / max(abs(yTotal)); % Normalize

% Play the sound
soundsc(yTotal, Fs);

% My 1st Convolution algorithm:
% Using Time domain convolution theorem: Time domain convolution is 
% equivalent to frequency domain multiplication.
% -------------------------------------------------------------------------
% function y = myconv(x,h)
%     x = x(:);
%     h = h(:);
%     xn = [x;zeros(length(h)-1,1)];
%     hn = [h;zeros(length(x)-1,1)];
%     y = ifft(fft(xn).*fft(hn));
% end



% My 2nd Convolution algorithm:
% Using the definition of convolution to implement convolution of sequence
% but the computational complexity is high and the computation time is 
% too long when the length of the sequence is large
% -------------------------------------------------------------------------
% function y = myconv2(x, h)
%     % x is the input signal
%     % h is the impulse response (filter)
%     % y will be the output signal, the result of the convolution
% 
%     xLen = length(x);
%     hLen = length(h);
%     yLen = xLen + hLen - 1; % Length of the output signal
%     y = zeros(1, yLen); % Initialize the output signal with zeros
% 
%     % Extend x and h with zeros to match the length of the output signal
%     xPad = [x, zeros(1, hLen - 1)];
%     hPad = [h', zeros(1, xLen - 1)];
% 
%     % Perform the convolution operation
%     for i = 1:yLen
%         for j = 1:i
%             if j <= xLen && (i-j+1) <= hLen
%                 y(i) = y(i) + xPad(j) * hPad(i-j+1);
%             end
%         end
%     end
% end
