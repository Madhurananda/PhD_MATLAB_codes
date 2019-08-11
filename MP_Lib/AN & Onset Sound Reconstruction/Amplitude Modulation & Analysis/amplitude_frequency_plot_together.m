function [yt_f, f, d] = amplitude_frequency_plot_together(filename)
% This function plots the amplitude and frequency of a given signal.
%   The input has to the signal and the time vector for that signal.
%   The frequency has been plotted around where it becomes maximum FIRST
%   TIME. So, if the frequency is highest at index 300 and 4500, the x-axis
%   for plotting frequency has been chosen as 0- 300*2 = 0:600. 

[yt, fs] = audioread(filename);

[m, n] =  size(yt);
if n>1
    yt = yt(:,1);
end

Tmax = length(yt)/fs;
t = (0+(1/fs)):(1/fs):Tmax;               % time vector


% Calculate the frequencies here - 
N = 2^nextpow2(length(t));
f=fs*(0:N-1)/N;
yt_f = (2/N)*abs(fft(yt,N));

% Plot the DSB-LC Signal
figure;
subplot(2,1,1);
plot(t,yt);
hold on;
xlabel('Time (sec)');
ylabel('Amplitude');
title('Signal Amplitude Plot');
grid on;
subplot(2,1,2);
d = find(yt_f==max(yt_f));
plot(f(1:d(1)*2),yt_f(1:d(1)*2));
% plot(f(1:1000),yt_f(1:1000));
% plot(f,yt_f);
xlabel('Frequency (Hz)');
ylabel('| Amplitude |');
title('Spectral Analysis (Single Sided PSD)');
grid on;

end

