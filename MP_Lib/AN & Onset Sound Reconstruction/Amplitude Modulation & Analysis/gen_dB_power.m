function [yp, ydb, y] = gen_dB_power(filename)
% This function generates the power from the desibels, which is generated 
% from the magnitidues of the given signal.

[yt, fs] = wavread(filename);

y = yt(:,1);

ydb=mag2db(y);

yp=db2pow(ydb);

Tmax = length(yt)/fs;
t = (0+(1/fs)):(1/fs):Tmax;               % time vector

figure;
subplot(2,1,1);
plot(t,ydb);
hold on;
xlabel('Time (sec)');
ylabel('Desibel');
title('Signal Desibel');
grid on;
subplot(2,1,2);
plot(t,yp);
hold on;
xlabel('Time (sec)');
ylabel('Power');
title('Signal Power');
grid on;


end

