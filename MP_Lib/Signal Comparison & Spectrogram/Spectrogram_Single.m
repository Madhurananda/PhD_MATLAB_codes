function [ OutPut ] = Spectrogram_Single( sig_name_mono )
% This function plots the spectrogram only for a single sound file.

[sig,fs] = audioread(sig_name_mono);
length_sig     = length(sig);
disp(' ')
disp(['Reading wave file: ' sig_name_mono])
disp(' ')

figure
[S,F,T,P] = spectrogram(sig,256,250,256,fs);
surf(T,F,10*log10(max(P, 10^-18)),'edgecolor','none'); 
% surf(T,F,10*log10(P),'edgecolor','none'); 
axis tight;  
view(0,90);
xaxisMax = (length_sig/fs)+ ((length_sig/fs)*1)/100;
axis([0 xaxisMax 0 8000]);
xlabel('Time (Seconds)'); ylabel('Frequency (Hz)');


% figure
% spectrogram(sig,100,80,100,fs,'yaxis')
% view(-77,72)
% shading interp
% colorbar off
% [s,f,t,p] = spectrogram(sig,100,80,100,fs);
% [q,nd] = max(10*log10(p));
% hold on
% plot3(t,f(nd),q,'r','linewidth',4)
% hold off
% colorbar
% view(2)


OutPut.sig=sig;
OutPut.fs=fs;
OutPut.length_sig=length_sig;
OutPut.S=S;
OutPut.F=F;
OutPut.T=T;
OutPut.P=P;

end

