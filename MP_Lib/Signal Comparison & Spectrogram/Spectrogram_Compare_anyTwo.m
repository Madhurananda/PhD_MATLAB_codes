function [ OutPut ] = Spectrogram_Compare_anyTwo( sig_name_mono, Re_sig_name_mono)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% sig_name_mono_Reconst = [sig_name_mono(1:end-4) '_NEW'];
% sig_name_mono_Reconst = [sig_name_mono_Reconst '.wav'];
% 
% str_check = sig_name_mono(end-2:end);
% switch strcmp(str_check,'wav')
%     case 1
%         % Read the wav file into a vector. Remember it will be normalised so the
%         % peak amplitudes are 1. Also read the sample rate and bit depth.
%         [sig,fs,nbits] = wavread(sig_name_mono);
%         [sig_Reconst] = wavread(sig_name_mono_Reconst);
%         length_sig = length(sig);
%         disp(' ')
%         disp(['Reading wave file: ' sig_name_mono])
%         disp(' ')
%         disp(['Reading wave file: ' sig_name_mono_Reconst])
%         disp(' ')
%         
%     case 0
%         error('The input file is not a .WAV file: Please convert it to .WAV and try again')
% end

[sig,fs,nbits] = wavread(sig_name_mono);
[sig_Reconst] = wavread(Re_sig_name_mono);
length_sig = length(sig);
disp(' ')
disp(['Reading wave file: ' sig_name_mono])
disp(' ')
disp(['Reading wave file: ' Re_sig_name_mono])
disp(' ')


t = 0:(1/fs):((length_sig/fs)-(1/fs));
[S,F,T,P] = spectrogram(sig,256,250,256,fs);
subplot(2,2,1)
surf(T,F,10*log10(P),'edgecolor','none'); axis tight; 
title([sig_name_mono ' Spectrogram'])
view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');

hold on

subplot(2,2,2)
plot(t,sig); axis tight;
title([sig_name_mono ' Amplitude'])
xlabel('Time (Seconds)'); ylabel('Amplitude');

grid on
hold on

[Sr,Fr,Tr,Pr] = spectrogram(sig_Reconst,256,250,256,fs);
subplot(2,2,3)
surf(Tr,Fr,10*log10(Pr),'edgecolor','none'); axis tight; 
title([Re_sig_name_mono ' Spectrogram'])
view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');


hold on

subplot(2,2,4)
plot(t,sig_Reconst); axis tight;
title([Re_sig_name_mono ' Amplitude'])
xlabel('Time (Seconds)'); ylabel('Amplitude');
grid on


OutPut.sig=sig;
OutPut.sig_Reconst=sig_Reconst;
OutPut.fs=fs;
OutPut.nbits=nbits;
OutPut.length_sig=length_sig;
OutPut.S=S;
OutPut.F=F;
OutPut.T=T;
OutPut.P=P;
OutPut.S_Reconst=Sr;
OutPut.F_Reconst=Fr;
OutPut.T_Reconst=Tr;
OutPut.P_Reconst=Pr;


end

