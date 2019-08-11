function [ OutPut ] = Spectrogram_Compare_anyFour( sig_name_mono1, sig_name_mono2, sig_name_mono3, sig_name_mono4)
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

[sig,fs,nbits] = wavread(sig_name_mono1);
[sig_2] = wavread(sig_name_mono2);
[sig_3] = wavread(sig_name_mono3);
[sig_4] = wavread(sig_name_mono4);
length_sig = length(sig);
disp(' ')
disp(['Reading wave file: ' sig_name_mono1])
disp(' ')
disp(['Reading wave file: ' sig_name_mono2])
disp(' ')
disp(['Reading wave file: ' sig_name_mono3])
disp(' ')
disp(['Reading wave file: ' sig_name_mono4])
disp(' ')


t = 0:(1/fs):((length_sig/fs)-(1/fs));

[S1,F1,T1,P1] = spectrogram(sig,256,250,256,fs);
subplot(2,2,1)
surf(T1,F1,10*log10(P1),'edgecolor','none'); axis tight; 
title([sig_name_mono1 ' Spectrogram'])
view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');
hold on

% subplot(2,2,2)
% plot(t,sig); axis tight;
% title([sig_name_mono ' Amplitude'])
% xlabel('Time (Seconds)'); ylabel('Amplitude');
% 
% grid on
% hold on

[S2,F2,T2,P2] = spectrogram(sig_2,256,250,256,fs);
subplot(2,2,2)
surf(T2,F2,10*log10(P2),'edgecolor','none'); axis tight; 
title([sig_name_mono2 ' Spectrogram'])
view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');
hold on

% subplot(2,2,4)
% plot(t,sig_Reconst); axis tight;
% title([Re_sig_name_mono ' Amplitude'])
% xlabel('Time (Seconds)'); ylabel('Amplitude');
% grid on

[S3,F3,T3,P3] = spectrogram(sig_3,256,250,256,fs);
subplot(2,2,3)
surf(T3,F3,10*log10(P3),'edgecolor','none'); axis tight; 
title([sig_name_mono3 ' Spectrogram'])
view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');
hold on

[S4,F4,T4,P4] = spectrogram(sig_4,256,250,256,fs);
subplot(2,2,4)
surf(T4,F4,10*log10(P4),'edgecolor','none'); axis tight; 
title([sig_name_mono4 ' Spectrogram'])
view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');
% hold on

OutPut.sig=sig;
OutPut.sig_2=sig_2;
OutPut.sig_3=sig_3;
OutPut.sig_4=sig_4;
OutPut.fs=fs;
OutPut.nbits=nbits;
OutPut.length_sig=length_sig;
OutPut.S1=S1;
OutPut.F1=F1;
OutPut.T1=T1;
OutPut.P1=P1;
OutPut.S2=S2;
OutPut.F2=F2;
OutPut.T2=T2;
OutPut.P2=P2;
OutPut.S3=S3;
OutPut.F3=F3;
OutPut.T3=T3;
OutPut.P3=P3;
OutPut.S4=S4;
OutPut.F4=F4;
OutPut.T4=T4;
OutPut.P4=P4;


end

