function [ Final_Signal,fs_O ] = De_convolve_twoSound( convovledSound,... 
    impulseSound)
% De-convolve_twoSound: - De-Convolves two sounds.

[h,fs_O] = audioread(convovledSound); % read the convolved sound
max_V = max(abs(h));

[d,fs_I] = audioread(impulseSound); % read the impulse sound
max_V = max(max(abs(h)), max_V);

% % Add an artificial non-zero number at the beggining - 
% d = vertcat(, d);

% Get the Original sound back by using de-convolution
[Final_Signal,r] = deconv(h,d);

% Normalise the convolved sound - 
Final_Signal = max_V*Final_Signal/max(abs(Final_Signal));

% save the new convolved sound - 
fileName = [convovledSound(1:end-4) '_DeConvlv.wav'];
audiowrite(fileName, Final_Signal, fs_O);

disp(' ')
disp(['The file has been saved as - ' fileName])
disp(' ')



end

