function [Y] = AU_to_WAV(a_sound)
% This function converts an AU sound file into a WAV file.
% You must input the file as - 'filename.AU' and the output file will be 'filename.WAV'

if strcmp(a_sound(end-2:end) ,'.AU')
    [signal,Fs,nbits] = auread(a_sound);
    % The signal has to be modified so that it is loud enough.
    signal = signal/max(abs(signal));
    sound_name = [a_sound(1:end-3) '.wav'];
    wavwrite(signal,Fs,nbits,sound_name);

    Y.signal=signal;
    Y.Fs=Fs;
    Y.nbits=nbits;
else
    disp('This is not a AU file. ')
end



end
