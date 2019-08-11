function [ Final_Signal,fs_O ] = convolve_twoSound( sound1, sound2 )
% convolve_twoSound: - Convolves two sounds.
%   'Sound1' should be the original sound.
%   'Sound2' should be the impulse sound, which has to be convolved with
%   the original sound 'Sound1'.

[h,fs_O] = audioread(sound1); % read the original sound
max_V = max(abs(h));

[d,fs_I] = audioread(sound2); % read the impulse sound
max_V = max(max(abs(h)), max_V);

if (fs_O == fs_I)
    y=(conv(d,h)); % Convolve the two

    % Normalise the convolved sound - 
    Final_Signal = max_V*y/max(abs(y));

    % save the new convolved sound - 
    fileName = [sound1(1:end-4) '_Convlv_' sound2(1:end-4) '.wav'];
    audiowrite(fileName, Final_Signal, fs_O);

    disp(' ')
    disp(['The file has been saved as - ' fileName])
    disp(' ')
else
    disp(' ')
    disp('Sorry, Covolving is not possible as the sampling rate of two input sounds are different.')
    disp(['Original Sound Sample Rate: ' num2str(fs_O)])
    disp(['Impulse Sound Sample Rate: ' num2str(fs_I)])
    disp(' ')
end



end

