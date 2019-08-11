function [ OutPut ] = Spectrogram_Compare( sig_name_mono,calling_mode,ch )
% This function compares both spectrogram and signal ploting for one single sound file.
%   Can be called as: Spectrogram_Compare('testfile.wav', 'sound_file', 5)


sound_mode = strcmp(calling_mode,'sound_file');
signal_mode = strcmp(calling_mode,'signal');
if sound_mode == 1 
    str_check = sig_name_mono(end-2:end);
    switch strcmp(str_check,'wav')
        case 1
            % Read the wav file into a vector. Remember it will be normalised so the
            % peak amplitudes are 1. Also read the sample rate and bit depth.
            [sig,fs,nbits] = wavread(sig_name_mono);
            length_sig     = length(sig);
            disp(' ')
            disp(['Reading wave file: ' sig_name_mono])
            disp(' ')

        case 0
            error('The input file is not a .WAV file: Please convert it to .WAV and try again')
    end
elseif signal_mode == 1
    sig=sig_name_mono;
    length_sig = length(sig);
else
    error('Error: incorrect assignment of calling mode, check your first input parameter (should be either "sound_file" or "signal"')
end



if sound_mode == 1 
    t = 0:(1/fs):((length_sig/fs)-(1/fs));
    [S,F,T,P] = spectrogram(sig,256,250,256,fs);
    subplot(2,1,1)
    surf(T,F,10*log10(P),'edgecolor','none'); axis tight; 
    view(0,90);
    xlabel('Time (Seconds)'); ylabel('Hz');
    ylim([0 13000])

    hold on

    subplot(2,1,2)
    plot(t,sig); axis tight; 
    xlabel('Time (Seconds)'); ylabel('Amplitude');
    grid on
elseif signal_mode == 1
    figure(ch)
    title('Channel' +num2str(ch))
    t = 0:(1/fs):((length_sig/fs)-(1/fs));
    [S,F,T,P] = spectrogram(sig,256,250,256,fs);
    subplot(2,1,1)
    surf(T,F,10*log10(P),'edgecolor','none'); axis tight; 
    view(0,90);
    xlabel('Time (Seconds)'); ylabel('Hz');

    hold on

    subplot(2,1,2)
    plot(t,sig); axis tight; 
    xlabel('Time (Seconds)'); ylabel('Amplitude');
    grid on
end

if sound_mode == 1 
    OutPut.sig=sig;
    OutPut.fs=fs;
    OutPut.nbits=nbits;
    OutPut.length_sig=length_sig;
    OutPut.S=S;
    OutPut.F=F;
    OutPut.T=T;
    OutPut.P=P;
end

end

