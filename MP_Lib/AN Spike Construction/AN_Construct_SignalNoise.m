function [AN,bmSig,OutPut] = AN_Construct_SignalNoise(sig_opt, ns_opt, SNRdB, n_senLavels, n_channels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function constructs the spikes for a Sound Signal in a noisy environment.

% 'sig_opt'     : Signal Option - which sound signal to use.
% 'ns_opt'      : Noise Option  - which way the noises will be generated 
% 'SNRdB'       : Signal To Ratio - in Decibel
% 'n_senLavels' : Number of Sensitivity Levels
% 'n_channels'  : Number of Channels



if sig_opt == 1
    [Sig_sig,fs_sig,nbits_sig] = wavread('powerwords.wav');
elseif sig_opt == 2
    
else
    disp('You must choose 1 or 2 for Male or Female Speech.')
end
PeakofSignal = max(abs(Sig_sig));

if ns_opt == 1
    [Sig_ns,fs_ns,nbits_ns] = wavread('StaticNoise.wav');
    if length(Sig_sig)>= length(Sig_ns)
        div = (length(Sig_sig)/length(Sig_ns));
        div = floor(div);
        % div_ex = rem(length(sig_sig),length(sig_ns));
        div_diff = (length(Sig_sig) - (div*length(Sig_ns)));
        for i=1:(div-1)
            Sig_ns = vertcat(Sig_ns,Sig_ns);
        end
        Sig_nsDiff = Sig_ns(1:div_diff);
        Sig_ns = vertcat(Sig_ns,Sig_nsDiff);
        if length(Sig_sig)~= length(Sig_ns)
            disp('The length of Signal and Noise is not same.')
        end
    end
elseif ns_opt == 2
    [Sig_ns,fs_ns,nbits_ns] = wavread('sea-shipSmall_Left.wav');
    if length(Sig_sig)>= length(Sig_ns)
        div = (length(Sig_sig)/length(Sig_ns));
        div = floor(div);
        % div_ex = rem(length(sig_sig),length(sig_ns));
        div_diff = (length(Sig_sig) - (div*length(Sig_ns)));
        for i=1:(div-1)
            Sig_ns = vertcat(Sig_ns,Sig_ns);
        end
        Sig_nsDiff = Sig_ns(1:div_diff);
        Sig_ns = vertcat(Sig_ns,Sig_nsDiff);
        if length(Sig_sig)~= length(Sig_ns)
            disp('The length of Signal and Noise is not same.')
        end
    else
        Sig_ns = Sig_ns(1:length(Sig_sig));
    end
elseif ns_opt == 3
    [Sig_ns,fs_ns,nbits_ns] = wavread('school-canteen_Left.wav');
    if length(Sig_sig)>= length(Sig_ns)
        div = (length(Sig_sig)/length(Sig_ns));
        div = floor(div);
        % div_ex = rem(length(sig_sig),length(sig_ns));
        div_diff = (length(Sig_sig) - (div*length(Sig_ns)));
        for i=1:(div-1)
            Sig_ns = vertcat(Sig_ns,Sig_ns);
        end
        Sig_nsDiff = Sig_ns(1:div_diff);
        Sig_ns = vertcat(Sig_ns,Sig_nsDiff);
        if length(Sig_sig)~= length(Sig_ns)
            disp('The length of Signal and Noise is not same.')
        end
    else
        Sig_ns = Sig_ns(1:length(Sig_sig));
    end
else
    disp('You must choose 1 or 2 or 3 for Static Noise or Sea Noise or Human bubble Noise.')
end

if (fs_sig ~= fs_ns)
    disp('The sampling rate of Signal and Noise is not same.')
else
    fs = fs_sig;
end

if (nbits_sig ~= nbits_ns)
    disp('The bit depth of Signal and Noise is not same.')
else
    nbits = nbits_sig;
end

noise_normaliser = PeakofSignal/(max(abs(Sig_ns)));
Sig_ns = Sig_ns.*noise_normaliser;
A_ns = max(abs(Sig_ns));
A_sig = 10^(SNRdB/20)*A_ns;
Sig_sig = Sig_sig.*(A_sig/PeakofSignal);
Sig_mix = Sig_sig + Sig_ns;
% Normalise it - 
Sig_mix = Sig_mix.*(0.9/max(abs(Sig_mix)));

sig_name = ['MixSound' '_' num2str(sig_opt) '_' num2str(ns_opt) '_' num2str(SNRdB) 'dBSNR' '.wav'];
wavwrite(Sig_mix,fs,nbits,sig_name);

[AN,bmSig,output_filename] = AN_1_SigGen_Mono_EditMP('internal',sig_name,1,1,0,0,1, n_senLavels, n_channels);

disp('The Construction work done.')

% OutPut = AN_Reconstruct(output_filename, Ramp, SpikeCap)
% OutPut = AN_Reconstruct(output_filename, 1, 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%