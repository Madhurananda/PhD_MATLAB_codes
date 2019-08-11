function [ ] = Signal_Compare_Ori_Recon( ori_sig, recon_sig)
% This function compares any two signals by plotting their graphs

[sig_ori, Fs] = wavread(ori_sig);
sig_recon = wavread(recon_sig);

if (length(sig_ori)==length(sig_recon))
    timeAxis = 0:(1/Fs):((length(sig_ori)/Fs)-(1/Fs));
    plot(timeAxis,sig_ori,'r')
    hold on
    plot(timeAxis,sig_recon,'b')
    grid on
else
    disp('The lengths are not equal.')
end

end