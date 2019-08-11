function [ ] = Add_Zeros_atStart( ori_sig, recon_sig)
% This function add zeros at the beginning of the given sound signal.
% This function is made to add zeros specially for the pure frequency sine
% waves. 

[sig_ori, Fs, nbits] = wavread(ori_sig);
% Zeros to Add - 
ZeroToAdd = zeros(10000,1);
% Concat zeros at the beginning and at the end.
new_sig_ori = vertcat(ZeroToAdd,sig_ori,ZeroToAdd);

wavwrite(new_sig_ori,Fs,nbits,recon_sig);

end