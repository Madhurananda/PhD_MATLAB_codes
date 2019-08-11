function [AN,bmSig,output_filename] = AN_Construct_MP(sig_name_mono, n_senLavels, n_channels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function constructs the spikes from the given sound signal.
% 'onset_Option': - '0' - AN spike construction
%                   '1' - AN Onset Spike Construction. Onsets are generated
%                   by considering the very first occurences of a
%                   continuous train of spikes. So, each onset represents
%                   the beginning of a spike train after certain
%                   discontinuity of spikes in signal. 
% 
% if 'onset_Option' = 0:--> The spikes are coded by both Channels and 
% Sensitivity Levels 
% Example call - AN_Construct_MP('vocals.wav', 16, 50, 0)

[AN,bmSig,output_filename] = AN_1_SigGen_Mono_EditMP('internal',sig_name_mono,1,1,0,0,1, n_senLavels, n_channels);

% switch onset_Option
%     case 0
%         [AN,bmSig,output_filename] = AN_1_SigGen_Mono_EditMP('internal',sig_name_mono,1,1,0,0,1, n_senLavels, n_channels);
%     case 1
%         disp('')
%         disp('------>')
%         disp('AN_Onset spike coding has been used.')
%         disp('>--------')
%         [AN,bmSig,output_filename] = AN_Onset_SpikeGen(sig_name_mono,1, n_senLavels, n_channels);
% end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%