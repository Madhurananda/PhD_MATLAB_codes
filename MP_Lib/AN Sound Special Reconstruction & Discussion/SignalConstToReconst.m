function SignalConstToReconst(Sig_Name)
% 
newSigName = [Sig_Name(1:end-4) '_Left.wav'];
Stereo_To_Mono(Sig_Name, newSigName, 'Left');

AN_Construct_MP([Sig_Name(1:end-4) '_Left.wav'], 32, 100)
AN_Reconstruct([Sig_Name(1:end-4) '_Left_32_100_AN.mat'],1,1)
AN_Reconstruct([Sig_Name(1:end-4) '_Left_32_100_AN.mat'],1,0)
AN_Reconstruct([Sig_Name(1:end-4) '_Left_32_100_AN.mat'],0,1)
AN_Reconstruct([Sig_Name(1:end-4) '_Left_32_100_AN.mat'],0,0)

end