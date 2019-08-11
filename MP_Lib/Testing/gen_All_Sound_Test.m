function gen_All_Sound_Test( sound, type )
% This is to reconstruct the sounds by all our techniques: - 
%   'AN' Sound Reconstuction
%   'AN' and 'Onset' Sound Reconstruction
%   Thomas's spike technique
%
%   A call of this function might be like - gen_All_Sound_Test('CONGA_TUMBA.wav', 3)
%   The input type -> 1.Male  2.Female  3.Musical

% We need to see that the sound is .wav or not. 

Stereo_To_Mono(sound,sound,'Left')

% Construct the AN spikes first - 
AN_Construct_MP(sound, 16, 50, 0);
% Then, regenerate the sound - 
AN_sound_file = [sound(1:end-4) '_16_50_AN.mat'];
AN_Reconstruct(AN_sound_file, 1, 0, 0);

% Next, generate the Onset spikes 
generateonsetspikes1_mono(AN_sound_file, 1, [500 1100 25], 0.0015, 1000, 1)
generateonsetspikes1_mono(AN_sound_file, 2, [500 1100 25], 0.0015, 1000, 1)
% Then, regenerate the sound from both AN and Onset spikes
ANOnset_sound_file = [AN_sound_file(1:end-4) '_ANOnset.mat'];
OriOnset_sound_file = [AN_sound_file(1:end-4) '_OriginalOnset.mat'];
AN_and_Onset_Reconstruct(AN_sound_file, ANOnset_sound_file, OriOnset_sound_file, type)


% Next, do the regeneration for Thomas's work
Thomas_code_Work(sound, 5);


disp(['RECONSTRUCTION FOR ' sound ' HAS BEEN COMPLETED.'])
clear all, close all, clc

end

