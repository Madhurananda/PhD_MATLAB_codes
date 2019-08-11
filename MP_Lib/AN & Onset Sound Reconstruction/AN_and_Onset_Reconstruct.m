function [OutPut] = AN_and_Onset_Reconstruct(AN_data_files, ... 
    Onset_data_files, Original_Onset_data_files, range) 
% This Re construction work reconstructes the signal from both the AN and
% Onset spikes. 
% It also require the values for 'range'. The values have to be  - 
%   '1' --> Male (channel - 15 to 40)
%   '2' --> Female (channel - 25 to 40)
%   '3' --> Musical (channel - 25 to 45)
%   It produces only one output file.
%   The output file saves all the important data used to produced the
%   signal.
%   
%   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the part, which actually generates the signals for each channel.

% For the very first time, no ramp and no spike cap has been considered.
[multiplied_zc_cell, new_zc_cell, delay_zc_cell_AN, delay_zc_cell_Onset, ... 
    max_Time_Value_AN, min_Time_Value_AN, max_Time_Value_Onset, ... 
    min_Time_Value_Onset, num_channels, length_sig, PeakofSignal, Fs] = ... 
    AN_and_Onset_Reconstruct_GenerateSignal(AN_data_files, ... 
    Onset_data_files, Original_Onset_data_files, range);

totalnoSpike = false;
% The output signal should be saved in appropriate name.
new_Name = [AN_data_files(1:end-7) '_NEW_AN_and_Onset'];
new_Name = [new_Name '.wav'];
output_filename = [AN_data_files(1:end-7) '_OutPut_AN_and_Onset.mat'];



% A big structure is needed to hold all the produced signals for different 
% channels
final_signal = zeros(num_channels,length_sig);
h4 = waitbar(0,'Please wait untill the final signal are produced...');
for ch = 1:num_channels
    time_value = multiplied_zc_cell{1,ch};
%     disp('Channel')
%     disp(ch)
%     disp('Length of the time vaule taken from the cell data')
%     disp(length(time_value))
    if ch == 3
        disp(ch)
    end
    for i = 1:length_sig
%         disp(i)
%         disp(ch)
        final_signal(ch,i) = time_value(i);
    end
%     % This work should be done to save the audio file for each channel
%     temp_signal = final_signal(ch,:)/max(abs(final_signal(ch,:)));
%     temp_signal = temp_signal*0.9;
%     new_Name = [an_data_files(1:end-7) '_NEW'];
%     new_Name = [new_Name num2str(ch)];
%     new_Name = [new_Name '.wav'];
%     wavwrite(temp_signal,44100,16,new_Name);
    waitbar(ch/num_channels)
end
close(h4)

%%%%%% The signals should be added up to produce the one signal.
% 1. This is to sum up all the signals for all channels
Final_Signal = sum(final_signal,1);
Final_Signal = 0.9*Final_Signal/max(abs(Final_Signal));
% The reconstructed signals should be normalised by multiplying the peak of
% the original signal. 
Final_Signal = Final_Signal*PeakofSignal;
Final_Signal = Final_Signal';
audiowrite(new_Name, Final_Signal,Fs);
disp(' ')
disp(' ')
disp('Re-Construction complete.');
disp('');
disp(['The audio file has been saved to ' new_Name]);
% End of 1

% Plot the original signal with the Reconstructed signal. 
figure
t = 0+(1/Fs):(1/Fs):length_sig/Fs;
plot(t, Final_Signal)
hold on
data  = audioread([AN_data_files(1:end-13) '.wav']);
plot(t,data, 'r')

% The output file should consist of important data
OutPut.new_zc_cell=new_zc_cell;
OutPut.delay_zc_cell_AN=delay_zc_cell_AN;
OutPut.delay_zc_cell_Onset=delay_zc_cell_Onset;
OutPut.multiplied_zc_cell=multiplied_zc_cell;
OutPut.Final_Signal=Final_Signal;
OutPut.max_Time_Value_AN=max_Time_Value_AN;
OutPut.min_Time_Value_AN=min_Time_Value_AN;
OutPut.max_Time_Value_Onset=max_Time_Value_Onset;
OutPut.min_Time_Value_Onset=min_Time_Value_Onset;


if (totalnoSpike == true)
    OutPut.TotalNoSpikes=TotalNoSpikes;
    OutPut.TotalSortedNoSpikes=TotalSortedNoSpikes;
end

% output_filename = [an_data_files(1:end-7) '_OutPut.mat'];
save(output_filename, 'OutPut')
disp(['The data has been saved to ' output_filename]); 



end
