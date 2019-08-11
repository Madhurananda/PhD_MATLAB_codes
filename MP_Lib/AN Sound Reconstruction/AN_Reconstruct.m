function [OutPut] = AN_Reconstruct(an_data_files, Ramp, SpikeCap,Jitrsp_tech)
% This function requires the values of 'Ramp' and 'SpikeCap' as 1 or 0. 

% This Re construction work reproduces the signals from the spike.
%   It has only one output file.
%   The output file saves all the important data used to produced the
%   signal.
%   Example Call: - AN_Reconstruct('CONGA_TUMBA_16_50_AN.mat', 1, 0, 0);

% The data should be loaded from the file
big.data = load(an_data_files);

% Pick out some useful information
% How many channels?
num_channels    = big(1).data.AN.channels;
% How many sen levels?
num_sen_levels  = big(1).data.AN.iterations;
% What is the length of the original input signal data?
length_sig = big(1).data.AN.datalength;
% What is the threshhold level values for sensitivity levels?
thresh_levels = big(1).data.AN.thresh_levels;
% What is the delay vector for this sound signal?
delayVector = big(1).data.AN.delayVector;
% What is the Sampling Rate of the original input sound signal?
Fs = big(1).data.AN.fs ;
% % What is the bit depth of the original input sound signal?
% nbits = big(1).data.AN.nbits ;
% What is the peak of the original signal?
PeakofSignal = big(1).data.AN.PeakofSignal;
% What are the center frequencies used for each channel by gammatone 
% filterbank?
cochCFs = big(1).data.AN.cochCFs;
% % Which way the spikes are constructed?
% Spike_Gen_Option = big(1).data.AN.Spike_Gen_Option;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is for taking all the data from that big cell array and assigning 
% them according to the number of channels 
% switch Spike_Gen_Option
%     case 0
%         new_zc_cell = AN_Reconstruct_AssignData_eachChannel(an_data_files);
%     case 1
%         new_zc_cell = big(1).data.AN.assigned_Spikes;
% end
new_zc_cell = big(1).data.AN.Spike_Assign_Channels; % This is according to 
% the new technique to reduce latency
% new_zc_cell = AN_Reconstruct_AssignData_eachChannel(an_data_files); 
% This was the initial way to assign spikes for each channel 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Action should be taken for delay vectors
% To deactivate the delay compensation, 'time_value' has to be taken from
% new_zc_cell, not from 'delay_zc_cell'. To do that search for 'time_value'
% and toggle the comments with the underneath.
delay_zc_cell=cell(1,num_channels);
for ch = 1:num_channels
    time_value = new_zc_cell{1,ch};
    new_time_value = time_value(:,1);
    for i = 1:length(new_time_value)
        new_time_value(i) = new_time_value(i)-((1.0)*delayVector(ch));
    end
    time_value(:,1) = new_time_value;
    delay_zc_cell{1,ch} = time_value;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The maximum value of the whole data set should be computed here. 
for ch = 1:num_channels
    % -- This is where the delay vector effects can be altered. -- %
    
%     time_value = new_zc_cell{1,ch};
%     delay_zc_cell{1,ch} = new_zc_cell{1,ch};
    
    time_value = delay_zc_cell{1,ch};
    
    new_time_value = time_value(:,1);
    if ch == 1
        max_Time_Value = max(new_time_value);
        min_Time_Value = min(new_time_value);
    end
    % Another array should be made to count the highest data in the whole 
    % set of data.
    max_new_timevalue = vertcat(max_Time_Value,new_time_value);
    min_new_timevalue = vertcat(min_Time_Value,new_time_value);
    max_Time_Value = max(max_new_timevalue);
    min_Time_Value = min(min_new_timevalue);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the part, which actually generates the signals for each channel.

% Is total no. of spike necessary (only necessary for maximum spike rate)?
totalnoSpike = false;

if Ramp== 1
    if SpikeCap == 1
        % 'new_channel_frequency_array' holds the maximum and minimum
        % frequncy of the sine waves generated between two consecutive
        % occurrences of spikes. It helps to understand the effect of spike
        % cap at the frequencies of each bit of sine wave.
        [multiplied_zc_cell, TotalNoSpikes, TotalSortedNoSpikes, ...
            new_channel_frequency_array] = ... 
            AN_Reconstruct_GenerateSignal_eachChannel(an_data_files, ... 
            delay_zc_cell, 200, Jitrsp_tech);
        totalnoSpike = true;
        % The output signal should be saved in appropriate name.
        new_Name = [an_data_files(1:end-7) '_NEW' '_JTR' num2str(Jitrsp_tech)];
        new_Name = [new_Name '.wav'];
        output_filename = [an_data_files(1:end-7) '_OutPut.mat'];
    elseif SpikeCap == 0
        [multiplied_zc_cell, new_channel_frequency_array] = ...
        AN_Reconstruct_GenerateSignal_eachChannel_NoSpikeCap(an_data_files, ... 
        delay_zc_cell);
        totalnoSpike = false;
        % The output signal should be saved in appropriate name.
        new_Name = [an_data_files(1:end-7) '_NEW_Ramp_NoSpikeCap'];
        % new_Name = [new_Name num2str(SNRdB) 'dBSNR'];
        new_Name = [new_Name '.wav'];
        output_filename = [an_data_files(1:end-7) '_OutPut_Ramp_NoSpikeCap.mat'];
    else
        disp('"SpikeRate" must have value 1 or 0.')
    end
    % The output file should consist of important data
    OutPut.new_channel_frequency_array=new_channel_frequency_array;
elseif Ramp == 0
    if SpikeCap == 1
        [multiplied_zc_cell, TotalNoSpikes, TotalSortedNoSpikes] = ... 
            AN_Reconstruct_GenerateSignal_eachChannel_NoRamp(an_data_files, ...
            delay_zc_cell, 200);
        totalnoSpike = true;
        % The output signal should be saved in appropriate name.
        new_Name = [an_data_files(1:end-7) '_NEW_NoRamp_SpikeCap'];
        new_Name = [new_Name '.wav'];
        output_filename = [an_data_files(1:end-7) '_OutPut_NoRamp_SpikeCap.mat'];
    elseif SpikeCap == 0
        multiplied_zc_cell = ... 
            AN_Reconstruct_GenerateSignal_eachChannel_NoRamp_NoSpikeCap ... 
            (an_data_files, delay_zc_cell);
        totalnoSpike = false;
        % The output signal should be saved in appropriate name.
        new_Name = [an_data_files(1:end-7) '_NEW_NoRamp_NoSpikeCap'];
        new_Name = [new_Name '.wav'];
        output_filename = ...
        [an_data_files(1:end-7) '_OutPut_NoRamp_NoSpikeCap.mat'];
    else
        disp('"SpikeRate" must have value 1 or 0.')
    end
else
    disp('"Ramp" must have value 1 or 0.')
end

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
    for i = 1:length_sig
        
        final_signal(ch,i) = time_value(i);
        % disp(length(time_value))
        
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
Final_Signal = Final_Signal/max(abs(Final_Signal));
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

% The output file should consist of important data
OutPut.new_zc_cell=new_zc_cell;
OutPut.delay_zc_cell=delay_zc_cell;
OutPut.multiplied_zc_cell=multiplied_zc_cell;
OutPut.Final_Signal=Final_Signal;
OutPut.max_Time_Value=max_Time_Value;
OutPut.min_Time_Value=min_Time_Value;


if (totalnoSpike == true)
    OutPut.TotalNoSpikes=TotalNoSpikes;
    OutPut.TotalSortedNoSpikes=TotalSortedNoSpikes;
end

% output_filename = [an_data_files(1:end-7) '_OutPut.mat'];
save(output_filename, 'OutPut')
disp(['The data has been saved to ' output_filename]); 

% Plot the original with the reconstructed signal
% Plot the original signal with the Reconstructed signal. 
figure
t = 0+(1/Fs):(1/Fs):length_sig/Fs;
plot(t, Final_Signal, 'r')
hold on
data  = audioread([an_data_files(1:end-13) '.wav']);
plot(t,data, 'b')

end
