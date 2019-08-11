%function [Final_Signal,multiplied_zc_cell,new_zc_cell,zc_cell] = AN_Reconstruct(an_data_files)
function [OutPut] = AN_Reconstruct_NoRamp_NoSpikeCap(an_data_files, channel)
%function [delay_zc_cell] = AN_Reconstruct(an_data_files, channel)

% This Re construction work reproduces the signals from the spike.
%   It has only one output file.
%   The output file saves all the important data used to produced the
%   signal.

% The data should be loaded from the file
big.data = load(an_data_files);

% Pick out some useful information
% How many channels?
num_channels    = big(1).data.AN.channels;
% How many sen levels?
num_sen_levels  = big(1).data.AN.iterations;
% What is the length of the original input signal data?
length_sig = big(1).data.AN.datalength;
% What is the threshhold level values?
thresh_levels = big(1).data.AN.thresh_levels;
% What is the delay vector for this sound signal?
delayVector = big(1).data.AN.delayVector;
% What is the frequency of the original input sound signal?
Fs = big(1).data.AN.fs ;
% What is the peak of the original signal?
PeakofSignal = big(1).data.AN.PeakofSignal;
% What are the center frequencies used for each channel by gammatone filterbank?
cochCFs = big(1).data.AN.cochCFs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This one is constructing a huge big cell array to store all the data.
% zc_cell = cell(num_sen_levels,num_channels);
% h = waitbar(0,'Please wait. First cell array is being created ...');
% for i = 1:num_sen_levels
%     for ch = 1:num_channels
%         zc_cell{i,ch}=AN_FindData(an_data_files,i,ch);
%     end
%     waitbar(i /num_sen_levels)
% end
% close(h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is for taking all the data from that big cell array and assigning them according to the number of channels 
new_zc_cell = AN_Reconstruct_AssignData_eachChannel(an_data_files);

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
    %time_value = new_zc_cell{1,ch};
    time_value = delay_zc_cell{1,ch};
    new_time_value = time_value(:,1);
    if ch == 1
        max_Time_Value = max(new_time_value);
        min_Time_Value = min(new_time_value);
    end
    % Another array should be made to count the highest data in the whole set of data.
    max_new_timevalue = vertcat(max_Time_Value,new_time_value);
    min_new_timevalue = vertcat(min_Time_Value,new_time_value);
    max_Time_Value = max(max_new_timevalue);
    min_Time_Value = min(min_new_timevalue);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the part, which actually generates the signals for each channel.

if nargin==2
    multiplied_zc_cell = AN_Reconstruct_GenerateSignal_eachChannel_NoRamp_NoSpikeCap(an_data_files, delay_zc_cell, channel);
else
    multiplied_zc_cell = AN_Reconstruct_GenerateSignal_eachChannel_NoRamp_NoSpikeCap(an_data_files, delay_zc_cell);
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now there is different signals for all the different channels.



% % I am doing this work to save memory. 18/01/2011
% sum = zeros(1,length_sig);
% h4 = waitbar(0,'Please wait untill the Final signal are produced...');
% for ch = 1:num_channels
%     if ch == 1
%         sum = multiplied_zc_cell{1,ch};
%     else
%         disp(length(sum))
%         disp(length(multiplied_zc_cell{1,ch}))
%         sum = sum + multiplied_zc_cell{1,ch};
%     end
%     waitbar(ch/num_channels)
% end
% close(h4)
% Final_Signal = sum;




% A big structure is needed to hold all the produced signals for different channels
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
% data_axis = 0:(1/Fs):((length_sig/Fs))-(1/Fs);
% figure(2)
% plot(data_axis,Final_Signal,'r')
% plot(data_axis,Final_Signal,'b')
% The output signal should be saved in appropriate name.
new_Name = [an_data_files(1:end-7) '_NEW_NoRamp_NoSpikeCap'];
new_Name = [new_Name '.wav'];
wavwrite(Final_Signal,44100,16,new_Name);
final_signal(1,:) = final_signal(1,:)/max(abs(final_signal(1,:)));
wavwrite(final_signal(1,:),44100,16,'output_1.wav');
final_signal(45,:) = final_signal(45,:)/max(abs(final_signal(45,:)));
wavwrite(final_signal(45,:),44100,16,'output_45.wav');
disp(' ')
disp(' ')
disp('Re-Construction complete.');
disp('');
disp(['The audio file has been saved to ' new_Name]);
% End of 1

% % 2. This is the way where not all the signals are added up.
% for i = 1:(num_channels-5)
% %for i = ((num_channels/2)+1):num_channels
%     temp_final_signal = final_signal(i,:);
%     if i == 1
%         Final_Signal = temp_final_signal;
%     else
%         Final_Signal = Final_Signal + temp_final_signal;
%     end
%     disp('bingo');
% end
% Final_Signal = Final_Signal/max(abs(Final_Signal));
% Final_Signal = Final_Signal*0.9;
% Final_Signal = Final_Signal';
% data_axis = 0:(1/Fs):((length_sig/Fs))-(1/Fs);
% plot(data_axis,Final_Signal,'r')
% %plot(data_axis,Final_Signal,'b')
% % The output signal should be saved in appropriate name.
% new_Name = [an_data_files(1:end-7) '_NEW_Handicap'];
% new_Name = [new_Name '.wav'];
% wavwrite(Final_Signal,44100,16,new_Name);
% disp('Re-Construction complete.');
% disp('');
% disp(['The audio file has been saved to ' new_Name]);
% % End of 2

% The output file should consist of important data
%OutPut.zc_cell=zc_cell;
OutPut.new_zc_cell=new_zc_cell;
OutPut.delay_zc_cell=delay_zc_cell;
OutPut.multiplied_zc_cell=multiplied_zc_cell;
OutPut.Final_Signal=Final_Signal;
OutPut.max_Time_Value=max_Time_Value;
OutPut.min_Time_Value=min_Time_Value;

output_filename = [an_data_files(1:end-7) '_OutPut.mat'];
save(output_filename, 'OutPut')
disp(['The data has been saved to ' output_filename]);



end
