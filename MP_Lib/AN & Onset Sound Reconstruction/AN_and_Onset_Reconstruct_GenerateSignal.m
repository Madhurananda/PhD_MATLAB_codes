function [multiplied_zc_cell, new_AN_zc_cell, delay_zc_cell_AN, ... 
    delay_zc_cell_Onset, max_Time_Value_AN, min_Time_Value_AN, ... 
    max_Time_Value_Onset, min_Time_Value_Onset, num_channels, ...
    length_sig, PeakofSignal, Fs] = AN_and_Onset_Reconstruct_GenerateSignal ... 
    (an_data_files, Onset_data_files, Original_Onset_data_files, range)
% This one is with Ramp technique, but no spike Cap applied. 

% This Re construction work reproduces the signals from the spike.
%   It has only one output file.
%   The output file saves all the important data used to produced the
%   signal.

% The data should be loaded from the file
ANbig.data = load(an_data_files);
Onsetbig.data = load(Onset_data_files);
ori_Onsetbig.data = load(Original_Onset_data_files);

% Pick out some useful information
% How many channels for AN spikes?
num_channels    = ANbig(1).data.AN.channels;
% How many sen levels for AN spikes?
num_sen_levels  = ANbig(1).data.AN.iterations;
% How many channels for Onset spikes?
num_Onset_channels    = Onsetbig(1).data.ANparams.channels;
% How many sen levels for Onset spikes?
num_Onset_sen_levels  = Onsetbig(1).data.ANparams.iterations;

if (num_channels ~= num_Onset_channels)
    disp('The number of Channels are not same in AN and Onset data.'); 
    disp('Reconstruction is not possible.')
    exit;
elseif (num_sen_levels ~= num_Onset_sen_levels)
    disp('The number of Sensitivity levels are not same in AN and Onset data.');
    disp('Reconstruction is not possible.')
    exit;
end

% What is the length of the original input signal data?
length_sig = ANbig(1).data.AN.datalength;
% What is the threshhold level values?
thresh_levels = ANbig(1).data.AN.thresh_levels;
% What is the delay vector for this sound signal?
delayVector = ANbig(1).data.AN.delayVector;
% What is the Sampling Rate of the original input sound signal?
Fs = ANbig(1).data.AN.fs ;
% % What is the bit depth of the original input sound signal?
% nbits = ANbig(1).data.AN.nbits ;
% What is the peak of the original signal?
PeakofSignal = ANbig(1).data.AN.PeakofSignal;
% What are the center frequencies used for each channel by gammatone 
% filterbank?
cochCFs = ANbig(1).data.AN.cochCFs;


% Get the AN Spike trains according to each CHANNELS for AN dataset.  
new_AN_zc_cell = ANbig(1).data.AN.Spike_Assign_Channels;
% Get the Onset Spike trains according to each CHANNELS for Onset dataset.
new_Onset_zc_cell = Onsetbig(1).data.onsetparams.OnsetSpike_Channel;
% Get the Original OnsetSpike trains according to each CHANNELS for Onset 
% dataset.
new_ori_Onset_zc_cell = ori_Onsetbig(1).data.onsetparams.OnsetSpike_Channel;

%% Investigations and comparisons between AN and Onset spikes 
% % Ploting of both those spike trains at any random channel (say at 15)
% Ch = 2;
% xAxis = 0:(1/Fs):((length_sig/Fs)-(1/Fs));
% for i = 1:length(new_AN_zc_cell{1,Ch})
%     if (i==1)
%         xAxis_temp = 0:(1/Fs):new_AN_zc_cell{1,Ch}(i,1);
%         yAxis_temp = zeros(1,length(xAxis_temp)-1);
%         yAxis = horzcat(yAxis_temp, new_AN_zc_cell{1,Ch}(i,2));
%     else
%         xAxis_temp = new_AN_zc_cell{1,Ch}(i-1,1)+(1/Fs):(1/Fs):new_AN_zc_cell{1,Ch}(i,1);
%         yAxis_temp = zeros(1,length(xAxis_temp)-1);
%         yAxis = horzcat(yAxis,yAxis_temp, new_AN_zc_cell{1,Ch}(i,2));
%     end
% end
% xAxis_temp = new_AN_zc_cell{1,Ch}((length(new_AN_zc_cell{1,Ch})),1)+(1/Fs):(1/Fs):(length_sig/Fs);
% yAxis_temp = zeros(1,length(xAxis_temp)-1);
% AN_yAxis = horzcat(yAxis,yAxis_temp);
% 
% for i = 1:length(new_Onset_zc_cell{1,Ch})
%     if (i==1)
%         xAxis_temp = 0:(1/Fs):new_Onset_zc_cell{1,Ch}(i,1);
%         yAxis_temp = zeros(1,length(xAxis_temp)-1);
%         yAxis = horzcat(yAxis_temp, new_Onset_zc_cell{1,Ch}(i,2));
%     else
%         xAxis_temp = new_Onset_zc_cell{1,Ch}(i-1,1)+(1/Fs):(1/Fs):new_Onset_zc_cell{1,Ch}(i,1);
%         yAxis_temp = zeros(1,length(xAxis_temp)-1);
%         yAxis = horzcat(yAxis,yAxis_temp, new_Onset_zc_cell{1,Ch}(i,2));
%     end
% end
% xAxis_temp = new_Onset_zc_cell{1,Ch}((length(new_Onset_zc_cell{1,Ch})),1)+(1/Fs):(1/Fs):(length_sig/Fs);
% yAxis_temp = zeros(1,length(xAxis_temp)-1);
% Onset_yAxis = horzcat(yAxis,yAxis_temp);
% 
% % Due to Matlab computation or data error, the lengths of 'AN_yAxis' and
% % 'Onset_yAxis' are not same. 
% % This has been made equal forcefully. 
% if (length_sig ~= length(AN_yAxis))
%     if (length_sig > length(AN_yAxis))
%         temp_zeros = zeros(1,(length_sig-length(AN_yAxis)));
%         AN_yAxis = horzcat(AN_yAxis, temp_zeros);
%     elseif (length_sig < length(AN_yAxis))
%         AN_yAxis = AN_yAxis(1:length_sig);
%     end
% end
% 
% if (length_sig ~= length(Onset_yAxis))
%     if (length_sig > Onset_yAxis)
%         temp_zeros = zeros(1,(length_sig-length(Onset_yAxis)));
%         Onset_yAxis = horzcat(Onset_yAxis, temp_zeros);
%     elseif (length_sig < length(Onset_yAxis))
%         Onset_yAxis = Onset_yAxis(1:length_sig);
%     end
% end
% 
% figure
% plot(xAxis, AN_yAxis*0.2)
% hold on
% stem(xAxis, Onset_yAxis*0.2, 'r')
% hold on
% sig = wavread([an_data_files(1:end-13) '.wav']);
% plot(xAxis, sig, 'g')

% yAxis = zeros(1,length_sig);
% h = waitbar(0,'Please wait...');
% for i = 1:length_sig
%     for j = 1:length(new_AN_zc_cell{1,Ch})
%         if (xAxis(i) == new_AN_zc_cell{1,Ch}(j,1))
%             yAxis(i) = new_AN_zc_cell{1,Ch}(j,2);
%         end
%     end
%     waitbar(i / length_sig)
% end
% close(h) 
% figure
% plot(xAxis, yAxis)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Action should be taken for delay vectors for AN Spikes
delay_zc_cell_AN=cell(1,num_channels);
for ch = 1:num_channels
    time_value = new_AN_zc_cell{1,ch};
    if (~isempty(time_value))
        new_time_value = time_value(:,1);
        for i = 1:length(new_time_value)
            new_time_value(i) = new_time_value(i)-((1.0)*delayVector(ch));
        end
        time_value(:,1) = new_time_value;
        delay_zc_cell_AN{1,ch} = time_value;
    else
        new_time_value = zeros(1,1);
    end
    
end

% Action should be taken for delay vectors for AN Onset Spikes
delay_zc_cell_Onset=cell(1,num_Onset_channels);
for ch = 1:num_Onset_channels
    time_value = new_Onset_zc_cell{1,ch};
    if (~isempty(time_value))
        new_time_value = time_value(:,1);
        for i = 1:length(new_time_value)
            new_time_value(i) = new_time_value(i)-((1.0)*delayVector(ch));
        end
        time_value(:,1) = new_time_value;
        delay_zc_cell_Onset{1,ch} = time_value;
    else
        new_time_value = zeros(1,1);
    end
    
end

% Action should be taken for delay vectors for original Onset Spikes
delay_zc_cell_ori_Onset=cell(1,num_Onset_channels);
for ch = 1:num_Onset_channels
    time_value = new_ori_Onset_zc_cell{1,ch};
    if (~isempty(time_value))
        new_time_value = time_value(:,1);
        for i = 1:length(new_time_value)
            new_time_value(i) = new_time_value(i)-((1.0)*delayVector(ch));
        end
        time_value(:,1) = new_time_value;
        delay_zc_cell_ori_Onset{1,ch} = time_value;
    else
        new_time_value = zeros(1,1);
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The maximum value of the whole data set should be calculated here. 
for ch = 1:num_channels
    % -- This is where the delay vector effects can be altered. -- %
    % To deactivate the delay compensation, 'time_value' has to be taken from
    % new_zc_cell, not from 'delay_zc_cell'.
    %time_value = new_zc_cell{1,ch};
    time_value = delay_zc_cell_AN{1,ch};
    if (~isempty(time_value))
        new_time_value = time_value(:,1);
        if ch == 1
            max_Time_Value = max(new_time_value);
            min_Time_Value = min(new_time_value);
        end
        % Another array should be made to count the highest data in the 
        % whole set of data.
        max_new_timevalue = vertcat(max_Time_Value,new_time_value);
        min_new_timevalue = vertcat(min_Time_Value,new_time_value);
        max_Time_Value_AN = max(max_new_timevalue);
        min_Time_Value_AN = min(min_new_timevalue);
    else
        max_Time_Value_AN = 0;
        min_Time_Value_AN = 0;
    end
end

for ch = 1:num_Onset_channels
    % -- This is where the delay vector effects can be altered. -- %
    % To deactivate the delay compensation, 'time_value' has to be taken from
    % new_zc_cell, not from 'delay_zc_cell'.
    %time_value = new_zc_cell{1,ch};
    time_value = delay_zc_cell_Onset{1,ch};
    if (~isempty(time_value))
        new_time_value = time_value(:,1);
        if ch == 1
            max_Time_Value = max(new_time_value);
            min_Time_Value = min(new_time_value);
        end
        % Another array should be made to count the highest data in the 
        % whole set of data.
        max_new_timevalue = vertcat(max_Time_Value,new_time_value);
        min_new_timevalue = vertcat(min_Time_Value,new_time_value);
        max_Time_Value_Onset = max(max_new_timevalue);
        min_Time_Value_Onset = min(min_new_timevalue);
    else
        max_Time_Value_Onset = 0;
        min_Time_Value_Onset = 0;
    end
end

for ch = 1:num_Onset_channels
    % -- This is where the delay vector effects can be altered. -- %
    % To deactivate the delay compensation, 'time_value' has to be taken from
    % new_zc_cell, not from 'delay_zc_cell'.
    %time_value = new_zc_cell{1,ch};
    time_value = delay_zc_cell_ori_Onset{1,ch};
    if (~isempty(time_value))
        new_time_value = time_value(:,1);
        if ch == 1
            max_Time_Value = max(new_time_value);
            min_Time_Value = min(new_time_value);
        end
        % Another array should be made to count the highest data in the 
        % whole set of data.
        max_new_timevalue = vertcat(max_Time_Value,new_time_value);
        min_new_timevalue = vertcat(min_Time_Value,new_time_value);
        max_Time_Value_Onset = max(max_new_timevalue);
        min_Time_Value_Onset = min(min_new_timevalue);
    else
        max_Time_Value_Onset = 0;
        min_Time_Value_Onset = 0;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the part, which actually generates the signals for each channel.
% Reconstruction will be done in three parts. First, for the low
% frequencies, AN spikes will be used. For middle frequencies, Amplitude
% Modulated Onset spikes will be used. For the higher frequencies. only
% Onset spikes will be used. 
% For now: Low Freq.-0-206Hz(1st-5th); Middle Freq.-238-5211Hz(6th-41th);
% High Freq.-5606-10000(42nd-50th).
% The list of the cochlear central frequency is - 
% 100Hz             (1th	Channel)
% 123.8965Hz		(2th	Channel)
% 149.5295Hz		(3th	Channel)
% 177.0253Hz		(4th	Channel)
% 206.5192Hz		(5th	Channel)
% 238.1565Hz		(6th	Channel)
% 272.0928Hz		(7th	Channel)
% 308.4954Hz		(8th	Channel)
% 347.5433Hz		(9th	Channel)
% 389.4288Hz		(10th	Channel)
% 434.3582Hz		(11th	Channel)
% 482.5527Hz		(12th	Channel)
% 534.2494Hz		(13th	Channel)
% 589.703Hz		    (14th	Channel)
% 649.1865Hz		(15th	Channel)
% 712.9926Hz		(16th	Channel)
% 781.4355Hz		(17th	Channel)
% 854.8523Hz		(18th	Channel)
% 933.6042Hz		(19th	Channel)
% 1018.0791Hz		(20th	Channel)
% 1108.6929Hz		(21th	Channel)
% 1205.8915Hz		(22th	Channel)
% 1310.1537Hz		(23th	Channel)
% 1421.9927Hz		(24th	Channel)
% 1541.959Hz		(25th	Channel)
% 1670.6434Hz		(26th	Channel)
% 1808.6793Hz		(27th	Channel)
% 1956.7463Hz		(28th	Channel)
% 2115.5735Hz		(29th	Channel)
% 2285.9427Hz		(30th	Channel)
% 2468.6927Hz		(31th	Channel)
% 2664.7233Hz		(32th	Channel)
% 2874.9995Hz		(33th	Channel)
% 3100.5566Hz		(34th	Channel)
% 3342.5051Hz		(35th	Channel)
% 3602.036Hz		(36th	Channel)
% 3880.4272Hz		(37th	Channel)
% 4179.0493Hz		(38th	Channel)
% 4499.3723Hz		(39th	Channel)
% 4842.9734Hz		(40th	Channel)
% 5211.5442Hz		(41th	Channel)
% 5606.8992Hz		(42th	Channel)
% 6030.9848Hz		(43th	Channel)
% 6485.889Hz		(44th	Channel)
% 6973.8513Hz		(45th	Channel)
% 7497.274Hz		(46th	Channel)
% 8058.7342Hz		(47th	Channel)
% 8660.9959Hz		(48th	Channel)
% 9307.0244Hz		(49th	Channel)
% 10000Hz           (50th	Channel)

disp(' ')
if (range == 1)
    freq_range = [15 40];
    disp('The Frequency Range has been assigned for MALE')
elseif (range == 2)
    freq_range = [25 40];
    disp('The Frequency Range has been assigned for FEMALE')
elseif (range == 3)
    freq_range = [25 45];
    disp('The Frequency Range has been assigned for MUSICAL')
else
    disp('The range value has to be 1 or 2 or 3')
end

disp(' ')
disp('The frequency range has been distributed like this - ')
disp(['The Lower range as - ' num2str(cochCFs(1)) 'Hz frequency to ' ... 
    num2str(cochCFs(freq_range(1))) 'Hz frequency'])
disp(['The Middle range as - ' num2str(cochCFs(freq_range(1)+1)) ... 
    'Hz frequency to ' num2str(cochCFs(freq_range(2))) 'Hz frequency'])
disp(['The Upper range as - ' num2str(cochCFs(freq_range(2)+1)) ... 
    'Hz frequency to ' num2str(cochCFs(num_channels)) 'Hz frequency'])
disp(' ')
multiplied_zc_cell = cell(1,num_channels);
% Calculating the multiplier from the threshhold levels
temp_mult = 1/thresh_levels(length(thresh_levels));
mult_array = zeros(1,length(thresh_levels));
for i = 1:length(thresh_levels)
    mult_array(i) = temp_mult*thresh_levels(i);
end
% There should be a initial silence for each channel.
initial_sl = 0.05;
div = true; % This is boolean var which determines that a signal will be 
% normalised or not. 
TotalNoSpike = 0; % This is the total; number of spikes used to regenerate 
% the sound. 
h3 = waitbar(0,'Please wait untill the sine curves are produced...');
for ch = 1:num_channels
    
    
    %% for the frequency which is considered as low: we generate normal 
    % Auditory Nerve signal
    if ch<=freq_range(1)
        % ----------- Only AN spikes are to be used here --------
        % -- This is where the delay vector effects can be altered. -- %
        % To deactivate the delay compensation, 'time_value' has to be 
        % taken from new_zc_cell, not from 'delay_zc_cell'.
        %x_time_value = new_zc_cell{1,ch};
        x_time_value = delay_zc_cell_AN{1,ch};
        TotalNoSpike = TotalNoSpike + length(x_time_value(:,1)); 
        % Is there any spike at all?
        if (~isempty(x_time_value))
            % There SHOULD BE a quiet at the beginning of sound. 
            ind=[];
            cc1 = x_time_value(:,1);
            cc2 = x_time_value(:,2);
            for i = 1:length(cc1)
                if (cc1(i)) <= initial_sl
                    ind = horzcat(ind,i);
                end
            end
            cc1(ind) = [];
            cc2(ind) = [];
            clear time_value;
            [m, n] = size(cc1);
            if (m==1 && n==0)
                cc1 = cc1';
                cc2 = cc2';
            end
            time_value(:,1) = cc1;
            time_value(:,2) = cc2;
            
            new_time_value = time_value(:,1);
            occurrence = time_value(:,2);
            % Get the center-frequency of this channel.
            channel_frequency = cochCFs(ch);
            % So, each sine wave cycle will take (1/frequency) sec.
            cycle_time=1/channel_frequency;
            
            % If there is only one spike - 
            if (isempty(new_time_value))
                disp(['There is no AN Spike at Channel - ' num2str(ch)])
                new_Y=zeros(1,length_sig);
                div = false;
            elseif (length(new_time_value)==1)
                % Find out the starting time of the first cycle.
                start_cycle_time = new_time_value-cycle_time;
                % Find out the end time of the second cycle.
                end_cycle_time = new_time_value+cycle_time;

                % This is for the first cycle
                x_1=start_cycle_time+(1/Fs):(1/Fs):new_time_value;
                z_1=0+((2*pi)/length(x_1)):(2*pi)/length(x_1):(2*pi);
                y_1=sin(z_1);
                mult_1 = 0:((mult_array(occurrence)-0)/((length(x_1))-1)):... 
                    mult_array(occurrence);
                new_y_1 = y_1.*mult_1;

                % This is for the second cycle
                x_2=new_time_value+(1/Fs):(1/Fs):end_cycle_time;
                z_2=0+((2*pi)/length(x_2)):(2*pi)/length(x_2):(2*pi);
                y_2=sin(z_2);
                mult_2 = mult_array(occurrence):(-((mult_array... 
                    (occurrence)-0)/((length(x_2))-1))):0;
                new_y_2 = y_2.*mult_2;

                % Now add some zeros at the Beggining.
                zero_array_x_1 = 0:(1/Fs):start_cycle_time;
                zero_array_Y_1 = zeros(1,length(zero_array_x_1));
                % Now add some zeros at the End.
                zero_array_x_2 = end_cycle_time+(1/Fs):(1/Fs):... 
                    ((length_sig)/Fs);
                zero_array_Y_2 = zeros(1,length(zero_array_x_2));

                % Now concate everything in this channel.
                X=horzcat(zero_array_x_1,x_1,x_2,zero_array_x_2);
                new_Y=horzcat(zero_array_Y_1,new_y_1,new_y_2,... 
                    zero_array_Y_2); % 'new_Y' is the final signal for this 
                % channel
                % plot(X,new_Y)
            % So for this channel there are plenty spikes. Move ahead with the
            % normal reconstruction procedure....
            else
                % We are not implementing the Max Spike Rate technique in 
                % this AN and Onset reconstruction work. 
                new_Y = []; % This is the generated signal for this channel.
                for i = 1:length(new_time_value)-1
                   
                    % What is the differnce between this spike and next spike?
                    temp_length = new_time_value(i+1)-new_time_value(i);
                    % This is to see the frequency of this spike interval.
                    % So, the frequency of the signal is necessary is
                    new_channel_frequency = 1/temp_length;
                    % This is to find out the maximum channel frequency.
                    if i ==1
                        max_new_channel_frequency = new_channel_frequency;
                        min_new_channel_frequency = new_channel_frequency;
                    else
                        if new_channel_frequency > max_new_channel_frequency
                            max_new_channel_frequency = new_channel_frequency;
                        elseif new_channel_frequency < ... 
                                min_new_channel_frequency
                            min_new_channel_frequency = new_channel_frequency;
                        end
                    end
                    new_channel_frequency_array(ch,1) = ... 
                        max_new_channel_frequency;
                    new_channel_frequency_array(ch,2) = ... 
                        min_new_channel_frequency;
                    
                    % Is this is the very first spike?
                    if (i==1 && i~=(length(new_time_value)-1))
                        % It is necessary to check that there is enough 
                        % room to add
                        % another sine wave. If not then, just add zeros...
                        if ((new_time_value(i)-temp_length)<0)
            %                 disp('This time here is no extra sine wave.')
                            zero_array_x = 0:(1/Fs):new_time_value(i);
                            zero_array_Y = zeros(1,length(zero_array_x));
%                             X=horzcat(zero_array_x,x);
%                             Y=horzcat(zero_array_Y,y);
                            % 'zero_array_Y' is the extra zeros..
                            new_Y=horzcat(zero_array_Y,new_Y); 
                        else
                            % It is necessary to add some zeros at the 
                            % beginning of each channel signal.
                            zero_array_x = 0:(1/Fs):(new_time_value(i)-... 
                                temp_length);
                            zero_array_y = zeros(1,length(zero_array_x));
                            %zero_array_Y = zero_array_y.*mult;
                            zero_array_Y = zeros(1,length(zero_array_x));
%                             X=horzcat(zero_array_x,X);
%                             Y=horzcat(zero_array_y,Y);
                            % new_Y = horzcat(zero_array_Y,new_Y);
                            
                            % One extra sine wave should be added at the 
                            % beggining.
                            temp_x = (new_time_value(i)-temp_length)+(1/Fs):... 
                                (1/Fs):new_time_value(i);
            %                 disp('This length has been added for EXTRA SINE WAVE')
            %                 disp(length(temp_x))
                            temp_z=0+((2*pi)/length(temp_x)):(2*pi)/... 
                                length(temp_x):(2*pi);
                            temp_y = sin(temp_z);
                            mult = 0:(mult_array(occurrence(i))/... 
                                ((length(temp_x))-1)):mult_array(occurrence(i));
                            temp_new_y = temp_y.*mult;
%                             X=horzcat(temp_x,x);
%                             Y=horzcat(temp_y,y);
                            % 'temp_new_y' is the extra wave
                            new_Y=horzcat(zero_array_Y, temp_new_y, new_Y); 
                        end
                    end
                    
                    % Now generate the Amplititude Modulated signal and 
                    % concate with 'new_Y', but check huge difference 
                    % between two spikes.
                    spike_time_diff = temp_length;
                    % What is a reasonable gap?
                    % It shold be the normal cyle time of this chennel 
                    % frequency
                    spike_gap = cycle_time;
                    % Check if there is a big gap or not.
                    % Do extra things if the diffrence is 2 times greater
                    % than the spike time diffrence.
                    if (spike_time_diff>(spike_gap*2))
                        % diff = spike_time_diff-(spike_gap);
                        % diff_a = 0:(1/Fs):spike_time_diff;
                        % Generate the first cycle at the beggining
                        % Find out the starting time of the first cycle.
                        start_cycle_time = new_time_value(i);
                        % Find out the end time of the first cycle.
                        end_cycle_time = new_time_value(i)+cycle_time;
                        % This is for the first cycle
                        x_1=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                        z_1=0+((2*pi)/length(x_1)):(2*pi)/length(x_1):(2*pi);
                        y_1=sin(z_1);
                        mult_1 = mult_array(occurrence(i)):... 
                            (-((mult_array(occurrence(i))-0)/((length(x_1))... 
                            -1))):0;
                        new_y_1 = y_1.*mult_1;

                        % Generate second cycle at the end
                        % Find out the starting time of the second cycle.
                        start_cycle_time = new_time_value(i+1)-cycle_time;
                        % Find out the end time of the second cycle.
                        end_cycle_time = new_time_value(i+1);
                        % This is for the first cycle
                        x_3=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                        z_3=0+((2*pi)/length(x_3)):(2*pi)/length(x_3):(2*pi);
                        y_3=sin(z_3);
                        mult_3 = 0:((mult_array(occurrence(i+1))-0)/... 
                            ((length(x_3))-1)):mult_array(occurrence(i+1));
                        new_y_3 = y_3.*mult_3;

                        % Generate some zeros at the middle if there is enough space
                        len_gen = length(new_y_1)+length(new_y_3);
                        len_shouldBe = new_time_value(i)+(1/Fs):(1/Fs):... 
                            new_time_value(i+1);
%                             start_zero_time=end_cycle_time;
%                             end_zero_time = new_time_value(i+1)-cycle_time;
%                             x_2=start_zero_time+(1/Fs):(1/Fs):end_zero_time;
                        new_y_2=zeros(1,(abs(len_gen - length(len_shouldBe))));

                        % Check the length - 
                        len_gen = length(new_y_1)+length(new_y_2)+... 
                            length(new_y_3);
                        if (len_gen ~= length(len_shouldBe))
                            disp('Something wrong')
                            disp(abs(len_gen - length(len_shouldBe)))
                            disp(' ')
                        end

%                             X=horzcat(X,x_1,x_2,x_3);
%                             Y=horzcat(Y,y_1,new_y_2,y_3);
                        new_Y=horzcat(new_Y,new_y_1,new_y_2,new_y_3);
                    else
                        x=new_time_value(i)+(1/Fs):(1/Fs):new_time_value(i+1);
                        z=0+((2*pi)/length(x)):(2*pi)/length(x):(2*pi);
                        y=sin(z);
                        if mult_array(occurrence(i)) == ... 
                                mult_array(occurrence(i+1))
                            mult = mult_array(occurrence(i+1));
                            new_y = y*mult;
                        else
                            if mult_array(occurrence(i)) < ... 
                                    mult_array(occurrence(i+1))
                                mult = mult_array(occurrence(i)):... 
                                    ((mult_array(occurrence(i+1))- ... 
                                    mult_array(occurrence(i)))/((length(x))... 
                                    -1)):mult_array(occurrence(i+1));
                                new_y = y.*mult;
                            else
                                mult = mult_array(occurrence(i)):... 
                                    (-((mult_array(occurrence(i))- ... 
                                    mult_array(occurrence(i+1)))/... 
                                    ((length(x))-1))):mult_array... 
                                    (occurrence(i+1));
                                new_y = y.*mult;
                            end
                        end
                        new_Y=horzcat(new_Y,new_y);
                    end
                    
                     % Is this is the very last spike?
                    if (i==(length(new_time_value)-1) && i~=1)
                        % Calculate how much space is available at the End
                        end_space = ((length_sig)/Fs)-new_time_value(i+1);
                        % If that spave is big enough, add the sine cycle 
                        % at the channel frequency
                        if (end_space>=cycle_time)
                            % Find out the starting time of the first cycle.
                            start_cycle_time = new_time_value(i+1);
                            % Find out the end time of the second cycle.
                            end_cycle_time = new_time_value(i+1)+cycle_time;
                            x=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                            z=0+((2*pi)/length(x)):(2*pi)/length(x):(2*pi);
                            y=sin(z);
                            mult = mult_array(occurrence(i+1)):...
                                (-((mult_array(occurrence(i+1))-0)/... 
                                ((length(x))-1))):0;
                            new_y = y.*mult;
%                             X=horzcat(X,x);
                            new_Y = horzcat(new_Y,new_y);
                            % Just add zeros upto end of the signal.
                            zero_array_x = length_sig-length(new_Y);
                            zero_array_y = zeros(1,zero_array_x);
                            new_Y = horzcat(new_Y,zero_array_y);
                        % If not, add the sine wave at the length equal to 
                        % the space left at the end.
                        else
                            disp('Hey, This is a rare case. Investigate through it.')
                            disp('Channel')
                            disp(ch)
                            new_cycle_time = end_space;
                            % Find out the starting time of the first cycle.
                            start_cycle_time = new_time_value(i+1);
                            % Find out the end time of the second cycle.
                            end_cycle_time = new_time_value(i+1)+new_cycle_time;
                            x=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                            z=0+((2*pi)/length(x)):(2*pi)/length(x):(2*pi);
                            y=sin(z);
                            mult = mult_array(occurrence(i+1)):... 
                                (-((mult_array(occurrence(i+1))-0)/...
                                ((length(x))-1))):0;
                            new_y = y.*mult;
%                             X=horzcat(X,x);
                            new_Y = horzcat(new_Y,new_y);
                        end
                    end
                   
                end
            end
            
            if (div)
                new_Y = 0.9*(new_Y/max(abs(new_Y)));
            else
                div = true;
            end
            
            % Find out the length of the produced signal and the length of the
            % original signal is same or not.
            compensate_zero=(length(new_Y)-length_sig);
            if (compensate_zero==0)
                % This is not necessary to see everytime.
                % disp('true')
            else
                disp('False. The signal length is not same.')
                disp('channel')
                disp(ch)
                disp('The missing number of values are .....')
                disp(compensate_zero)
                % zero_array_x = 0:(1/Fs):new_time_value(i);
                zero_array_Y = zeros(1,abs(compensate_zero));
                % X=horzcat(zero_array_x,x);
                % Y=horzcat(zero_array_Y,y);
                new_Y=horzcat(zero_array_Y,new_Y);
                disp(' ')
                disp('-----------------------------------------------------')
                disp(' ')
            end
            
        else
            new_Y=zeros(1,length_sig);
        end
        
        
        
    %% for the frequency which is considered as middle: we generate 
    % Amplitude Modulated Signal
    elseif ch>=(freq_range(1)+1) && ch<=freq_range(2)
        % ----------- Only AN_Onset spikes are to be used here --------
        % -- This is where the delay vector effects can be altered. -- %
        % To deactivate the delay compensation, 'time_value' has to be 
        % taken from new_Onset_zc_cell, not from 'delay_zc_cell_Onset'.
        % x_time_value = new_Onset_zc_cell{1,ch};
        x_time_value = delay_zc_cell_Onset{1,ch};
        TotalNoSpike = TotalNoSpike + length(x_time_value(:,1)); 
        
        if (~isempty(x_time_value))
            % There SHOULD BE a quiet at the beginning of sound. 
            ind=[];
            cc1 = x_time_value(:,1);
            cc2 = x_time_value(:,2);
            for i = 1:length(cc1)
                if (cc1(i)) <= initial_sl
                    ind = horzcat(ind,i);
                end
            end
            cc1(ind) = [];
            cc2(ind) = [];
            clear time_value;
            [m, n] = size(cc1);
            if (m==1 && n==0)
                cc1 = cc1';
                cc2 = cc2';
            end
            time_value(:,1) = cc1;
            time_value(:,2) = cc2;

            % Get the Onset times
            new_time_value = time_value(:,1);
            % Get the Onset occurrences
            occurrence = time_value(:,2);
            % Get the center-frequency of this channel.
            channel_frequency = cochCFs(ch);
            % So, each sine wave cycle will take (1/frequency) sec.
            cycle_time=1/channel_frequency;
            %     disp('Channel')
            %     disp(ch)


            % If there is only one Onset, then ....
            if (length(new_time_value)==1 || isempty(new_time_value))
                disp(['only One or Zero onset. There will not be any signal generated at Channel - ' num2str(ch)])
                new_Y=zeros(1,length_sig);
                div = false;
            % So for this channel there are many onset spikes. Move ahead with the
            % normal reconstruction procedure....    
            else
                % We are not implementing the Max Spike Rate technique in this AN
                % and Onset reconstruction work. 

                new_Y = []; % This is the generated signal for this channel.
                for i = 1:length(new_time_value)-1
                    % What is the differnce between this spike and next spike?
                    temp_length = new_time_value(i+1)-new_time_value(i);
                    % This is to see the frequency of this spike interval.
                    % So, the frequency of the signal is necessary in this case is
                    new_channel_frequency = 1/temp_length;
                    % This is to find out the maximum channel frequency.
                    if i ==1
                        max_new_channel_frequency = new_channel_frequency;
                        min_new_channel_frequency = new_channel_frequency;
                    else
                        if new_channel_frequency > max_new_channel_frequency
                            max_new_channel_frequency = new_channel_frequency;
                        elseif new_channel_frequency < min_new_channel_frequency
                            min_new_channel_frequency = new_channel_frequency;
                        end
                    end
                    new_channel_frequency_array(ch,1) = ... 
                        max_new_channel_frequency;
                    new_channel_frequency_array(ch,2) = ... 
                        min_new_channel_frequency;

                    % Modulated frequency :
                    % Assuming the gap between two onsets is 8 ms, the base
                    % modulated frequency comes upto 125 Hz. But, it has to
                    % be adjusted so the the regerated signal ends at the
                    % beginning of the next spike. 
                    base_mod_freq = 125; 
                    % Carrier frequency is this channel's centre frequency. 
                    base_crr_freq = (channel_frequency);
                    mod_freq = (round(base_mod_freq*temp_length))/temp_length;
                    crr_freq = (round(base_crr_freq*temp_length))/temp_length;
                    a1 = crr_freq;
                    a2 = crr_freq - mod_freq;

                    % Is this is the very first spike?
%                     if (i==1 && i~=(length(new_time_value)-1))
                    if (i==1)
                        % It is necessary to check that there is enough room to add
                        % another AM signal. If not then, just add zeros...
                        if ((new_time_value(i)-cycle_time)<0)
                            % Then, It is necessary to add some signal at 
                            % the Beginning.
                            zero_array_x = 0:(1/Fs):new_time_value(i);
                            zero_array_Y = zeros(1,length(zero_array_x));
                            % 'zero_array_Y' is the extra zeros..
                            new_Y=horzcat(zero_array_Y,new_Y); 
                        else
                            % It is necessary to add some zeros at the 
                            % beginning of each channel signal.
                            zero_array_x = 0:(1/Fs):(new_time_value(i)- ... 
                                cycle_time);
                            zero_array_y = zeros(1,length(zero_array_x));
                            zero_array_Y = zeros(1,length(zero_array_x));
                            
                            % One extra sine wave should be added at the beggining.
                            temp_x = (new_time_value(i)-cycle_time)+(1/Fs):...
                                (1/Fs):new_time_value(i);
                            temp_z=0+((2*pi)/length(temp_x)):(2*pi)/...
                                length(temp_x):(2*pi);
                            temp_y = sin(temp_z);
                            mult = 0:(mult_array(occurrence(i))/... 
                                ((length(temp_x))-1)):mult_array(occurrence(i));
                            temp_new_y = temp_y.*mult;
                            % 'temp_new_y' is the extra wave
                            new_Y=horzcat(zero_array_Y, temp_new_y, new_Y); 
                        end
                    end
                    
                    
                    % Now generate the Amplititude Modulated signal and 
                    % concate with 'new_Y', but
                    % remember to check huge difference between two spikes. 
                    spike_time_diff = temp_length;
                    % What is a reasonable gap?
                    % It shold be the normal cycle time of this chennel frequency
                    spike_gap = 1/base_mod_freq;
                    % Check if there is a big gap or not.
                    if (spike_time_diff>(spike_gap*6))
                        % diff = spike_time_diff-(spike_gap);
                        % diff_a = 0:(1/Fs):spike_time_diff;
                        % Generate the first cycle at the beggining
                        % Find out the starting time of the first cycle.
                        start_cycle_time = new_time_value(i);
                        % Find out the end time of the first cycle.
                        end_cycle_time = new_time_value(i)+spike_gap*3;
                        % This is for the first cycle
                        x_1=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                        a1 = crr_freq;
                        a2 = crr_freq - mod_freq;
                        t = (0):(1/Fs):((length(x_1)/Fs)-(1/Fs));
                        x1 = 0.8*sin(2*pi*a1*t);
                        x2 = 0.4*sin(2*pi*a2*t);
                        y_1 = x1-x2;

                        mult_1 = mult_array(occurrence(i)):... 
                            (-((mult_array(occurrence(i))-0)/... 
                            ((length(x_1))-1))):0;
                        new_y_1 = y_1.*mult_1;

                        % Generate second cycle at the end
                        % Find out the starting time of the second cycle.
                        start_cycle_time = new_time_value(i+1)-3*spike_gap;
                        % Find out the end time of the second cycle.
                        end_cycle_time = new_time_value(i+1);
                        % This is for the last cycle
                        x_3=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                        a1 = crr_freq;
                        a2 = crr_freq - mod_freq;
                        t = 0:(1/Fs):((length(x_3)/Fs)-(1/Fs));
                        x1 = 0.8*sin(2*pi*a1*t);
                        x2 = 0.4*sin(2*pi*a2*t);
                        y_3 = x1-x2;
                        mult_3 = 0:((mult_array(occurrence(i+1))-0)/... 
                            ((length(x_3))-1)):mult_array(occurrence(i+1));
                        new_y_3 = y_3.*mult_3;
                        
                        
                        % Generate some zeros at the middle if there is enough space
                        len_gen = length(new_y_1)+length(new_y_3);
                        len_shouldBe = new_time_value(i)+(1/Fs):(1/Fs):... 
                            new_time_value(i+1);
                        new_y_2=zeros(1,(abs(len_gen - length(len_shouldBe))));

                        % Check the length - 
                        len_gen = length(new_y_1)+length(new_y_2)+ ... 
                            length(new_y_3);
                        if (len_gen ~= length(len_shouldBe))
                            disp('Something wrong')
                            disp(abs(len_gen - length(len_shouldBe)))
                            disp(' ')
                        end
                        new_Y=horzcat(new_Y,new_y_1,new_y_2,new_y_3);
                        
                    else
                        x=new_time_value(i)+(1/Fs):(1/Fs):new_time_value(i+1);
                        t = 0+(1/Fs):(1/Fs):(length(x)/Fs);
                        x1 = 0.8*sin(2*pi*a1*t);
                        x2 = 0.4*sin(2*pi*a2*t);
                        y = x1-x2;

                        % z=0+((2*pi)/length(x)):(2*pi)/length(x):(2*pi);
                        % y=sin(z);
                        if mult_array(occurrence(i)) == ... 
                                mult_array(occurrence(i+1))
                            mult = mult_array(occurrence(i+1));
                            new_y = y*mult;
                        else
                            if mult_array(occurrence(i)) < ... 
                                    mult_array(occurrence(i+1))
                                mult = mult_array(occurrence(i)):... 
                                    ((mult_array(occurrence(i+1))- ... 
                                    mult_array(occurrence(i)))/((... 
                                    length(x))-1)):mult_array(occurrence(i+1));
                                new_y = y.*mult;
                            else
                                mult = mult_array(occurrence(i)):... 
                                    (-((mult_array(occurrence(i))- ... 
                                    mult_array(occurrence(i+1)))/(( ... 
                                    length(x))-1))):mult_array(occurrence(i+1));
                                new_y = y.*mult;
                            end
                        end
                        new_Y=horzcat(new_Y,new_y);

                    end


                    % Is this is the very last spike?
                    if (i==(length(new_time_value)-1))
                        
                        % Calculate how much space is available at the End
                        end_space = ((length_sig)/Fs)-new_time_value(i+1);
                        % If that spave is big enough, add the sine cycle 
                        % at the channel frequency
                        if (end_space>=cycle_time)
                            % Find out the starting time of the first cycle.
                            start_cycle_time = new_time_value(i+1);
                            % Find out the end time of the second cycle.
                            end_cycle_time = new_time_value(i+1)+cycle_time;
                            x=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                            z=0+((2*pi)/length(x)):(2*pi)/length(x):(2*pi);
                            y=sin(z);
                            mult = mult_array(occurrence(i+1)):(-((... 
                                mult_array(occurrence(i+1))-0)/(( ... 
                                length(x))-1))):0;
                            new_y = y.*mult;
                            new_Y = horzcat(new_Y,new_y);
                            % Just add zeros upto end of the signal.
                            zero_array_x = length_sig-length(new_Y);
                            zero_array_y = zeros(1,zero_array_x);
                            new_Y = horzcat(new_Y,zero_array_y);
                        else
                            disp('Hey, This is a rare case. Investigate through it.')
                            disp('Channel')
                            disp(ch)
                            new_cycle_time = end_space;
                            % Find out the starting time of the first cycle.
                            start_cycle_time = new_time_value(i+1);
                            % Find out the end time of the second cycle.
                            end_cycle_time = new_time_value(i+1)+new_cycle_time;
                            x=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                            z=0+((2*pi)/length(x)):(2*pi)/length(x):(2*pi);
                            y=sin(z);
                            mult = mult_array(occurrence(i+1)):(-(( ... 
                                mult_array(occurrence(i+1))-0)/(( ... 
                                length(x))-1))):0;
                            new_y = y.*mult;
%                             X=horzcat(X,x);
                            new_Y = horzcat(new_Y,new_y);
                        end
                        
                    end

                end
            end
            
            if (div)
                new_Y = 0.5*(new_Y/max(abs(new_Y)));
            else
                div = true;
            end
            
            
            % Find out the length of the produced signal and the length of the
            % original signal is same or not.
            compensate_zero=(length(new_Y)-length_sig);
            if (compensate_zero==0)
                % This is not necessary to see everytime.
                % disp('true')
            else
                disp('False. The signal length is not same.')
                disp('channel')
                disp(ch)
                disp('The missing number of values are .....')
                disp(compensate_zero)
                % zero_array_x = 0:(1/Fs):new_time_value(i);
                zero_array_Y = zeros(1,compensate_zero);
                % X=horzcat(zero_array_x,x);
                % Y=horzcat(zero_array_Y,y);
                new_Y=horzcat(zero_array_Y,new_Y);
                disp(' ')
                disp('------------------------------------------------------')
                disp(' ')
            end
        else
            disp('no onset at all. Investigate .. ')
            disp ('Channel')
            disp(ch)
            new_Y=zeros(1,length_sig);
        end
        
        
    
    
    %% for the frequency which is considered as high:
    % The strategy for the high frequency is that to generate some
    % white noise for short period of time (one cycle of centre frequncy)
    elseif ch>=(freq_range(2)+1) && ch<=num_channels
        
        % ----------- Only Original Onset spikes are to be used here --------
        % -- This is where the delay vector effects can be altered. -- %
        % To deactivate the delay compensation, 'time_value' has to be 
        % taken from new_Onset_zc_cell, not from 'delay_zc_cell_Onset'.
        % x_time_value = new_ori_Onset_zc_cell{1,ch};
        x_time_value = delay_zc_cell_ori_Onset{1,ch};
        TotalNoSpike = TotalNoSpike + length(x_time_value(:,1)); 
        
        if (~isempty(x_time_value))
            % There SHOULD BE a quiet at the beginning of sound. 
            ind=[];
            cc1 = x_time_value(:,1);
            cc2 = x_time_value(:,2);
            for i = 1:length(cc1)
                if (cc1(i)) <= initial_sl
                    ind = horzcat(ind,i);
                end
            end
            cc1(ind) = [];
            cc2(ind) = [];
            clear time_value;
            [m, n] = size(cc1);
            if (m==1 && n==0)
                cc1 = cc1';
                cc2 = cc2';
            end
            time_value(:,1) = cc1;
            time_value(:,2) = cc2;

            % Get the Onset times
            new_time_value = time_value(:,1);
            % Get the Onset occurrences
            occurrence = time_value(:,2);
            % Get the center-frequency of this channel.
            channel_frequency = cochCFs(ch);
            % So, each sine wave cycle will take (1/frequency) sec.
            cycle_time=1/channel_frequency;
            %     disp('Channel')
            %     disp(ch)


            % If there is only one Onset, then ....
            if (length(new_time_value)==1 || isempty(new_time_value))
                disp('only One or Zero onset');
                disp(['There will not be any signal generated at Channel - ' num2str(ch)]);
                new_Y=zeros(1,length_sig);
                div = false;

            % So for this channel there are many onset spikes. Move ahead with the
            % normal reconstruction procedure....    
            else
                % We are not implementing the Max Spike Rate technique in this AN
                % and Onset reconstruction work. 

                new_Y = []; % This is the generated signal for this channel.
                for i = 1:length(new_time_value)-1
                    % What is the differnce between this spike and next spike?
                    temp_length = new_time_value(i+1)-new_time_value(i);
                    % This is to see the frequency of this spike interval.
                    % So, the frequency of the signal is necessary in this case is
                    new_channel_frequency = 1/temp_length;
                    % This is to find out the maximum channel frequency.
                    if i ==1
                        max_new_channel_frequency = new_channel_frequency;
                        min_new_channel_frequency = new_channel_frequency;
                    else
                        if new_channel_frequency > max_new_channel_frequency
                            max_new_channel_frequency = new_channel_frequency;
                        elseif new_channel_frequency < min_new_channel_frequency
                            min_new_channel_frequency = new_channel_frequency;
                        end
                    end
                    new_channel_frequency_array(ch,1) = max_new_channel_frequency;
                    new_channel_frequency_array(ch,2) = min_new_channel_frequency;

%                     % Modulated frequency :
%                     % Assuming the gap between two onsets is 8 ms, the base
%                     % modulated frequency comes upto 125 Hz. But, it has to
%                     % be adjusted so the the regerated signal ends at the
%                     % beginning of the next spike. 
%                     base_mod_freq = 125; 
%                     % Carrier frequency is this channel's centre frequency. 
%                     base_crr_freq = (channel_frequency);
%                     mod_freq = (round(base_mod_freq*temp_length))/temp_length;
%                     crr_freq = (round(base_crr_freq*temp_length))/temp_length;
%     %                     if (mod_freq>crr_freq)
%     %                         disp(ch)
%     %                         disp(temp_length)
%     %                         disp('xxxxxxxxxxxx Modulated frequency is greater than carrier. xxxxxxxxxx')
%     %                     end
%                     a1 = crr_freq;
%                     a2 = crr_freq - mod_freq;

                    % Is this is the very first spike?
                    if (i==1)
                        % It is necessary to check that there is enough 
                        % room to add another AM signal. If not then, just 
                        % add zeros...
                        if ((new_time_value(i)-cycle_time)<0)
                            % Then, It is necessary to add some signal at 
                            % the Beginning.
                            zero_array_x = 0:(1/Fs):new_time_value(i);
                            zero_array_Y = zeros(1,length(zero_array_x));
%                             X=horzcat(zero_array_x,x);
%                             Y=horzcat(zero_array_Y,y);
                            % 'zero_array_Y' is the extra zeros..
                            new_Y=horzcat(zero_array_Y,new_Y); 
                        else
                            % It is necessary to add some zeros at the beginning of each channel signal.
                            zero_array_x = 0:(1/Fs):(new_time_value(i)- ... 
                                cycle_time);
            %                 disp('This length of ZEROS has been added before EXTRA SINE WAVE')
            %                 disp(length(zero_array_x))
                            zero_array_y = zeros(1,length(zero_array_x));
                            %zero_array_Y = zero_array_y.*mult;
                            zero_array_Y = zeros(1,length(zero_array_x));
%                             X=horzcat(zero_array_x,X);
%                             Y=horzcat(zero_array_y,Y);
                            % new_Y = horzcat(zero_array_Y,new_Y);
                            
                            % One extra sine wave should be added at the beggining.
                            temp_x = (new_time_value(i)-cycle_time)+(1/Fs)... 
                                :(1/Fs):new_time_value(i);
            %                 disp('This length has been added for EXTRA SINE WAVE')
            %                 disp(length(temp_x))
                            temp_z=0+((2*pi)/length(temp_x)):(2*pi)/ ... 
                                length(temp_x):(2*pi);
                            temp_y = sin(temp_z);
                            mult = 0:(mult_array(occurrence(i))/(( ... 
                                length(temp_x))-1)):mult_array(occurrence(i));
                            temp_new_y = temp_y.*mult;
%                             X=horzcat(temp_x,x);
%                             Y=horzcat(temp_y,y);
                            % 'temp_new_y' is the extra wave
                            new_Y=horzcat(zero_array_Y, temp_new_y, new_Y); 
                        end
                    end
                    
                   
                    
                    % Now generate the white noise signal and concate with 'new_Y', but
                    % remember to check huge difference between two spikes. 
                    spike_time_diff = temp_length;
                    % What is a reasonable gap?
                    % It shold be the normal cycle time of this chennel frequency
                    spike_gap = cycle_time;
                    % Check if there is a big gap or not.
                    if (spike_time_diff>(spike_gap*3))
                        % diff = spike_time_diff-(spike_gap);
                        % diff_a = 0:(1/Fs):spike_time_diff;
                        % Generate the first cycle at the beggining
                        % Find out the starting time of the first cycle.
                        start_cycle_time = new_time_value(i);
                        % Find out the end time of the first cycle.
                        end_cycle_time = new_time_value(i)+spike_gap;
                        % This is for the first cycle
                        x_1=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                        
                        % Additive white Gaussian noise Coefficients
                        c=randn(length(x_1),1); 
                        a=1;
                        b=[1 1];
                        new_y=filter(b,a,c);
                        y_1 = new_y';
                        
%                         a1 = crr_freq;
%                         a2 = crr_freq - mod_freq;
%                         t = (0):(1/Fs):((length(x_1)/Fs)-(1/Fs));
%                         x1 = 0.8*sin(2*pi*a1*t);
%                         x2 = 0.4*sin(2*pi*a2*t);
%                         y_1 = x1-x2;

                        mult_1 = mult_array(occurrence(i)):(-((mult_array ... 
                            (occurrence(i))-0)/((length(x_1))-1))):0;
                        new_y_1 = y_1.*mult_1;

                        % Generate second cycle at the end
                        % Find out the starting time of the second cycle.
                        start_cycle_time = new_time_value(i+1)-spike_gap;
                        % Find out the end time of the second cycle.
                        end_cycle_time = new_time_value(i+1);
                        % This is for the last cycle
                        x_3=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                        % Additive white Gaussian noise Coefficients
                        c=randn(length(x_3),1); 
                        a=1;
                        b=[1 1];
                        new_y=filter(b,a,c);
                        y_3 = new_y';
                        mult_3 = 0:((mult_array(occurrence(i+1))-0)/((... 
                            length(x_3))-1)):mult_array(occurrence(i+1));
                        new_y_3 = y_3.*mult_3;
                        
                        
                        % Generate some zeros at the middle if there is enough space
                        len_gen = length(new_y_1)+length(new_y_3);
                        len_shouldBe = new_time_value(i)+(1/Fs):(1/Fs):... 
                            new_time_value(i+1);
%                             start_zero_time=end_cycle_time;
%                             end_zero_time = new_time_value(i+1)-cycle_time;
%                             x_2=start_zero_time+(1/Fs):(1/Fs):end_zero_time;
                        new_y_2=zeros(1,(abs(len_gen - length(len_shouldBe))));

                        % Check the length - 
                        len_gen = length(new_y_1)+length(new_y_2)+length(new_y_3);
                        if (len_gen ~= length(len_shouldBe))
                            disp('Something wrong')
                            disp(abs(len_gen - length(len_shouldBe)))
                            disp(' ')
                        end

%                             X=horzcat(X,x_1,x_2,x_3);
%                             Y=horzcat(Y,y_1,new_y_2,y_3);
                        new_Y=horzcat(new_Y,new_y_1,new_y_2,new_y_3);
                        
                    else
                        
                         % Now generate the white noise for the length of 
                         % '1/cochlea_center_frequency'
                        x=new_time_value(i)+(1/Fs):(1/Fs):new_time_value(i+1);
                        % Additive white Gaussian noise Coefficients
                        c=randn(length(x),1); 
                        a=1;
                        b=[1 1];
                        new_y=filter(b,a,c);
                        y = new_y';

                        % z=0+((2*pi)/length(x)):(2*pi)/length(x):(2*pi);
                        % y=sin(z);
                        if mult_array(occurrence(i)) == mult_array(occurrence(i+1))
                            mult = mult_array(occurrence(i+1));
                            new_y = y*mult;
                        else
                            if mult_array(occurrence(i)) < mult_array(occurrence(i+1))
                                mult = mult_array(occurrence(i)):((... 
                                    mult_array(occurrence(i+1))- ... 
                                    mult_array(occurrence(i)))/((... 
                                    length(x))-1)):mult_array(occurrence(i+1));
                                new_y = y.*mult;
                            else
                                mult = mult_array(occurrence(i)):(-((... 
                                    mult_array(occurrence(i))-mult_array(... 
                                    occurrence(i+1)))/((length(x))-1))):... 
                                    mult_array(occurrence(i+1));
                                new_y = y.*mult;
                            end
                        end
                        new_Y=horzcat(new_Y,new_y);

                    end


                    % Is this is the very last spike?
                    if (i==(length(new_time_value)-1))
                        
                        % Calculate how much space is available at the End
                        end_space = ((length_sig)/Fs)-new_time_value(i+1);
                        % If that spave is big enough, add the sine cycle 
                        % at the channel frequency
                        if (end_space>=cycle_time)
                            % Find out the starting time of the first cycle.
                            start_cycle_time = new_time_value(i+1);
                            % Find out the end time of the second cycle.
                            end_cycle_time = new_time_value(i+1)+cycle_time;
                            x=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                            z=0+((2*pi)/length(x)):(2*pi)/length(x):(2*pi);
                            y=sin(z);
                            mult = mult_array(occurrence(i+1)):(-((... 
                                mult_array(occurrence(i+1))-0)/((length(x))... 
                                -1))):0;
                            new_y = y.*mult;
%                             X=horzcat(X,x);
                            new_Y = horzcat(new_Y,new_y);
                            % Just add zeros upto end of the signal.
                            zero_array_x = length_sig-length(new_Y);
                            zero_array_y = zeros(1,zero_array_x);
                            new_Y = horzcat(new_Y,zero_array_y);
                        % If not, add the sine wave at the length equal to 
                        % the space left at the end.
                        else
                            disp('Hey, This is a rare case. Investigate ...')
                            disp('Channel')
                            disp(ch)
                            new_cycle_time = end_space;
                            % Find out the starting time of the first cycle.
                            start_cycle_time = new_time_value(i+1);
                            % Find out the end time of the second cycle.
                            end_cycle_time = new_time_value(i+1)+new_cycle_time;
                            x=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                            z=0+((2*pi)/length(x)):(2*pi)/length(x):(2*pi);
                            y=sin(z);
                            mult = mult_array(occurrence(i+1)):(-((... 
                                mult_array(occurrence(i+1))-0)/((length(x))... 
                                -1))):0;
                            new_y = y.*mult;
%                             X=horzcat(X,x);
                            new_Y = horzcat(new_Y,new_y);
                        end
                        
                    end

                end
            end
            
            if(div)
                new_Y = 0.2*(new_Y/max(abs(new_Y)));
            else
                div = true;
            end
            
            % Find out the length of the produced signal and the length of the
            % original signal is same or not.
            compensate_zero=(length(new_Y)-length_sig);
            if (compensate_zero==0)
                % This is not necessary to see everytime.
                % disp('true')
            else
                disp('False. The signal length is not same.')
                disp('channel')
                disp(ch)
                disp('The missing number of values are .....')
                disp(compensate_zero)
                % zero_array_x = 0:(1/Fs):new_time_value(i);
                zero_array_Y = zeros(1,compensate_zero);
                % X=horzcat(zero_array_x,x);
                % Y=horzcat(zero_array_Y,y);
                new_Y=horzcat(zero_array_Y,new_Y);
                disp(' ')
                disp('-----------------------------------------------------')
                disp(' ')
            end
        else
            disp('no onset at all. Investigate .. ')
            disp ('Channel')
            disp(ch)
            new_Y=zeros(1,length_sig);
        end
        
        
    end
    
    % I do not need to normalise here. It will be done at the end. 
    % Now the signals for all the channels have been generated. It is the
    % time to normalise all of them. 
    % new_Y = 0.9*(new_Y/max(abs(new_Y)));
    
    % This is to see the spectrogram of individual signal.
    %     Spectrogram_Compare(new_Y,'signal',Fs, ch);

    % This is to try the power spectrum of each individual channel.
    %     figure(ch)
    %     Hs=spectrum.welch;
    %     psd(Hs,new_Y,'Fs',Fs)
    
%     % if (ch == 43) 
%     % This is to plot the generated signal for this channel with the spikes
%     % This is optional. 
%     xAxis = 0:(1/Fs):((length_sig/Fs)-(1/Fs));
%     disp(' ')
%     figure
%     if (isempty(new_time_value)) 
%         Onset_yAxis = zeros(1,length(xAxis));
%     else
%         disp(['The number of spikes at ' num2str(cochCFs(ch)) 'Hz frequency  (' num2str(ch) 'th Channel) is ' num2str(length(new_time_value))])
%         for i = 1:length(new_time_value)
%             if (i==1)
%                 xAxis_temp = 0:(1/Fs):new_time_value(i);
%                 yAxis_temp = zeros(1,length(xAxis_temp)-1);
%                 yAxis = horzcat(yAxis_temp, occurrence(i));
%             else
%                 xAxis_temp = new_time_value(i-1)+(1/Fs):(1/Fs):new_time_value(i);
%                 yAxis_temp = zeros(1,length(xAxis_temp)-1);
%                 yAxis = horzcat(yAxis,yAxis_temp, occurrence(i));
%             end
%         end
%         xAxis_temp = new_time_value(i)+(1/Fs):(1/Fs):(length_sig/Fs);
%         yAxis_temp = zeros(1,length(xAxis_temp)-1);
%         Onset_yAxis = horzcat(yAxis,yAxis_temp);
%     end
% 
%     plot(xAxis((1:min(length(xAxis), length(Onset_yAxis)))), Onset_yAxis((1:min(length(xAxis), length(Onset_yAxis)))), 'o')
%     % plot the sound signal as well.
%     hold on
%     plot(xAxis((1:min(length(xAxis), length(Onset_yAxis)))), 15*new_Y((1:min(length(xAxis), length(Onset_yAxis)))), 'r')
%     grid on
%         %ylim([0 18])
%     title(['Channel Frequency: - ' num2str(cochCFs(ch)) ' Hz'])
%      % end
    
    
    % At last, add the generated signal to the cell array
    multiplied_zc_cell{1,ch} = new_Y;
    waitbar(ch /num_channels)
    
    
    
%     % This is to Plot the reconstructed signal for each channel
%     figure
%     t = 0+(1/Fs):(1/Fs):length_sig/Fs;
%     data = audioread([an_data_files(1:end-13) '.wav']); % change the wav to something else if necessary
%     plot(t,data, 'y')
%     hold on
%     plot(t, new_Y)
%     hold on
%     % Now plot the spikes
% %     j=1;
% %     spike_plot = zeros(1,length(t));
% %     for i = 1:length(t)
% %         if t(i) == time_value(j,1)
% %             spike_plot(i)=time_value(j,2);
% %             j = j+1;
% %         else
% %             spike_plot(i)=0;
% %         end
% %         
% %     end
%     h = stem(time_value(:,1), (time_value(:,2)*0.05625), 'fill','green');
%     set(h,'MarkerFaceColor','green')
%     % grid on
    
    
    
    
end
close(h3)
   
disp('Generating Signals are done.')
disp(['Total nuber of Spikes used in this Decoding is: ' num2str(TotalNoSpike)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end
