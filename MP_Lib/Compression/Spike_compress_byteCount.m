function [ total_KilloBytes ] = Spike_compress_byteCount( data_files, type, spike_rate)
%Spike_compress_byteCount: 
%   This function calculates how many bytes are required to store all the
%   spikes in a sound file. 
%   The 'data_files' is the name of the data file. 
%   The 'type' is the AN spike or Onset spike. 
%   'spike_rate' is for applying maximum spiking rate or not. 
%   The example call of this function - 
%       'Spike_compress_byteCount('duduk_16_50_AN.mat', 1, 0);'

% %disp(an_data_files);
% % How many files?
% num_files       = length(an_data_files);
% %disp(num_files);
% 
% % Now loop over the number of files, loading them into structures one by one
% for i=1:num_files
%     big(i).data = load(char(an_data_files(i)));        
% end

big.data = load(data_files); 

% Now the actual spike rate is 1/NoofSpike per second 
max_spike_rate = 200;
max_spike_rate = (1/max_spike_rate);


if type == 1
    % % Pick out some useful information
    % % How many channels?
    num_channels    = big(1).data.AN.channels;
    cochCFs = big(1).data.AN.cochCFs;
    % % How many sen levels?
    % num_sen_levels  = big(1).data.AN.iterations;
    length_sig = big(1).data.AN.datalength;
    fs = big(1).data.AN.fs;
    AN_data = big(1).data.AN.Spike_Assign_Channels;
    output_str = 'AN';
elseif type == 2
    % How many channels for Onset spikes?
    num_channels    = big(1).data.ANparams.channels;
    cochCFs = big(1).data.ANparams.cochCFs;
    % How many sen levels for Onset spikes?
    num_Onset_sen_levels  = big(1).data.ANparams.iterations;
    length_sig = big(1).data.ANparams.datalength;
    fs = big(1).data.ANparams.fs;
    AN_data = big(1).data.onsetparams.OnsetSpike_Channel;
    output_str = 'Onset';
else
    disp('The type has to be either 1 or 2')
end

sound_length = (length_sig/fs);

% Now construct the final time value arrary 
% The maximum length of this array will be the number of spikes in the
% first sensitivity level. 
% 
% final_time_values = zeros(1, length(AN_data(1).list));
% 
% for i = 1:num_sen_levels-1
%     [final_time_values, iA, iB] = intersect(AN_data(i).list, AN_data(i+1).list);
% end
% 


freq_mult = 10;
% Take the spike occurence data from the big structure
% 
total_bits = 0;

ch_array = zeros(1, num_channels); 
spike_length_array = zeros(1, num_channels); 
byte_array = zeros(1, num_channels); 
TotalNoSpike = 0;

for ch = 1:num_channels
    % What is the center frequency of this chennel?
    channel_frequency = cochCFs(ch);
    time_value = AN_data{1,ch};
    
    if time_value ~= 0
        
        new_time_value = time_value(:,1);
        occurrence = time_value(:,2);


        if spike_rate == 1
            %% This is to implement the maximum spiking rate, i.e. 200. 
            % To avoid, just comment these bit out. 
            % So, each sine wave cycle will take (1/frequency) sec.
            cycle_time=1/channel_frequency;
            NoOfSpike=length(AN_data{1,ch});
            Sorted_new_time_value=0;    % The sorted spike time values.
            Sorted_new_occurence=0;     % The sorted multipliers.
            if (cycle_time<max_spike_rate)
                % Do the spiking rate work.
                Sorted_new_time_value=0;    % The sorted spike time values.
                Sorted_new_occurence=0;     % The sorted multipliers.
                j=1;
                for i = 1:NoOfSpike
                    if (i==1)
                        Sorted_new_time_value(j)=new_time_value(i);
                        Sorted_new_occurence(j)=occurrence(i);
                        j=j+1;
                        next_counter = new_time_value(1);
                    end
                    if (new_time_value(i)>=(next_counter+max_spike_rate))
                        Sorted_new_time_value(j)=new_time_value(i);
                        Sorted_new_occurence(j)=occurrence(i);
                        next_counter = Sorted_new_time_value(j);
                        j=j+1;
                    end
                end
                spike_length = length(Sorted_new_time_value);
            else
                spike_length = length(AN_data{1,ch});
            end
            % The end of implementing spiking rate. 
        elseif spike_rate == 0
            spike_length = length(AN_data{1,ch});
        else
            disp('The value for "spike_rate" has to be either 1 or 0')
        end

        % Now to avoid error we multiply the frequency by the freq_multiplier
        channel_frequency = channel_frequency*freq_mult;
        expo = ceil(log10(channel_frequency)); % This is number of zeros, the multiplier should have. 
        mult = 10^expo; % This is the level of accuracy, depends on the channel number. So, for lower 
            % channel, the accuray is low, but for the higher channel, the
            % accuracy has to be high as there are more spikes. 
        
        new_sound_length = sound_length*mult; % This is the adjusted sound length with accuracy. 
        bit_length = ceil(log2(new_sound_length)); % log2 of the adjusted sound length will the required bits
        % The total number of bits required is - 
        total_bits = total_bits + (bit_length*spike_length);

        
        ch_array(ch) = ch;
        spike_length_array(ch) = spike_length;
        byte_array(ch) = ((bit_length*spike_length)/(8*1024)); % 1 Byte = 8 bits; 1 Killobyte = 1024 Bytes
        
    else
        bit_length=0;
        spike_length = 0;
        ch_array(ch) = ch;
        spike_length_array(ch) = 0;
        byte_array(ch) = 0;
        
    end
    
    % To display the results for each channel
    disp('')
    disp(['Spikes are - : ' num2str(spike_length) ' at channel: ' num2str(ch) ' of - ' num2str((bit_length*spike_length)/(8*1024)) ' Kilobytes'])
    TotalNoSpike = TotalNoSpike + spike_length; 
end

total_KilloBytes = total_bits/(8*1024);

disp(' ')
disp(' ------------------------------------------------ ')
disp(['The total number of spikes: - ' num2str(TotalNoSpike)])
disp(['Bytes required for ' output_str ' spikes of a ' num2str(sound_length) ' sec sound file - '])
disp([num2str(total_KilloBytes), ' Kilobytes'])
disp(' ------------------------------------------------ ')
disp(' ')
disp(' ')

% figure
% plot(ch_array, spike_length_array)
% figure
% plot(ch_array, byte_array)

figure
 % Now draw the graphs for channels and the spike numers and their bytes
[haxes,hline1,hline2] = plotyy(ch_array, spike_length_array, ch_array, byte_array,'semilogy', 'plot');
ylabel(haxes(1),'Number of Spikes') % label left y-axis
ylabel(haxes(2),'Required KiloBytes') % label right y-axis
xlabel(haxes(2),'Channels') % label x-axis
set(hline1,'LineStyle','--','LineWidth',2);
set(hline2,'LineWidth',2);
grid(haxes(1),'on');
grid(haxes(2),'on');

end

