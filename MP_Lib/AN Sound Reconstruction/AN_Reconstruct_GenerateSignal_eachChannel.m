%function [Final_Signal,multiplied_zc_cell,new_zc_cell,zc_cell] = AN_Reconstruct(an_data_files)
function [multiplied_zc_cell, TotalNoSpikes, TotalSortedNoSpikes,new_channel_frequency_array] = AN_Reconstruct_GenerateSignal_eachChannel(an_data_files,delay_zc_cell,max_spike_rate,Jitrsp_tech)
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

% How many spike are altogether present in this spike based coding system?
% The answer is necessary to estimate the level of loosen spikes by
% Maximum spiking rate.
TotalNoSpikes = 0;

% How many spikes are being used for this sound to reconstruct?
TotalSortedNoSpikes = 0; 

% Now the actual spike rate is 1/NoofSpike per second 
max_spike_rate = (1/max_spike_rate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the part, which actually generates the signals for each channel.

multiplied_zc_cell = cell(1,num_channels);
% Calculating the multiplier from the threshhold levels
temp_mult = 1/thresh_levels(length(thresh_levels));
mult_array = zeros(1,length(thresh_levels));
for i = 1:length(thresh_levels)
    mult_array(i) = temp_mult*thresh_levels(i);
end
h3 = waitbar(0,'Please wait untill the sine curves are produced...');
for ch = 1:num_channels
    %time_value = new_zc_cell{1,ch};
    time_value = delay_zc_cell{1,ch};
    new_time_value = time_value(:,1);
    occurrence = time_value(:,2);
    % Get the center-frequency of this channel.
    channel_frequency = cochCFs(ch);
    % So, each sine wave cycle will take (1/frequency) sec.
    cycle_time=1/channel_frequency;
    NoOfSpike=length(new_time_value);
%     disp('Channel')
%     disp(ch)
    
    % Now the legth of the 'new_time_value' is necessary to consider. 
    % If that is empty, just add some zeros....
    if (isempty(new_time_value))
        X = 0:(1/Fs):((length_sig)/Fs);
        new_Y = zeros(1,length(X));
    % If there is only one spike, generate two cycles same as the center frequency of this channel.
    elseif (length(new_time_value)==1)
        % Find out the starting time of the first cycle.
        start_cycle_time = new_time_value-cycle_time;
        % Find out the end time of the second cycle.
        end_cycle_time = new_time_value+cycle_time;
        
        % This is for the first cycle
        x_1=start_cycle_time+(1/Fs):(1/Fs):new_time_value;
        z_1=0+((2*pi)/length(x_1)):(2*pi)/length(x_1):(2*pi);
        y_1=sin(z_1);
        mult_1 = 0:((mult_array(occurrence)-0)/((length(x_1))-1)):mult_array(occurrence);
        new_y_1 = y_1.*mult_1;
        
        % This is for the second cycle
        x_2=new_time_value+(1/Fs):(1/Fs):end_cycle_time;
        z_2=0+((2*pi)/length(x_2)):(2*pi)/length(x_2):(2*pi);
        y_2=sin(z_2);
        mult_2 = mult_array(occurrence):(-((mult_array(occurrence)-0)/((length(x_2))-1))):0;
        new_y_2 = y_2.*mult_2;
        
        % Now add some zeros at the Beggining.
        zero_array_x_1 = 0:(1/Fs):start_cycle_time;
        zero_array_Y_1 = zeros(1,length(zero_array_x_1));
        % Now add some zeros at the End.
        zero_array_x_2 = end_cycle_time+(1/Fs):(1/Fs):((length_sig)/Fs);
        zero_array_Y_2 = zeros(1,length(zero_array_x_2));
        
        % Now concate everything in this channel.
        X=horzcat(zero_array_x_1,x_1,x_2,zero_array_x_2);
        new_Y=horzcat(zero_array_Y_1,new_y_1,new_y_2,zero_array_Y_2); % 'new_Y' is the final signal for this channel
        % plot(X,new_Y)
    % If there is more than one spike, go for the normal procedure....
    % So for this channel there are plenty spikes. Go Ahead...
    else
        % First check that this channel has the spike rates less than 200Hz
        % THE REASON:- Our Auditory nerve can generate spikes in a certain
        % rate. The maximum spiking rate is 200 per second. There are too
        % many spikes, which can not be produced inside auditory nerve, at 
        % close interval at the higher frequency channels. 
        % It is expected that this lossy technique will not have a huge
        % impact on the reconstructed signal, as our ear is not very
        % sensitive for very high frequencies. 
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
            
            % disp('the sorting work done')
            % Now do some graphical ploting to see the sorted and original values.
%             figure(ch)
%             plot_x_axix = zeros(1,length(new_time_value));
%             plot(new_time_value,plot_x_axix,'x')
%             plot_y_axix = zeros(1,length(new_time_value));
%             for i=1:length(plot_x_axix)
%                 for j=1:length(Sorted_new_time_value)
%                     if (new_time_value(i)==Sorted_new_time_value(j))
%                         plot_y_axix(i)=new_time_value(i);
%                     end
%                 end
%             end
%             hold on
%             plot(plot_y_axix,plot_x_axix,'or')
%             grid on
            
%             % This is to show how many spikes were in this channel and how
%             % many spikes are now in this channel after sorting.
            NoOfSortedSpike = length(Sorted_new_time_value);
%             disp('Channel ')
%             disp(ch)
%             disp('Channel Frequency ')
%             disp(channel_frequency)
%             disp('The Number of the spikes were - ')
%             disp(NoOfSpike)
%             disp('The Number of the spikes are NOW - ')
%             disp(NoOfSortedSpike)
            TotalNoSpikes = TotalNoSpikes+NoOfSpike;
            TotalSortedNoSpikes = TotalSortedNoSpikes+NoOfSortedSpike;
%             disp('Total Number of Spike ')
%             disp(TotalNoSpikes)
%             disp('Total Number of Sorted Spike ')
%             disp(TotalSortedNoSpikes)
%             disp(' ')
%             disp(' ')
            
            % Some extra code to check that the sorting is correct.
            for i=1:(length(Sorted_new_time_value)-1)
                inner_diff = Sorted_new_time_value(i+1)-Sorted_new_time_value(i);
                if inner_diff<(1/200)
                    disp('The sorting technique is not right.')
                    disp('The difference between two spikes are : ')
                    disp(inner_diff)
                end
            end
            for i = 1:length(Sorted_new_time_value)-1
                % What is the differnce between this spike and next spike?
                temp_length = Sorted_new_time_value(i+1)-Sorted_new_time_value(i);
                x=Sorted_new_time_value(i)+(1/Fs):(1/Fs):Sorted_new_time_value(i+1);
                % Count how many cycles will be in between these two sorted spikes.
                temp_cycles = temp_length*channel_frequency;
                % Now this number of cycles should be a integer
                temp_cycles = floor(temp_cycles);
                % Now the frequency of the signal is important in this case.
                if ch == 2
                    if i == 9
                        disp(' ')
                    end
                end
                % For 0, there is no change in frequency - 
                if Jitrsp_tech == 0
                    new_channel_frequency = temp_cycles/temp_length;
                    z = 0+(1/Fs):(temp_length-(1/Fs))/(length(x)-1):temp_length;
                    y = sin(2*pi*z*new_channel_frequency);
                elseif Jitrsp_tech == 1
                    % To reduce the frequency effect of Spike Cap technique,
                    % the frequency values has been adjusted. The frequency
                    % range is - 
                    if ch == 1
                        channel_frequency_range = (cochCFs(2)-cochCFs(1))/2;
                    else
                        channel_frequency_range = (cochCFs(ch)-cochCFs(ch-1))/2;
                    end
                    new_channel_frequency = random('Normal',channel_frequency,channel_frequency_range,1,1);
                    % Count how many cycles will be in between these two sorted spikes.
                    temp_cycles = temp_length*new_channel_frequency;
                    % Now this number of cycles sholud be a integer
                    temp_cycles = floor(temp_cycles);
                    new_channel_frequency = temp_cycles/temp_length;
                    z = 0+(1/Fs):(temp_length-(1/Fs))/(length(x)-1):temp_length;
                    y = sin(2*pi*z*new_channel_frequency);
                elseif Jitrsp_tech == 2
                    start_time = Sorted_new_time_value(i);
                    this_temp_length = temp_length/temp_cycles;
                    for j = 1:temp_cycles
                        if j == temp_cycles
                            end_time = start_time+this_temp_length;
                        else
                            % The idea of jittering cycle time is applied here.
                            % The cycle end time is randomly changed by this prcentage - 
                            percent_change = 10;
                            jittered_length = random('Normal',0,((1/channel_frequency)*(percent_change/100)),1,1);
                            end_time = start_time+this_temp_length+jittered_length;
                        end
                        x1=start_time+(1/Fs):(1/Fs):end_time;
                        if ch == 1
                            channel_frequency_range = (cochCFs(2)-cochCFs(1))/100;
                        else
                            channel_frequency_range = (cochCFs(ch)-cochCFs(ch-1))/100;
                        end
                        %channel_frequency_range = 0;
                        new_channel_frequency = random('Normal',channel_frequency,channel_frequency_range,1,1);
                        % Count how many cycles will be in between these two sorted spikes.
                        temp_cycles = this_temp_length*new_channel_frequency;
                        % Now this number of cycles should be a integer
                        temp_cycles = floor(temp_cycles);
                        if temp_cycles == 0
                            temp_cycles = 1;
                        end
                        new_channel_frequency = temp_cycles/this_temp_length;
                        z1 = 0+(1/Fs):(this_temp_length-(1/Fs))/(length(x1)-1):this_temp_length;
                        y1 = sin(2*pi*z1*new_channel_frequency);
                        if j==1
                            y = y1;
                        else
                            y = horzcat(y,y1);
                        end
                        start_time = end_time;
                    end
                else
                    disp('The values for Jittered Spike technique must be 0, 1 or 2.')
                end
                % This work to find out the maximum channel frequency.
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
                if mult_array(Sorted_new_occurence(i)) == mult_array(Sorted_new_occurence(i+1))
                    mult = mult_array(Sorted_new_occurence(i+1));
                    new_y = y*mult;
                else
                    if mult_array(Sorted_new_occurence(i)) < mult_array(Sorted_new_occurence(i+1))
                        mult = mult_array(Sorted_new_occurence(i)):((mult_array(Sorted_new_occurence(i+1))-mult_array(Sorted_new_occurence(i)))/((length(x))-1)):mult_array(Sorted_new_occurence(i+1));                          
                    else
                        mult = mult_array(Sorted_new_occurence(i)):(-((mult_array(Sorted_new_occurence(i))-mult_array(Sorted_new_occurence(i+1)))/((length(x))-1))):mult_array(Sorted_new_occurence(i+1));
                    end
                    % This is to check the length of signal and multipliers
                    if (length(y) == length(mult))
                        % Do nothing. 
                    else
                        if (length(y)-length(mult))>0
                            zerosToAdd = zeros(1,length(y)-length(mult));
                            mult = horzcat(mult,zerosToAdd);
                        else
                            zerosToAdd = zeros(1,length(mult)-length(y));
                            y = horzcat(y,zerosToAdd);
                        end
                        
                    end
                    new_y = y.*mult;
                end
                
                % It is necessary to add some signal at the Beginning.
                if (i==1)
                    % It is necessary to check that there is enough room to add
                    % another sine wave. If not then, just add zeros...
                    if ((Sorted_new_time_value(i)-cycle_time)<0)
                        disp('This time here is no extra sine wave.')
                        zero_array_x = 0:(1/Fs):Sorted_new_time_value(i);
                        zero_array_Y = zeros(1,length(zero_array_x));
                        X=horzcat(zero_array_x,x);
                        Y=horzcat(zero_array_Y,y);
                        new_Y=horzcat(zero_array_Y,new_y); % 'zero_array_Y' is the extra zeros..
                    else
                        % One extra sine wave should be added at the beggining.
                        temp_x = (Sorted_new_time_value(i)-cycle_time)+(1/Fs):(1/Fs):Sorted_new_time_value(i);
        %                 disp('This length has been added for EXTRA SINE WAVE')
        %                 disp(length(temp_x))
                        temp_z=0+((2*pi)/length(temp_x)):(2*pi)/length(temp_x):(2*pi);
                        temp_y = sin(temp_z);
                        mult = 0:(mult_array(Sorted_new_occurence(i))/((length(temp_x))-1)):mult_array(Sorted_new_occurence(i));
                        temp_new_y = temp_y.*mult;
                        X=horzcat(temp_x,x);
                        Y=horzcat(temp_y,y);
                        new_Y=horzcat(temp_new_y,new_y); % 'temp_new_y' is the extra wave

                        % It is necessary to add some zeros at the beginning of each channel signal.
                        zero_array_x = 0:(1/Fs):(Sorted_new_time_value(i)-cycle_time);
        %                 disp('This length of ZEROS has been added before EXTRA SINE WAVE')
        %                 disp(length(zero_array_x))
                        zero_array_y = zeros(1,length(zero_array_x));
                        %zero_array_Y = zero_array_y.*mult;
                        zero_array_Y = zeros(1,length(zero_array_x));
                        X=horzcat(zero_array_x,X);
                        Y=horzcat(zero_array_y,Y);
                        new_Y = horzcat(zero_array_Y,new_Y);
                    end
                % For the end cycle some work to do.........
                elseif (i==(length(Sorted_new_time_value)-1))
                    % Calculate how much space is available at the End
                    end_space = ((length_sig)/Fs)-Sorted_new_time_value(i+1);
                    % If that spave is big enough, add the sine cycle at the channel frequency
                    if (end_space>=cycle_time)
                        % Find out the starting time of the first cycle.
                        start_cycle_time = Sorted_new_time_value(i+1);
                        % Find out the end time of the second cycle.
                        end_cycle_time = Sorted_new_time_value(i+1)+cycle_time;
                        x=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                        z=0+((2*pi)/length(x)):(2*pi)/length(x):(2*pi);
                        y=sin(z);
                        mult = mult_array(Sorted_new_occurence(i+1)):(-((mult_array(Sorted_new_occurence(i+1))-0)/((length(x))-1))):0;
                        new_y = y.*mult;
                        X=horzcat(X,x);
                        new_Y = horzcat(new_Y,new_y);
                        % Just add zeros upto end of the signal.
                        zero_array_x = length_sig-length(new_Y);
                        % zero_array_x = end_cycle_time+(1/Fs):(1/Fs):((length_sig)/Fs);
                        zero_array_y = zeros(1,zero_array_x);
    %                     disp('This length of ZEROS has been added AT THE END')
    %                     disp(length(zero_array_x))
    %                     disp(length(zero_array_Y))
                        new_Y = horzcat(new_Y,zero_array_y);
                    % If not, add the sine wave at the length equal to the space left at the end.
                    else
                        disp('Hey, This is a rare case. Investigate through it.')
                        disp('Channel')
                        disp(ch)
                        new_cycle_time = end_space;
                        % Find out the starting time of the first cycle.
                        start_cycle_time = Sorted_new_time_value(i+1);
                        % Find out the end time of the second cycle.
                        end_cycle_time = Sorted_new_time_value(i+1)+new_cycle_time;
                        x=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                        z=0+((2*pi)/length(x)):(2*pi)/length(x):(2*pi);
                        y=sin(z);
                        mult = mult_array(Sorted_new_occurence(i+1)):(-((mult_array(Sorted_new_occurence(i+1))-0)/((length(x))-1))):0;
                        new_y = y.*mult;
                        X=horzcat(X,x);
                        new_Y = horzcat(new_Y,new_y);
                    end
                % If it is not at the first or last case, do the concate, but
                % remember to check and do stuffs for each big difference between spikes.     
                else
                    spike_time_diff = temp_length;
                    % What is a reasonable gap?
                    % It shold be the normal cyle time of this chennel frequency
                    % spike_gap = cycle_time;
                    % Check if there is a big gap or not.
                    % Do extra things if the diffrence is 2 times greater
                    % than the spike time diffrence.
                    if (spike_time_diff>(max_spike_rate*2))
                        % diff = spike_time_diff-(spike_gap);
                        % diff_a = 0:(1/Fs):spike_time_diff;
                        % Generate the first cycle at the beggining
                        % Find out the starting time of the first cycle.
                        start_cycle_time = Sorted_new_time_value(i);
                        % Find out the end time of the first cycle.
                        end_cycle_time = Sorted_new_time_value(i)+cycle_time;
                        % This is for the first cycle
                        x_1=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                        z_1=0+((2*pi)/length(x_1)):(2*pi)/length(x_1):(2*pi);
                        y_1=sin(z_1);
                        mult_1 = mult_array(Sorted_new_occurence(i)):(-((mult_array(Sorted_new_occurence(i))-0)/((length(x_1))-1))):0;
                        new_y_1 = y_1.*mult_1;

                        % Generate some zeros at the middle if there is enough space
                        start_zero_time=end_cycle_time;
                        end_zero_time = Sorted_new_time_value(i+1)-cycle_time;
                        x_2=start_zero_time+(1/Fs):(1/Fs):end_zero_time;
                        new_y_2=zeros(1,length(x_2));

                        % Generate second cycle at the end
                        % Find out the starting time of the second cycle.
                        start_cycle_time = Sorted_new_time_value(i+1)-cycle_time;
                        % Find out the end time of the second cycle.
                        end_cycle_time = Sorted_new_time_value(i+1);
                        % This is for the first cycle
                        x_3=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                        z_3=0+((2*pi)/length(x_3)):(2*pi)/length(x_3):(2*pi);
                        y_3=sin(z_3);
                        mult_3 = 0:((mult_array(Sorted_new_occurence(i+1))-0)/((length(x_3))-1)):mult_array(Sorted_new_occurence(i+1));
                        new_y_3 = y_3.*mult_3;

                        new_Y=horzcat(new_Y,new_y_1,new_y_2,new_y_3);
                    else
                        X=horzcat(X,x);
                        Y=horzcat(Y,y);
                        new_Y=horzcat(new_Y,new_y);
                    end
                end
            end
        else
%             % This is to show how many spikes were in this channel and how
%             % many spikes are now in this channel after sorting.
%             % NoOfSpike=length(new_time_value);
%             % NoOfSortedSpike = length(Sorted_new_time_value);
%             disp('Channel ')
%             disp(ch)
%             disp('Channel Frequency ')
%             disp(channel_frequency)
%             disp('The Number of the spikes are - ')
%             disp(NoOfSpike)
%             % disp('The Number of the spikes are NOW - ')
%             % disp(NoOfSortedSpike)
            TotalNoSpikes = TotalNoSpikes+NoOfSpike;
            TotalSortedNoSpikes = TotalSortedNoSpikes+NoOfSpike;
%             disp('Total Number of Spike ')
%             disp(TotalNoSpikes)
%             disp('Total Number of Sorted Spike ')
%             disp(TotalSortedNoSpikes)
%             disp(' ')
%             disp(' ')
            for i = 1:length(new_time_value)-1
                % What is the differnce between this spike and next spike?
                temp_length = new_time_value(i+1)-new_time_value(i);
                
                % This work is to see the frequency of this spike interval.
                % So, the frequency of the signal is necessary in this case is
                new_channel_frequency = 1/temp_length;
                % This work to find out the maximum channel frequency.
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
                
                x=new_time_value(i)+(1/Fs):(1/Fs):new_time_value(i+1);
                z=0+((2*pi)/length(x)):(2*pi)/length(x):(2*pi);
                y=sin(z);
                if mult_array(occurrence(i)) == mult_array(occurrence(i+1))
                    mult = mult_array(occurrence(i+1));
                    new_y = y*mult;
                else
                    if mult_array(occurrence(i)) < mult_array(occurrence(i+1))
                        mult = mult_array(occurrence(i)):((mult_array(occurrence(i+1))-mult_array(occurrence(i)))/((length(x))-1)):mult_array(occurrence(i+1));
                        new_y = y.*mult;
                    else
                        mult = mult_array(occurrence(i)):(-((mult_array(occurrence(i))-mult_array(occurrence(i+1)))/((length(x))-1))):mult_array(occurrence(i+1));
                        new_y = y.*mult;
                    end
                end
                % It is necessary to add some signal at the Beginning.
                if (i==1)
                    % It is necessary to check that there is enough room to add
                    % another sine wave. If not then, just add zeros...
                    if ((new_time_value(i)-temp_length)<0)
        %                 disp('This time here is no extra sine wave.')
                        zero_array_x = 0:(1/Fs):new_time_value(i);
                        zero_array_Y = zeros(1,length(zero_array_x));
                        X=horzcat(zero_array_x,x);
                        Y=horzcat(zero_array_Y,y);
                        new_Y=horzcat(zero_array_Y,new_y); % 'zero_array_Y' is the extra zeros..
                    else
                        % One extra sine wave should be added at the beggining.
                        temp_x = (new_time_value(i)-temp_length)+(1/Fs):(1/Fs):new_time_value(i);
        %                 disp('This length has been added for EXTRA SINE WAVE')
        %                 disp(length(temp_x))
                        temp_z=0+((2*pi)/length(temp_x)):(2*pi)/length(temp_x):(2*pi);
                        temp_y = sin(temp_z);
                        mult = 0:(mult_array(occurrence(i))/((length(temp_x))-1)):mult_array(occurrence(i));
                        temp_new_y = temp_y.*mult;
                        X=horzcat(temp_x,x);
                        Y=horzcat(temp_y,y);
                        new_Y=horzcat(temp_new_y,new_y); % 'temp_new_y' is the extra wave

                        % It is necessary to add some zeros at the beginning of each channel signal.
                        zero_array_x = 0:(1/Fs):(new_time_value(i)-temp_length);
        %                 disp('This length of ZEROS has been added before EXTRA SINE WAVE')
        %                 disp(length(zero_array_x))
                        zero_array_y = zeros(1,length(zero_array_x));
                        %zero_array_Y = zero_array_y.*mult;
                        zero_array_Y = zeros(1,length(zero_array_x));
                        X=horzcat(zero_array_x,X);
                        Y=horzcat(zero_array_y,Y);
                        new_Y = horzcat(zero_array_Y,new_Y);
                    end
                % For the end cycle some work to do.........
                elseif (i==(length(new_time_value)-1))
                    % Calculate how much space is available at the End
                    end_space = ((length_sig)/Fs)-new_time_value(i+1);
                    % If that spave is big enough, add the sine cycle at the channel frequency
                    if (end_space>=cycle_time)
                        % Find out the starting time of the first cycle.
                        start_cycle_time = new_time_value(i+1);
                        % Find out the end time of the second cycle.
                        end_cycle_time = new_time_value(i+1)+cycle_time;
                        x=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                        z=0+((2*pi)/length(x)):(2*pi)/length(x):(2*pi);
                        y=sin(z);
                        mult = mult_array(occurrence(i+1)):(-((mult_array(occurrence(i+1))-0)/((length(x))-1))):0;
                        new_y = y.*mult;
                        X=horzcat(X,x);
                        new_Y = horzcat(new_Y,new_y);
                        % Just add zeros upto end of the signal.
                        zero_array_x = length_sig-length(new_Y);
                        % zero_array_x = end_cycle_time+(1/Fs):(1/Fs):((length_sig)/Fs);
                        zero_array_y = zeros(1,zero_array_x);
    %                     disp('This length of ZEROS has been added AT THE END')
    %                     disp(length(zero_array_x))
    %                     disp(length(zero_array_Y))
                        new_Y = horzcat(new_Y,zero_array_y);
                    % If not, add the sine wave at the length equal to the space left at the end.
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
                        mult = mult_array(occurrence(i+1)):(-((mult_array(occurrence(i+1))-0)/((length(x))-1))):0;
                        new_y = y.*mult;
                        X=horzcat(X,x);
                        new_Y = horzcat(new_Y,new_y);
                    end
                % If it is not at the first or last case, do the concate, but
                % remember to check and do stuffs for each big difference between spikes. 
                else
                    spike_time_diff = temp_length;
                    % What is a reasonable gap?
                    % It shold be the normal cyle time of this chennel frequency
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
                        mult_1 = mult_array(occurrence(i)):(-((mult_array(occurrence(i))-0)/((length(x_1))-1))):0;
                        new_y_1 = y_1.*mult_1;

                        % Generate some zeros at the middle if there is enough space
                        start_zero_time=end_cycle_time;
                        end_zero_time = new_time_value(i+1)-cycle_time;
                        x_2=start_zero_time+(1/Fs):(1/Fs):end_zero_time;
                        new_y_2=zeros(1,length(x_2));

                        % Generate second cycle at the end
                        % Find out the starting time of the second cycle.
                        start_cycle_time = new_time_value(i+1)-cycle_time;
                        % Find out the end time of the second cycle.
                        end_cycle_time = new_time_value(i+1);
                        % This is for the first cycle
                        x_3=start_cycle_time+(1/Fs):(1/Fs):end_cycle_time;
                        z_3=0+((2*pi)/length(x_3)):(2*pi)/length(x_3):(2*pi);
                        y_3=sin(z_3);
                        mult_3 = 0:((mult_array(occurrence(i+1))-0)/((length(x_3))-1)):mult_array(occurrence(i+1));
                        new_y_3 = y_3.*mult_3;

                        new_Y=horzcat(new_Y,new_y_1,new_y_2,new_y_3);
                    else
                        X=horzcat(X,x);
                        Y=horzcat(Y,y);
                        new_Y=horzcat(new_Y,new_y);
                    end
                end
            end
        end
    end
    
    % Find out the length of the produced signal and the length of the
    % original signal is same or not.
    compensate_zero=abs(length(new_Y)-length_sig);
    if (compensate_zero==0)
        % This is not necessary to see everytime.
        % disp('true')
    else
        disp('False. The signal length is not same.')
        disp('The missing number of values are .....')
        disp('channel')
        disp(ch)
        disp(compensate_zero)
        % zero_array_x = 0:(1/Fs):new_time_value(i);
        zero_array_Y = zeros(1,compensate_zero);
        % X=horzcat(zero_array_x,x);
        % Y=horzcat(zero_array_Y,y);
        new_Y=horzcat(zero_array_Y,new_Y);
        disp(' ')
        disp('-----------------------------------------------------------------')
        disp(' ')
    end
    
    % This is to see the spectrogram of individual signal.
%     Spectrogram_Compare(new_Y,'signal',Fs, ch);
    
    % This is to try the power spectrum of each individual channel.
%     figure(ch)
%     Hs=spectrum.welch;
%     psd(Hs,new_Y,'Fs',Fs)
    
    multiplied_zc_cell{1,ch} = new_Y;
    waitbar(ch /num_channels)
end
close(h3)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%