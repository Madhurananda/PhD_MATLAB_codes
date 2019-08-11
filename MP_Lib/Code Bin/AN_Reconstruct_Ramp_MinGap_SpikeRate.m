%function [Final_Signal,multiplied_zc_cell,new_zc_cell,zc_cell] = AN_Reconstruct(an_data_files)
function [OutPut] = AN_Reconstruct_Ramp_MinGap_SpikeRate(an_data_files, channel)
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
new_zc_cell = cell(1,num_channels);
log = false;
h2 = waitbar(0,'Please wait. Second cell array is being created ...');
for ch = 1:num_channels
    for i = 1:num_sen_levels

        Time_Value = AN_FindData(an_data_files,(num_sen_levels-i+1),ch);
        %zc_cell{i,ch}=AN_FindData(an_data_files,i,ch);
        %Time_Value = zc_cell{(num_sen_levels-i+1),ch};
        Time_Value(Time_Value==0)=[];
        if (isempty(Time_Value))
            log = true;
            Copy_Time_Value = 0;
        else
            if i==1
                % Do nothing.
%                 if (isempty(Time_Value))
%                     Copy_Time_Value = 0;
%                 end
            else
%                 disp(ch)
%                 disp(i)
                Time_Value = setdiff(Time_Value,Copy_Time_Value);
                
            end
            levels = zeros(length(Time_Value),1);
            for j = 1:length(Time_Value)
                levels(j,1) = (num_sen_levels-i+1);
            end
            if (isempty(Time_Value))
                new_Time_Value = [];
            else
                new_Time_Value = cat(2,Time_Value,levels);
            end
            if i==1
                final_Time_Value = new_Time_Value;
            elseif log
                final_Time_Value = new_Time_Value;
                log = false;
            else
                final_Time_Value = cat(1,final_Time_Value,new_Time_Value);
            end
            % It is necessary to copy the data inside Time_Value at this
            % instance , so that it can be compared to the next Time_Value.
            if (isempty(final_Time_Value))
                
            else
                Copy_Time_Value = final_Time_Value(:,1);
            end
        end
    
    end
    final_Time_Value = sortrows(final_Time_Value,1);
    new_zc_cell{1,ch} = final_Time_Value;
    waitbar(ch/num_channels)
end
close(h2)

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
%     disp('Channel')
%     disp(ch)
    
    % Now the legth of the 'new_time_value' is necessary to consider. 
    % If that is empty, just add some zeros....
    % If there is only one spike, generate two cycles same as the center frequency of this channel.
    % If there is more than one spike, go for the normal procedure....
    if (isempty(new_time_value))
        X = 0:(1/Fs):((length_sig)/Fs);
        new_Y = zeros(1,length(X));
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
    % So for this channel there are plenty spikes. Go Ahead...
    else
        % First check that this channel has the spike rates less than 200Hz
        % THE REASON:- Our Auditory nerve can generate spikes in a certain
        % rate. The maximum spiking rate is 200 per second. There are too
        % many spikes, which can not be produced inside auditory nerve, at 
        % close interval at the higher frequency channels. 
        % It is expected that this lossy technique will not have a huge
        % impact on the reconstructed signal, as our ear is not very
        % sensitive for too high frequencies. 
        max_spike_rate = (1/200);
        if (cycle_time<max_spike_rate)
            % Do the spiking rate work.
%             disp(['Channel' num2str(ch)])
            % It is necessary to sort the data according to need.
            NoOfSpike=length(new_time_value);
            End_SpikeTimeDiff = new_time_value(NoOfSpike)-new_time_value(1);
            % Count how many loops are necessary - 
            loop = (floor(End_SpikeTimeDiff/max_spike_rate)+1);
            % So, what will be the max loop counter value?
            max_loop_value = loop*max_spike_rate;
            % Sorted_new_time_value=zeros(1,loop);
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
            
            % Some extra code to check that the sorting is correct.
            for i=1:(length(Sorted_new_time_value)-1)
                inner_diff = Sorted_new_time_value(i+1)-Sorted_new_time_value(i);
                if inner_diff<=(1/200)
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
                % Now this number of cycles sholud be a integer
                temp_cycles = floor(temp_cycles);
                % So, the frequency of the signal is necessary in this case is
                new_channel_frequency = temp_cycles/temp_length;
                %z=0+((2*pi)/length(x)):(2*pi)/length(x):(2*pi);
                %y=sin(z);
                if ch == 20
                    disp(ch);
                end
                % z = 0:(1/length(x)):1-(1/length(x));
                % z = 0+(1/Fs):(1/Fs):temp_length;
                z = 0+(1/Fs):(temp_length-(1/Fs))/(length(x)-1):temp_length;
%                 if (ch==14)
%                     if (i==22)
%                         disp(i);
%                     end
%                 end
%                 if (length(x)==length(z))
%                     
%                 elseif (length(x)>length(z))
%                     z = 0:(1/Fs):temp_length;
%                 elseif (length(x)<length(z))
%                     z = 0:(1/Fs):temp_length+(1/Fs);
%                 end
                y=sin(2*pi*z*new_channel_frequency);
                % y=sin(2*pi*x*new_channel_frequency);
                if mult_array(Sorted_new_occurence(i)) == mult_array(Sorted_new_occurence(i+1))
                    mult = mult_array(Sorted_new_occurence(i+1));
                    new_y = y*mult;
                else
                    if mult_array(Sorted_new_occurence(i)) < mult_array(Sorted_new_occurence(i+1))
                        mult = mult_array(Sorted_new_occurence(i)):((mult_array(Sorted_new_occurence(i+1))-mult_array(Sorted_new_occurence(i)))/((length(x))-1)):mult_array(Sorted_new_occurence(i+1));
%                         disp(ch)
%                         disp(i)
%                         disp(length(y))
%                         disp(length(mult))
                        new_y = y.*mult;
                    else
                        mult = mult_array(Sorted_new_occurence(i)):(-((mult_array(Sorted_new_occurence(i))-mult_array(Sorted_new_occurence(i+1)))/((length(x))-1))):mult_array(Sorted_new_occurence(i+1));
%                         disp(length(y))
%                         disp(length(mult))
                        new_y = y.*mult;
                    end
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
            for i = 1:length(new_time_value)-1
                % What is the differnce between this spike and next spike?
                temp_length = new_time_value(i+1)-new_time_value(i);
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
                % This is to see the signal for the first channel.
        %         if ch == 1
        %             figure(1)
        %             plot(X,new_Y)
        %         end

            end
            
            % Now replace the original values by the new sorted values.
            % new_time_value = Sorted_new_time_value;
        end
        
        

        
        
        % This is to get the figure for each channel signal.
%         figure(ch)
%         plot(X,new_Y)

        if nargin==2
            if ch == channel
                figure(1)
                plot(X,new_Y)
                hold on

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % This is to plot all the spikes at the given channel and channel.
                m=1;
                for i = 1:length(new_time_value)
                    plot_value(m,1) = new_time_value(i);
                    plot_value(m,2) = 0;
                    m=m+1;
                end
                figure(1)
                plot(plot_value(:,1),plot_value(:,2),'og')
                hold on
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % This is to plot the BMSig (Basilar Membrane Signal)
                % Take the data from the big structure
                temp_an_data_files = [an_data_files(1:end-7) '_bmSig.mat'];
                BM_data =  load(temp_an_data_files);
                BM_data = BM_data.BM.bmSig;
                BM_data = BM_data(ch,:);
                BM_data_axis = 0:(1/Fs):((length_sig/Fs))-(1/Fs);
                figure(1)
                plot(BM_data_axis,BM_data,'.y')

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end

    %     disp('length of original signal')
    %     disp(length_sig)
    %     disp('Length of this signal')
    %     disp(length(new_Y))

        
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
% save('multiplied_zc_cell.mat', 'multiplied_zc_cell')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





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
new_Name = [an_data_files(1:end-7) '_NEW_SpikeCap'];
new_Name = [new_Name '.wav'];
wavwrite(Final_Signal,44100,16,new_Name);
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
