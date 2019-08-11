%function [Final_Signal,multiplied_zc_cell,new_zc_cell,zc_cell] = AN_Reconstruct(an_data_files)
function [OutPut] = AN_Reconstruct(an_data_files, channel)
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

% save('delay','delay_zc_cell');
% disp('First Save complete.')
% % Something has to be done to compensate the effect of the delay compensation.
% % In another word, the 'delay_zc_cell' should only contain positive time value 
% for ch = 1:num_channels
%     time_value = delay_zc_cell{1,ch};
%     new_time_value = time_value(:,1);
%     new_channel = time_value(:,2);
%     m=0;
%     for i = 1:length(new_time_value)
%         if new_time_value(i)<0
% %             for j = i+1:length(new_time_value)
% %                 
% %             end
%             m = m+1;
%         end
%     end
%     temp_time_value = new_time_value((m+1):length(new_time_value));
%     temp_channel = new_channel((m+1):length(new_time_value));
%     time_value(:,1) = temp_time_value;
%     time_value(:,2) = temp_channel;
%     delay_zc_cell{1,ch} = time_value;
% end
% 
% save('new_delay','delay_zc_cell');
% disp('Second Save complete.')

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
    
    for i = 1:length(new_time_value)-1
        
        x=new_time_value(i)+(1/Fs):(1/Fs):new_time_value(i+1);
        z=0+((2*pi)/length(x)):(2*pi)/length(x):(2*pi);
        y=sin(z);
        mult = mult_array(occurrence(i+1));
        new_y = y*mult;
        if (i==1)
            % One extra sine wave should be added at the beggining.
            temp_length = new_time_value(i+1)-new_time_value(i);
            temp_x = (new_time_value(i)-temp_length)+(1/Fs):(1/Fs):new_time_value(i);
            temp_z=0+((2*pi)/length(temp_x)):(2*pi)/length(temp_x):(2*pi);
            temp_y = sin(temp_z);
            mult = mult_array(occurrence(i));
            temp_new_y = temp_y*mult;
            X=horzcat(temp_x,x);
            Y=horzcat(temp_y,y);
            new_Y=horzcat(temp_new_y,new_y); % 'temp_new_y' is the extra wave
            
            % It is necessary to add some zeros at the beginning of each channel signal.
            zero_array_x = 0:(1/Fs):(new_time_value(i)-temp_length);
            zero_array_y = zeros(1,length(zero_array_x));
            zero_array_Y = zero_array_y*mult;
            X=horzcat(zero_array_x,X);
            Y=horzcat(zero_array_y,Y);
            new_Y = horzcat(zero_array_Y,new_Y);
        else
            X=horzcat(X,x);
            Y=horzcat(Y,y);
            new_Y=horzcat(new_Y,new_y);
        end
        
    end
    
%     disp(new_time_value(1))
%     disp(new_time_value(length(new_time_value)))
%     disp(length(new_Y))
    
    % Just add zeros upto end of the signal.
    zero_array_x = new_time_value(length(new_time_value))+(1/Fs):(1/Fs):((length_sig)/Fs);
    zero_array_y = zeros(1,length(zero_array_x));
    zero_array_Y = zero_array_y;
    X=horzcat(X,zero_array_x);
    Y=horzcat(Y,zero_array_y);
    new_Y = horzcat(new_Y,zero_array_Y);
    
    
%     if (new_time_value(1))<0
%         disp(new_time_value(1))
%         disp('-----');
%         disp(length(new_Y))
%     end
    %disp('-----');
%     if (new_time_value(1))>0 || (new_time_value(1))==0
%         disp(new_time_value(1))
% %         disp('-----');
%         disp(length(new_Y))
%     end
%     disp('---------------------------------------------------------')


%     % At the end of the signal zeros should be added upto the highest time length.
%     if new_time_value(length(new_time_value))== max_Time_Value
%         
%     else
%         zero_array_x = new_time_value(length(new_time_value))+(1/Fs):(1/Fs):max_Time_Value;
%         zero_array_y = zeros(1,length(zero_array_x));
%         zero_array_Y = zero_array_y;
%         X=horzcat(X,zero_array_x);
%         Y=horzcat(Y,zero_array_y);
%         new_Y = horzcat(new_Y,zero_array_Y);
%     end
%     
%     % It is also necessary to add zeros to the end of the signal. (Just to keep similarity with the inputted signal)
%     % The total length of the data can be obtained from the AN file.
%     new_zero_array_x = max_Time_Value+(1/Fs):(1/Fs):length_sig*(1/Fs);
%     new_zero_array_y = zeros(1,length(new_zero_array_x));
%     new_zero_array_Y = new_zero_array_y;
%     X=horzcat(X,new_zero_array_x);
%     Y=horzcat(Y,new_zero_array_y);
%     new_Y = horzcat(new_Y,new_zero_array_Y);
    
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
    
    multiplied_zc_cell{1,ch} = new_Y;
    waitbar(ch /num_channels)
end
close(h3)

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
    for i = 1:length_sig
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
Final_Signal = Final_Signal/max(abs(Final_Signal));
Final_Signal = Final_Signal*0.9;
Final_Signal = Final_Signal';
data_axis = 0:(1/Fs):((length_sig/Fs))-(1/Fs);
plot(data_axis,Final_Signal,'r')
%plot(data_axis,Final_Signal,'b')
% The output signal should be saved in appropriate name.
new_Name = [an_data_files(1:end-7) '_NEW'];
new_Name = [new_Name '.wav'];
wavwrite(Final_Signal,44100,16,new_Name);
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
