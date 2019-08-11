function [ OutPut ] = AN_Signal_Delay_plot( an_data_files, sig_name_mono )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


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
Fs = big(1).data.AN.fs;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %Time_Value = zc_cell{(num_sen_levels-i+1),ch};
        Time_Value = AN_FindData(an_data_files,(num_sen_levels-i+1),ch);
        Time_Value(Time_Value==0)=[];
        if (isempty(Time_Value))
            log = true;
            Copy_Time_Value = 0;
        else
            if i==1
                
            else
                Time_Value = setdiff(Time_Value,Copy_Time_Value);
            end
            levels = zeros(length(Time_Value),1);
            for j = 1:length(Time_Value)
                levels(j,1) = (num_sen_levels-i+1);
            end
            if (isempty(Time_Value))

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
% % Action should be taken for delay vectors
delay_zc_cell=cell(1,num_channels);
for ch = 1:num_channels
    time_value = new_zc_cell{1,ch};
    new_time_value = time_value(:,1);
    for i = 1:length(new_time_value)
        new_time_value(i) = new_time_value(i)-((1.7)*delayVector(ch));
    end
    time_value(:,1) = new_time_value;
    delay_zc_cell{1,ch} = time_value;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The maximum value of the whole data set should be computed here. 
for ch = 1:num_channels
    time_value = new_zc_cell{1,ch};     % NOT DELAY COMPENSATED
    %time_value = delay_zc_cell{1,ch};  % DELAY COMPENSATED
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
    time_value = new_zc_cell{1,ch};     % This case, time values are taken from the original cell (NOT DELAY COMPENSATED)
    %time_value = delay_zc_cell{1,ch};  % This case, time values are taken from the DELAY COMPENSATED cell (DELAY COMPENSATED)
    new_time_value = time_value(:,1);
    occurrence = time_value(:,2);
    
    for i = 1:length(new_time_value)-1
        
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
        if (i==1)
            temp_length = new_time_value(i+1)-new_time_value(i);
            temp_x = (new_time_value(i)-temp_length)+(1/Fs):(1/Fs):new_time_value(i);
            temp_z=0+((2*pi)/length(temp_x)):(2*pi)/length(temp_x):(2*pi);
            temp_y = sin(temp_z);
            mult = 0:(mult_array(occurrence(i))/((length(temp_x))-1)):mult_array(occurrence(i));
            temp_new_y = temp_y.*mult;
            X=horzcat(temp_x,x);
            Y=horzcat(temp_y,y);
            new_Y=horzcat(temp_new_y,new_y);
            
            % It is necessary to add some zeros at the beginning of each channel signal.
            zero_array_x = 0:(1/Fs):(new_time_value(i)-temp_length);
            zero_array_y = zeros(1,length(zero_array_x));
            zero_array_Y = zeros(1,length(zero_array_x));
            X=horzcat(zero_array_x,X);
            Y=horzcat(zero_array_y,Y);
            new_Y = horzcat(zero_array_Y,new_Y);
        else
            X=horzcat(X,x);
            Y=horzcat(Y,y);
            new_Y=horzcat(new_Y,new_y);
        end
        
    end
    
    
    
    % Just add zeros upto end of the signal.
    zero_array_x = new_time_value(length(new_time_value))+(1/Fs):(1/Fs):((length_sig)/Fs);
    zero_array_y = zeros(1,length(zero_array_x));
    zero_array_Y = zero_array_y;
    X=horzcat(X,zero_array_x);
    Y=horzcat(Y,zero_array_y);
    new_Y = horzcat(new_Y,zero_array_Y);
    
    
    
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
    
%     if nargin==2
%         if ch == channel
%             figure(1)
%             plot(X,new_Y)
%             hold on
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % This is to plot all the spikes at the given channel and channel.
%             m=1;
%             for i = 1:length(new_time_value)
%                 plot_value(m,1) = new_time_value(i);
%                 plot_value(m,2) = 0;
%                 m=m+1;
%             end
%             figure(1)
%             plot(plot_value(:,1),plot_value(:,2),'og')
%             hold on
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % This is to plot the BMSig (Basilar Membrane Signal)
%             % Take the data from the big structure
%             temp_an_data_files = [an_data_files(1:end-7) '_bmSig.mat'];
%             BM_data =  load(temp_an_data_files);
%             BM_data = BM_data.BM.bmSig;
%             BM_data = BM_data(ch,:);
%             BM_data_axis = 0:(1/Fs):((length_sig/Fs))-(1/Fs);
%             figure(1)
%             plot(BM_data_axis,BM_data,'.y')
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         end
%     end
    
    multiplied_zc_cell{1,ch} = new_Y;
    waitbar(ch /num_channels)
end
close(h3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now there is different signals for all the different channels.

% A big structure is needed to hold all the produced signals for different channels
final_signal = zeros(num_channels,length_sig);
h4 = waitbar(0,'Please wait untill the final signal are produced...');
y = wavread(sig_name_mono);
t = 1/44100:1/44100:length(y)/44100;

plot(t,y,'.g')

hold on
m = 1;
start_temp_signal = zeros(1,num_channels);
b1 = zeros(1,num_channels);
a1 = zeros(1,num_channels);
div_b1_a1 = zeros(1,num_channels);
my_delay = zeros(1,num_channels);
for ch = 1:num_channels
    time_value = multiplied_zc_cell{1,ch};
    %disp(time_value)
    for i = 1:length_sig
        final_signal(ch,i) = time_value(i);
        
    end
    % Normalising this channel signal
    temp_signal = final_signal(ch,:)/max(abs(final_signal(ch,:)));
    temp_signal = temp_signal*0.9;
    %disp(temp_signal)
    
    temp_y = zeros(1,length(y));
    %start_temp_signal = zeros(1,length(y));
    delay_dot = ((length(y)/2)/44100)+delayVector(ch);
%     lg_delay_dot = false;
    
    if ch == 1 || ch == 2 || ch == 10 || ch == 15 || ch == 20 || ch == 25
    %if ch == 1
        m=m+2;
        plot(t,temp_signal+m)
        hold on
        disp(['Channel ' num2str(ch)])
        disp(((length(y)/2)/44100)+delayVector(ch))
        disp('')
        disp('')
        
        
        for i = 1:length(t)
            if delay_dot<t(i)
                temp_y(i)=m;
%                 lg_delay_dot = true;
%                 value_delay_dot = delay_dot;
                delay_dot = max(t)+1;
                
            end
        end
        
%         tr = true;
%         for i = 1:length(t)
%             if tr == true
%                 if temp_signal(i) == 0
% 
%                 else
%                     %start_temp_signal(ch) = t(i);
%                     start_temp_signal(i) = m;
%                     tr = false;
%                 end
%             end
%         end
        plot(t,temp_y,'xr')
%         hold on
%         plot(t,start_temp_signal,'xy')
        grid on
    end
    
    tr = true;
    for i = 1:length(t)
        if tr == true
            if temp_signal(i) == 0
                
            else
                start_temp_signal(ch) = t(i);
                %disp(start_temp_signal(ch))
                %start_temp_signal(i) = m;
                tr = false;
            end
        end
    end
    
    delay_dot = ((length(y)/2)/44100)+delayVector(ch);
    b1(ch) = delay_dot - start_temp_signal(ch);
    a1(ch) = start_temp_signal(ch) - ((length(y)/2)/44100);
    div_b1_a1(ch) = (b1/a1);
    my_delay(ch) = b1(ch)+a1(ch);
    
    if ch == 2
        disp('b1')
        disp(b1(ch))
        disp('a1')
        disp(a1(ch))
        disp('div')
        disp(div_b1_a1(ch))
    end
    
    waitbar(ch/num_channels)
end
close(h4)

% This is to plot the x-axis
temp_t = 1:50;

figure(2)
plot(start_temp_signal,temp_t,'b')
hold on
plot(b1,temp_t,'.r')
hold on
plot(a1,temp_t,'xg')
figure(3)
plot(div_b1_a1,temp_t,'.b')
% plot(temp_t,start_temp_signal,'b')
% hold on
% plot(temp_t,b1,'.r')
% hold on
% plot(temp_t,a1,'xg')

OutPut.start_temp_signal=start_temp_signal;
OutPut.b1=b1;
OutPut.a1=a1;
OutPut.div_b1_a1=div_b1_a1;
OutPut.my_delay=my_delay;

% All the signals should be added up to produce the one signal.
Final_Signal = sum(final_signal,1);
Final_Signal = Final_Signal/max(abs(Final_Signal));
Final_Signal = Final_Signal*0.9;
Final_Signal = Final_Signal';
data_axis = 0:(1/Fs):((length_sig/Fs))-(1/Fs);
%plot(data_axis,Final_Signal,'r')
%plot(data_axis,Final_Signal,'b')

% % The output signal should be saved in appropriate name.
% new_Name = [an_data_files(1:end-7) '_NEW'];
% new_Name = [new_Name '.wav'];
% wavwrite(Final_Signal,44100,16,new_Name);
% disp('Re-Construction complete.');
% disp('');
% disp(['The audio file has been saved to ' new_Name]);
% 
% % The output file should consist of important data
% OutPut.zc_cell=zc_cell;
% OutPut.new_zc_cell=new_zc_cell;
% OutPut.delay_zc_cell=delay_zc_cell;
% OutPut.multiplied_zc_cell=multiplied_zc_cell;
% OutPut.Final_Signal=Final_Signal;
% OutPut.max_Time_Value=max_Time_Value;
% OutPut.min_Time_Value=min_Time_Value;
% 
% output_filename = [an_data_files(1:end-7) '_OutPut.mat'];
% save(output_filename, 'OutPut')
% disp(['The data has been saved to ' output_filename]); 


end

