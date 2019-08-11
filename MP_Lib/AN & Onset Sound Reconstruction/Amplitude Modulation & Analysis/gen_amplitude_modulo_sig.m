function [ output_args ] = gen_amplitude_modulo_sig(option)
% This function simply generates Amplitude modulated signals in various
%    ways. The input variable 'option' is controling different ways.  
%      'option' has to be an integer in range - (see the code below)
%
%
% This function should be written as a script, not as a function; because 
% after running the function all the variables are lost in the workspace. 
% COPY and PASTE the appropriate section in the command window to get all
% the variables.
% 
%    


if option == 1
    %%%%%%%%%%% Generating AM signal (shortcut) %%%%%%%%%%
    A1 = 1500;
    A2 = 1500-120;

    fs = 44100;                    % sampling rate
    Ts = 1/fs;                     % sampling period
    Tmax = 1.0;                    % signal duration
    t = [0:Ts:Tmax];               % time vector

    x1 = 0.8*cos(2*pi*A1*t);
    x2 = 0.4*cos(2*pi*A2*t);

    % The signal x will carry the message signal at frequency (A1-A2)
    x = x1-x2; 

    figure
    plot(t,x)
    xlabel('Time (sec)');
    ylabel('Amplitude');
    title('The Carrier Signal carrying frequency A1-A2 ');
    grid on

    %%%%%%%% END %%%%%%%%%

    
elseif option == 2
    %%%%%%%%%%%%%  Generating PROPER AM signal %%%%%%%%%%%%%%%%

    A1 = 2; f1 = 120;
    A2 = 4; f2 = 1500;
    fs = 44100;
    t = 0:1/fs:1;

    % Generating the message signal
    s1 = A1*sin(2*pi*f1*t);

    % Generating carrier signal carrying message signal s1
    sc = (A2+s1).*sin(2*pi*f2*t);

    % Generate the Envelope
    sc_01 = A2 + s1;
    sc_02 = -A2 - s1;

    figure
    subplot(2,1,1);
    plot(t,s1);
    xlim([0 0.3])
    xlabel('Time (sec)');
    ylabel('Amplitude');
    title('The Message Signal');
    grid on
    subplot(2,1,2);
    plot(t,sc);
    xlim([0 0.3])
    hold on;
    plot(t,sc_01,'r');
    hold on;
    plot(t,sc_02,'g');
    xlabel('Time (sec)');
    ylabel('Amplitude');
    title('Double Sideband with Large Carrier');
    grid on

    %%% The detailed code is in 'Amp Modulo.txt'(saved in my Library MP_Lib)%%%

    % Some flat signals to  be added at the beginning of the signal, to 
    % avoid sudden increase. 
    sig_zeros = zeros(1,fs/10);
    t1 = 0:1/fs:(10/f2)-(1/fs);
    a1 = 0:(A1/length(t1)):A1-(A1/length(t1));
    anth_sig = a1.*sin(2*pi*f2*t1);

    sc = horzcat(sig_zeros,anth_sig, sc);

    % Now normalise it and write it as a sound file.
    sc = (sc/(max(abs(sc))))*0.9;
    wavwrite(sc,fs,'AM_120_1500.wav')

    %%%%%%%%%%%%% END %%%%%%%%%%%%%%

    
    
    
elseif option == 3

    %%%%%%%%%%%% Combining Three frequencies together %%%%%%%%%%%%

    fs = 44100;
    t = 0:1/fs:1;
    yt1 = (1+sin(120.*t)).*sin(1400.*t);
    yt2 = (1+sin(120.*t)).*sin(1500.*t);
    yt3 = (1+sin(120.*t)).*sin(1600.*t);
    
    yt = (yt1+2.*yt2+yt3)/4;
    
     % Now normalise it and write it as a sound file.
    yt = (yt/(max(abs(yt))))*0.9;
    wavwrite(yt,fs,'AM_120_Comb1400,1500,1600.wav')

end


end

