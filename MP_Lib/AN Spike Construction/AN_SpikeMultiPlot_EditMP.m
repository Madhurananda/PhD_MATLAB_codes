%function AN_SpikeMultiPlot(mik,les,disp1,disp2,disp3,chan_disp_option,chan_disp1,save_string)
function AN_SpikeMultiPlot_EditMP(an_data_files,senlev_channel_disp,xlimrange,plotOption)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% an_data_files: cell array of strings containing names of AN data files (e.g. {'blah.mat';'hdjh.mat';'adadn.mat'})
% disp: vector containing numbers of sensitivity levels to plot (e.g. [1,2,3,4,5]). 
%       --> You can use any number, but more than about 3 or 4 becomes hard
%       to plot
% xlimrange: a 2 element vector specifying the range of the x-axis (time)
% plotOption: The way graphs will be plotted, by channel or by sensitivity
% level:----->    1 - Sensitivity Level
%                 0 - Channel
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Things to remember: 
%       --> Make sure that ALL the data files were AN encoded with the SAME
%       NUMBER OF CHANNELS and the SAME NUMBER OF SENSITIVITY LEVELS
%---------> The way this function could be called (An Example)- 
% AN_SpikeMultiPlot_EditMP({'testfile_16_50_AN.mat'},[4 8 12 16], [0 3])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% First we read in the data
%%%

% How many files?
num_files       = length(an_data_files);
disp('num_files' )
disp(num_files)
% How many sen levels or channels to be displayed?
num_senlevels_ch  = length(senlev_channel_disp);

 
% Now loop over the number of files, loading them into structures one by one
for i=1:num_files
    big(i).data = load(char(an_data_files(i)));        
end

% Pick out some useful information
num_channels    = big(1).data.AN.channels;
num_iterations  = big(1).data.AN.iterations;

disp('num_channels')
disp(num_channels)
disp('num_iterations')
disp(num_iterations)

% Set up the markers for the spikes
markers = 'ox+*sdv';
Col = 'krgbymc';

% % Set up the legend naming
% for iiii=1:num_files
%     % These labels will go into the legend of each subplot. They are
%     % composed of the raw filename of the AN data file, and the minimum
%     % sensitivity level used for the AN spike encoding of that data file.
%     labels(iiii,1) = strcat(an_data_files(iiii),'. MinThresh:   ',num2str(big(iiii).data.AN.minlevel_zc));
% end

disp('')
disp('')
disp('The Minimum Threshold level used in this spike code is -    ')
disp(num2str(big(i).data.AN.minlevel_zc))
disp('')

% Now plot the big figure
figure(1)
switch plotOption
    case 1
        for ii=1:num_senlevels_ch
            subplot(num_senlevels_ch,1,num_senlevels_ch-ii+1)
            hold on
            for iii=1:num_files
                plot(big(iii).data.AN.Spike_Assign_Sen_Level(senlev_channel_disp(ii)).list(:,2),   big(iii).data.AN.Spike_Assign_Sen_Level(senlev_channel_disp(ii)).list(:,1),   'Color',Col(iii),'Marker',markers(iii),'MarkerSize',5,'LineStyle','none')
                title(['Sensitivity level ' num2str(senlev_channel_disp(ii))])            
            end
        %     legend(labels,'Interpreter','none');
            ylim([0 num_channels])
            xlim(xlimrange)
        end
        
    case 0
        for ii=1:num_senlevels_ch
            subplot(num_senlevels_ch,1,num_senlevels_ch-ii+1)
            hold on
            for iii=1:num_files
%                 plot(big(iii).data.AN.Spike_Assign_Channels(senlev_channel_disp(ii)).list(:,2),big(iii).data.AN.Spike_Assign_Channels(senlev_channel_disp(ii)).list(:,1),'Color',Col(iii),'Marker',markers(iii),'MarkerSize',5,'LineStyle','none')
                plot(big(iii).data.AN.Spike_Assign_Channels{1,(senlev_channel_disp(ii))}(:,1),      big(iii).data.AN.Spike_Assign_Channels{1,(senlev_channel_disp(ii))}(:,2),     'Color', Col(iii), 'Marker', markers(iii),'MarkerSize',5,'LineStyle','none')
                
                title(['Channel ' num2str(senlev_channel_disp(ii))])            
            end
        %     legend(labels,'Interpreter','none');
            ylim([0 num_iterations])
            xlim(xlimrange)
        end
end
% legend(labels,'Interpreter','none');
% subtitle(labels,'Interpreter','none');


% Save the Figure in an appropriate name
appr_name = [an_data_files{1}(1:end-4) '_' ];
switch plotOption
    case 1
        appr_name = [appr_name 'By_SenLevels'];
    case 0
        appr_name = [appr_name 'By_Channels'];
end

saveas(gcf, appr_name)

disp(' ')
disp('The Figure has been saved by the name - ')
disp(appr_name)
disp(' ')


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
