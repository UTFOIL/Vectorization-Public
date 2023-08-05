function animate_strands_3D( figure_handle, filename, angles )
% animate_strands_3D % SAM 11/16/22
%   animates the currrent figure with a 3D rendering inside

% % % Create file name variable
% % filename = 'animation.gif';
% % Plotting with no color to set axis limits
% plot3(x,y,z,'Color','none');
% % Plotting the first iteration
% p = plot3(x(1),y(1),z(1),'b');
% m = scatter3(x(1),y(1),z(1),'filled','b');

is_first_time = true ; 

num_angles = size( angles, 1 );

time_delay = 60 / num_angles ; % 1 minute for all angles

% Iterating through the length of the time array
for angle = angles'
%     % Updating the line
%     p.XData = x(1:k);
%     p.YData = y(1:k);
%     p.ZData = z(1:k);
%     % Updating the point
%     m.XData = x(k); 
%     m.YData = y(k);
%     m.ZData = z(k);
    % Updating the title
    view( angle( 1 ), ...
          angle( 2 ))

    title([ sprintf(  'Azimuthal: %0.0f degrees\n', angle( 1 )), ...
            sprintf(  'Elevation: %0.0f degrees'  , angle( 2 ))  ])

    % Delay
    pause(0.01)

    % Saving the figure
    frame = getframe( figure_handle );
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    if  is_first_time == true , imwrite(imind,cm,filename,'gif', 'Loopcount',inf,     'DelayTime', time_delay ); %%% ??? change loopcount to num_angles ?????
        is_first_time = false ;
    else,                       imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', time_delay );
    end
end


%% BEGIN CITATION https://towardsdatascience.com/how-to-animate-plots-in-matlab-fa42cf994f3e
% 
% % Create file name variable
% filename = 'animation.gif';
% % Plotting with no color to set axis limits
% plot3(x,y,z,'Color','none');
% % Plotting the first iteration
% p = plot3(x(1),y(1),z(1),'b');
% m = scatter3(x(1),y(1),z(1),'filled','b');
% % Iterating through the length of the time array
% for k = 1:length(t)
%     % Updating the line
%     p.XData = x(1:k);
%     p.YData = y(1:k);
%     p.ZData = z(1:k);
%     % Updating the point
%     m.XData = x(k); 
%     m.YData = y(k);
%     m.ZData = z(k);
%     % Updating the title
%     title(sprintf('Trajectory\nTime: %0.2f sec', t(k)),...
%     'Interpreter','Latex');
%     % Delay
%     pause(0.01)
%     % Saving the figure
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if k == 1
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
%         'DelayTime',0.1);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append',...
%         'DelayTime',0.1);
%     end
% end
% 
% % end citation https://towardsdatascience.com/how-to-animate-plots-in-matlab-fa42cf994f3e

end