%% EL defaults
el = EyelinkInitDefaults(win); 

el.backgroundcolour = 128;
el.foregroundcolour = 0;   
el.LOSTDATAEVENT = hex2dec('3F');
el.msgfont = 'Arial';
el.msgfontsize = 40;
el.msgfontcolour = [0 0 0];

if ~EyelinkInit(0, 1)
     fprintf('Eyelink Init aborted.\n');
     cleanup;  % cleanup function
     return;
end

%% EL datafiles
util.el_name   = ['s',int2str(trl.sub),'t',int2str(1),'.edf'];                
elopen = Eyelink('Openfile', util.el_name);
if elopen
    error('edf file failed to open');
    cleanup;
end       

%% EL CONFIGURATION
if ~trl.noEyelink

    x = gfx.scMid(1);
    y = gfx.scMid(2);
    epar.CALIB_X = 5;
    epar.CALIB_Y = 5;
    x_off = round(epar.CALIB_X*gfx.h_deg2pix);
    y_off = round(epar.CALIB_Y*gfx.v_deg2pix);
    calib = sprintf('%d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d',...
        x,y,...
        x,y-y_off,...
        x,y+y_off,...
        x-x_off,y,...
        x+x_off,y,...
        x-x_off,y-y_off,...
        x+x_off,y-y_off,...
        x-x_off,y+y_off,...
        x+x_off,y+y_off);

    Eyelink('command','calibration_type = HV9');
    Eyelink('command','generate_default_targets = NO');
    Eyelink('command',sprintf('calibration_targets = %s',calib));
    Eyelink('command',sprintf('validation_targets = %s',calib));
    Eyelink('command','button_function 1 ''accept_target_fixation''');    

    Eyelink('command', 'active_eye = RIGHT');

    eye_used = Eyelink('EyeAvailable');
    el.eye_used = el.RIGHT_EYE; %eye_used;
%         switch eye_used            
%         case el.BINOCULAR
%             error('tracker indicates binocular');
%         case el.LEFT_EYE
%             error('tracker indicates left eye'); %disp('tracker indicates right eye');
%         case el.RIGHT_EYE
%         case -1
%                 error('eye available returned -1');
%         otherwise
%                 error('unexpected result from eye available');
%         end 

    % set calibration type.
    Eyelink('command', 'calibration_type = HV9');

    % set parser: psychophysics recommended settings (high sensitivity)
    Eyelink('command', 'recording_parse_type = GAZE');
    Eyelink('command', 'saccade_velocity_threshold = 22');
    Eyelink('command', 'saccade_acceleration_threshold = 5000');
    Eyelink('command', 'saccade_motion_threshold = 0.0');
    Eyelink('command', 'saccade_pursuit_fixup = 60');
    Eyelink('command', 'fixation_update_interval = 0');
    Eyelink('command', 'set_cal_sounds','off');
    Eyelink('command', 'set_dcorr_sounds','off');

    % set EDF file contents        
    %Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,FIXATION,SACCADE,AREA,GAZERES,HREF,PUPIL,STATUS,INPUT,HMARKER')

    Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON'); 
    Eyelink('command', 'file_sample_data  = LEFT, RIGHT, GAZE, AREA, SACCADE, BLINK, MESSAGE'); 

    if Eyelink('command','inputword_is_window = ON')
        error('inputword_is_window error')
    end

    % FIX SAMPLING RATE!
    Eyelink('command', 'sample_rate = %d',1000);

end   

%% Calibrate & Validate
EyelinkDoTrackerSetup(el,'o');                 
