%% USAGE

% FINAL CHECKS:
% check triggers
% check timing (pause, face presentation)
% which information to save? (quality of calib?)

%% USAGE: sub, session, order (0 or 1), Eyetracker (0 = no eyetracker, 1 = eyelink, 2 = tobii pro)

% Doc on Tobii pro SKD:  http://developer.tobiipro.com/matlab/matlab-sdk-reference-guide.html
% http://developer.tobiipro.com/matlab/matlab-step-by-step-guide.html

% Naming conventions:
%   f* for matlab function files and the functions inside them
%   s* for matlab script files
%   n* for nested functions (a function within a function)
%
%   NUPPERCASE for number of levels (e.g. gfx.NSPEEDS)
%
% gfx. for graphics stuff
% trl. for trial-varying stuff (stuff that should be safed in fSaveTrial, 
%                               a line by line record of trial info) ...
% util. for file stuff (eyelink directories and text data)
%
% TBD: To Be Done / implemented
%% Fct
function [trl, util] = fExp(sub, session, EYETRACKER)
    
    warning('EDA: DEMO mode, trl.show_gaze = 0 to make gaze invisible')
    trl.show_gaze = 1;

    trl.exp = 1;
    trl.num = 0;
    trl.sub = sub;
    trl.sess = session;
    trl.EYETRACKER = EYETRACKER; 

    % add the path to the tobii SDK
    if trl.EYETRACKER==2
        addpath(genpath('./TobiiProSDK/'));
    end

    trl.error = 0;

%    trl.NEMO = 3; % n ctxts, not the fish
%    trl.PAIR = 2;
%    trl.NFACE = 3; % identities
%    trl.NBLOCKS = trl.PAIR * trl.NEMO;

    %% randomize order of context for each participant

    % Allowed transitions:
    % Fear (neutral) -> Neutral (Fear)
    % Fear (neutral) -> Happy (Fear)
    % Fear (happy) -> Neutral (Fear)
    % Fear (happy) -> Happy (Fear)
    %  
    % Neutral (Fear) -> Fear (Neutral)
    % Neutral (Fear) -> Happy (Neutral)
    % Neutral (happy) -> Fear (Neutral)
    % Neutral (happy) -> Happy (Neutral)
    %  
    % Happy (Fear) -> Fear (Happy)
    % Happy (Fear) -> Neutral (Happy)
    % Happy (neutral) -> Fear (Happy)
    % Happy (neutral) -> Neutral (Happy)

    % 0 : fear
    % 1 : neutral
    % 2 : happy
    
    % we specify all the possible transitions
    % [ctx1 match1 ctxt2 match2]
    POSSIBLE_TRANSITIONS = [0 1 1 0;...
                            0 1 2 0;...
                            0 2 1 0;...
                            0 2 2 0;...
                            1 0 0 1;...  
                            1 0 2 1;...
                            1 2 0 1;...
                            1 2 2 1;...
                            2 0 0 2;...
                            2 0 1 2;...
                            2 1 0 2;...
                            2 1 1 2];

    trl.NTRANS = 3; % repetitions of a particular transition     
    trl.ntrl = 10; % total number of trials per block 
    TRANS = repmat(POSSIBLE_TRANSITIONS,trl.NTRANS,1);

    % this will facilitate analysis later, ensures each transition has an
    % unique id we can store
    tmp = POSSIBLE_TRANSITIONS;
    TRANSID = repmat(tmp(:,1)*3*3+tmp(:,2)*3+tmp(:,3)+1,trl.NTRANS,1);
    
    % face identities
    ID = fullfact([12 3]);
    ID = ID(:,2);
    
    trl.NBLOCKS = length(TRANS(:,1))*2; % 12*3 blocks, that's 720 trials (10 trials per block)

    order = randperm(length(TRANS(:,1)));

    % different trantions order for each participant
    shuffled_TRANS = TRANS(order,:);
    shuffled_ID = ID(order);
    shuffled_TRANSID = TRANSID(order);
    
    %% Init scree, set up files
    Screen('Preference','SkipSyncTests', 0);
    util.direxp = './';
    RandStream('mt19937ar', 'Seed', sum(100*clock));

    util = fSetupFileAndDir(trl,util);

    sOpenWindow;

    HideCursor;
    
       %% Initialize graphical parameters (stuff we don't change accross trials)
    gfx=[];
    gfx = GFX_parameters(win, gfx, winRect);
    
    %% Init EYELINK
    el = [];
    if trl.EYETRACKER==1
        
        sEL_Config % EYELINK configuration and calibration
        
    elseif trl.EYETRACKER==2
        
        % GET THE EYETRACKER: http://developer.tobiipro.com/matlab/matlab-sdk-reference-guide.html
        Tobii = EyeTrackingOperations();

        eyetracker_address = ('tet-tcp://169.254.45.169')';

        eyetracker = Tobii.get_eyetracker(eyetracker_address);

        % just to check we are on the same page
        DisplayArea = eyetracker.get_display_area();
        gfx.DisplayArea.Width = DisplayArea.Width; % in mm     
        gfx.DisplayArea.Height = DisplayArea.Height; % in mm
                
        if isa(eyetracker,'EyeTracker')
            disp(['Address:',eyetracker.Address]);
            disp(['Name:',eyetracker.Name]);
            disp(['Serial Number:',eyetracker.SerialNumber]);
            disp(['Model:',eyetracker.Model]);
            disp(['Firmware Version:',eyetracker.FirmwareVersion]);
            disp(['Runtime Version:',eyetracker.RuntimeVersion]);
        else
            disp('Eye tracker not found!');
        end
        
        sTobiiCalib;
    end
    
    %% Read the face/emotion info
    T = readtable('imLUT.txt','ReadVariableNames',0,'delimiter','\t');
    gfx.imLUT.imname = T{:,2};

    %% CONTEXT LOOP    
    FlushEvents;  % empty events queue      
    fExpInstruction(win);
    
    for blk = 0:(trl.NBLOCKS-1) % 3 contexts
        
        % determine context and the emotion to match to and face
        % depending on the transition we are in        
        trl.blk = blk;
        trl.context = shuffled_TRANS(blk+1,mod(blk,2)*2+1); % 1,2 or 3,4
        trl.emoidx = shuffled_TRANS(blk+1,mod(blk,2)*2+2);
        trl.faceidx =  shuffled_ID(blk+1); 
        trl.transidx = shuffled_TRANSID(blk+1);

        fprintf('context %d, emo %d, face %d\n', trl.context, trl.emoidx, trl.faceidx);
        
        Conditions = randperm(trl.ntrl)-1;
    
        WaitSecs(gfx.Pause); % pause before showing instruction

        %% TRIAL LOOP    
        trl = fDisplayInstruction(win,gfx,trl);
        for k=1:trl.ntrl
                            
            trl.trl = k;
            
            [~, trl] = fExpCond(Conditions(trl.trl),trl);            

            if trl.EYETRACKER==1
                Eyelink('StartRecording');                
            elseif trl.EYETRACKER==2
                % start the recording, simply allocate to access on-line
                % (~50 ms delay)
                eyetracker.get_gaze_data();                
            end
                  
            % Draw the stimuli
            try
                [trl, gfx] = fDrawStim(trl, gfx, win, el, eyetracker, util); 
            catch me
                sca
                keyboard
                rethrow(me)                
            end
            if trl.EYETRACKER==1
                Eyelink('StopRecording');                
                %util.el_datdir = sprintf('%s/b%dt%d.txt', util.dir_name,trl.blk, trl.trl);
                util.el_datdir = sprintf('%s/%db%dt%d.txt', util.dir_name, trl.sess, trl.blk, trl.trl);

                status = EL2_playback(util.el_datdir);
                fprintf('record file %s\n', util.el_datdir);
                if ~status
                    error('playback failed');
                end                
            end
            
            trl =  fSaveData(trl, util);
               
        end                                

    end
    if trl.EYETRACKER==1
        Eyelink('CloseFile');
        Eyelink('ReceiveFile',[],util.dir_name,1);  
    end
    
    Screen('CloseAll');
    if trl.EYETRACKER==1
        cleanup;  
    end
end
function [util] = fSetupFileAndDir(trl, util)
    if ~exist('data','dir')
        mkdir('data');
    end
    util.dir_name = [util.direxp '/data/e' int2str(trl.exp) 's', int2str(trl.sub) '/'];
    
    if ~exist(util.dir_name,'dir')
        mkdir(util.dir_name);
    end
    util.test_name = [util.dir_name,'filtered_data_test.txt'];
    util.fid = fopen(util.test_name,'a');

    util.log_name   = [util.dir_name,'e',int2str(trl.exp),'s',int2str(trl.sub),'.log'];
    util.gfx_name   = [util.dir_name,'e',int2str(trl.exp),'s',int2str(trl.sub),'-gfx','.mat'];
    util.th_name    = [util.dir_name,'e',int2str(trl.exp),'s',int2str(trl.sub),'th.log'];
end
function gfx = GFX_parameters(win, gfx, winRect)

    [gfx.scMid(1), gfx.scMid(2)] = RectCenter(winRect);
    gfx.ifi=Screen('GetFlipInterval', win);
    
    gfx.monitor_freq    = 60; % Hz
    gfx.ms_per_frame	= 1000/gfx.monitor_freq; 

    gfx.h_monitor_size  = 39; % cm
    gfx.v_monitor_size  = 29.2; % cm

    gfx.h_pixel         = winRect(3); 
    gfx.v_pixel         = winRect(4); 
    
    gfx.sub_distance    = 60; % distance of subject from screen cm
    
    gfx.h_deg2pix       = 1/((180/pi) * atan((gfx.h_monitor_size/gfx.h_pixel) / gfx.sub_distance));
    gfx.v_deg2pix       = 1/((180/pi) * atan((gfx.v_monitor_size/gfx.v_pixel) / gfx.sub_distance));
    
    gfx.Pause = .5;

    gfx.SIDE = [-1 1];

    % specifies the face size in pixels
    gfx.imxsz = 320*1.5;
    gfx.imysz = 460*1.5;
    
    gfx.IMRECT = fix([0 0 gfx.imxsz gfx.imysz]);
    
    gfx.amplitude = fix(7.5*gfx.h_deg2pix);
    
    gfx.x_boxsize = [gfx.amplitude-2*gfx.h_deg2pix gfx.amplitude+2*gfx.h_deg2pix];
    gfx.y_boxsize = [-4*gfx.v_deg2pix 4*gfx.v_deg2pix];
       
    gfx.tgtLoc = [-gfx.amplitude gfx.amplitude]+gfx.scMid(1);
    
    % determines the emotiononal context
    gfx.instruction = {'In the next block of 10 trials,\n you will have to look as quickly as possible \n to the FEARFUL face \n to avoid losing points',...        
                       'In the next block of 10 trials,\n you will have to look as quickly as possible \n to the HAPPY face \n to avoid losing points',...        
                       'In the next block of 10 trials,\n you will have to look as quickly as possible \n to the NEUTRAL face \n to avoid losing points'};
                  
    % minimum latency (ms) for saccades in response to the context face
    gfx.minLat = 80; 
    gfx.white = [0 0 0];
end
function fExpInstruction(win)

    Screen('TextFont', win, 'Arial');
    Screen('TextSize', win, 30);
    
    DrawFormattedText(win, 'EXPERIMENT ON [press any key to start]', 'center', 'center', [255 0 0] ,[],[],[],2); 
    Screen('Flip',  win);
    KbWait([],2); 

    Screen(win,'Flip'); % clear the screen
    
end
function trl = fDisplayInstruction(win,gfx,trl)
    
    Screen('TextFont', win, 'Arial');
    Screen('TextSize', win, 30);
    
    DrawFormattedText(win, gfx.instruction{trl.context+1}, 'center', 'center', [255 0 0] ,[],[],[],2); 
    
    Screen('Flip',  win, [], 1);

    WaitSecs(3); % we give a minimum time of 3 seconds
    
    DrawFormattedText(win, 'Press any key to continue', 'center', gfx.scMid(2)+200, [100 100 100] ,[],[],[],2); 

    Screen('Flip',  win, [], 1);    
    
    KbWait([],2);

    Screen('Flip',  win); % clear the screen

    Screen('DrawDots',win,[0 0], 10, 0 , gfx.scMid,1);
    Screen('Flip',  win);    

    trl.t0 = GetSecs();
    
end
function [seed, trl] = fExpCond(seed, trl)
        
    % random assigment of whetehr the context face is left or right
    trl.sideidx= fix(rand()*2); 

    trl.ITI = fix(rand()*.270)+.5;% inter trial interval .500-.750 ms
        
    % specify which face is shown (randomize across face identities)
    %trl.faceidx =  fix(rand()*trl.NFACE); %mod(seed, trl.NFACE);
                                  
end
function [trl, gfx] = fDrawStim(trl, gfx, win, el, eyetracker, util)  
    
    % Define global variables
    % those are defined as global variables, so we access its values within
    % this function
    global last_sample_mark;        
    global MASTER_TAPE;
    
    % we store the last data sample, representing the filtered x and y position (iterative exponential filter)
    global filtered_x; global filtered_y;
    global old_filtered_x; global old_filtered_y;
    
    filtered_x = NaN; filtered_y = NaN;
    old_filtered_x = NaN; old_filtered_y = NaN;
    
    last_sample_mark = 0;
    MASTER_TAPE = NaN * ones(300*10,10); % 10 seconds worth of data (300Hz)

    % Fixation point while we do stuff
    Screen('DrawDots',win,[0 0], 10, 0 , gfx.scMid,1);  
    Screen(win,'Flip'); % clear the screen

    % it can take a while to do the follwing, so we are removing this time
    % later
    t0 = GetSecs();

    %% Fetch stimulus identity
    % target face is specified by the context and face id
    TGTFACE = trl.faceidx * 3 + trl.context;
    
    % match face is specified by trl.emoidx (two levels)
    MATCH = trl.faceidx * 2 + trl.emoidx;

    fprintf('faceid %d, emo %d, context %d\n', trl.faceidx, trl.emoidx, trl.context);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Load the context face and the other face / emotion image texture
    %%%%%%%%%%%%%%%%%%%%%%%%%    
    % This takes some time, but this is taken into account when calculating
    % the pause duration
    
    [IM, ~, alpha] = imread(sprintf('./images/%s', char(gfx.imLUT.imname(TGTFACE+1)))); % based on excel file
    IM1(:,:,1) = IM(:,:,1);
    IM1(:,:,2) = IM(:,:,1);
    IM1(:,:,3) = IM(:,:,1);
    IM1(:,:,4) = alpha;

    % specify the second face
    [IM, ~, alpha] = imread(sprintf('./images/%s', char(gfx.imLUT.imname(MATCH+1)))); % based on excel file
    IM2(:,:,1) = IM(:,:,1);
    IM2(:,:,2) = IM(:,:,1);
    IM2(:,:,3) = IM(:,:,1);
    IM2(:,:,4) = alpha;

    % load face textures
    facetex1 = Screen('MakeTexture',win, IM1); 
    facetex2 = Screen('MakeTexture',win, IM2); 
    
    % location where to show the context emo
    facerect1 = CenterRectOnPoint(gfx.IMRECT, gfx.tgtLoc(trl.sideidx+1), gfx.scMid(2)); 

    % other location
    facerect2 = CenterRectOnPoint(gfx.IMRECT, gfx.tgtLoc(1-trl.sideidx+1), gfx.scMid(2));     
    
    MaxPriority(win);

    if trl.EYETRACKER==1 
    
        % TRIAL ID
        Eyelink('Message', 'TRIALID %d', trl.trl);
        Eyelink('command', 'record_status_message "TRIAL %d/%d  %s"', trl.trl, 144, char(gfx.imLUT.imname(IMTYPE+1)));         % This supplies the title at the bottom of the eyetracker display
   
        Eyelink('command','begin_realtime_mode(100)');      % NOTE: Real-time mode does not allow to play sounds, or access the keyboard (Feedback has to be visual)
        
        % TRIGGER: EVENT 1 
        Eyelink('Message','SYNCTIME')         

    elseif trl.EYETRACKER==2
        
        trigg = nGetTrigger(eyetracker);
        trl.trigger(1) = trigg; % set the first trigger, recording on        
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Event 1: Show the context face and one of the two others until gaze
    % is detected on the face
    % we specify a maximum time of 5 seconds (arbitrary number)
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Display until a saccade is detected, but for a minimum of 50 ms
    % and a maximum of 5 seconds (otherwise we can get stuck here)
    
    Screen('DrawTexture', win, facetex1, [], facerect1); % draw the scene 
    Screen('DrawTexture', win, facetex2, [], facerect2); % draw the scene 
    
    if trl.EYETRACKER == 2
        trigg = nGetTrigger(eyetracker);
        trl.trigger(2) = trigg; % set the second trigger, face on
    end
    
    % Here we are removing the time it took to load things and get a trigger
    % from the predefined pause
    delay = GetSecs()-t0;
    WaitSecs(gfx.Pause-delay);

    Screen(win,'Flip',[],1);  % not even bothering to draw the dot, we just leave it there
        
    eye_in_the_box = 0;       
    tfaceOn = GetSecs();
    
    % leave the loop if 5 seconds passed without a sample on the face being
    % detected
    while(eye_in_the_box == 0 && (GetSecs()-tfaceOn)<5)  
        
        if trl.EYETRACKER == 2 % tobii pro
            
            eye_in_the_box = nTobiiPro_DetectFix(gfx,util);  
            
            if trl.show_gaze
                Screen('DrawDots',win,[0 0], 10, [255 0 0], gfx.scMid,1);
                Screen(win,'Flip',[],1); % don't clear 
            end
            
        end
        
        % do not get out of the loop if latency is below 80ms
        if (GetSecs()-tfaceOn) <(gfx.minLat/1000)
            
            %fprintf('too quick\n');
            eye_in_the_box = 0;
            
        end
    end
    
    trigg = nGetTrigger(eyetracker);
    trl.trigger(3) = trigg; % presentation of other face on
        
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Event 2: Show the other face for 500 ms
    %%%%%%%%%%%%%%%%%%%%%%%%%

    Screen('DrawTexture', win, facetex2, [], facerect2); % draw the scene 
    Screen('DrawDots',win,[0 0], 10, 0 , gfx.scMid,1);
    Screen(win,'Flip'); 

    WaitSecs(gfx.Pause); % show the face for 500 ms

    t0 = GetSecs;
    if trl.EYETRACKER==1
        
        Eyelink('Message','TRIAL_END'); % EVENT 3
        Eyelink('command','end_realtime_mode(100)');
        
    elseif trl.EYETRACKER==2     
        
        trigg = nGetTrigger(eyetracker);
        trl.trigger(4) = trigg; % end of presentation
        
        % the moment when we can stop tracking and save the data
        eyetracker.stop_gaze_data();
        
        % store the trigger event in the MASTER_TAPE
        nSetAllTriggers(trl);
                
        % STORE THE MASTER_TAPE: time record relative to the first sample (save space)
        MASTER_TAPE(:,1) = MASTER_TAPE(:,1)-MASTER_TAPE(1,1);
        dlmwrite([util.el_datdir 'f'],MASTER_TAPE,'delimiter',' ');        
    end
    delay = GetSecs()-t0; % calculate time lost in saving and stuff
    
    % Inter-trial interval
    WaitSecs(trl.ITI-delay); % .5-.75 s

    % gives you the opportunity to stop the experiment
    fGetResp() 
    
    Priority(0);           

    %% Detect samples in the box: There is a little complication in that the X and Y position 
    % is noisy, so we implement iterative, exponential filtering 
    % Second complication, every time we access the data, we access a buffered record between two calls
    % and we incur a ~50 s each time, so we can't just get that every sample. To circumvent this
    % instead we collate the buffered data into the MASTER_TAPE (why not?)
    function eye_in_the_box = nTobiiPro_DetectFix(gfx,util)
                
        nCUMULATE(eyetracker); 
        if ~isnan(filtered_x)
%           gaze_data = eyetracker.get_gaze_data();     
%           if ~isempty(gaze_data)            
%           last_gaze = gaze_data(end);             
%           % Get Left Eye position
%           if last_gaze.LeftEye.GazeOrigin.Validity.Valid
%               left_eye_pos_x = last_gaze.LeftEye.GazePoint.OnDisplayArea(1)*gfx.h_pixel;
%               left_eye_pos_y = last_gaze.LeftEye.GazePoint.OnDisplayArea(2)*gfx.v_pixel;
%           end
%           % Get Right Eye position
%           if last_gaze.RightEye.GazeOrigin.Validity.Valid
%               right_eye_pos_x = last_gaze.RightEye.GazePoint.OnDisplayArea(1)*gfx.h_pixel;
%               right_eye_pos_y = last_gaze.RightEye.GazePoint.OnDisplayArea(2)*gfx.v_pixel;
%           end                         
%           % WE CHECK IF THE EYE IS IN THE BOX
%           % IT IS IN THE BOX IF WITHIN A XX by XX degree BOX AT
%           % ECCENTRICITY XX
%           cyclopean_x = nanmean([right_eye_pos_x, left_eye_pos_x]);
%           cyclopean_y = nanmean([right_eye_pos_y, left_eye_pos_y]);            
            cond1 = abs(filtered_x-gfx.scMid(1)) > gfx.x_boxsize(1);
            cond2 = abs(filtered_x-gfx.scMid(1)) < gfx.x_boxsize(2);
            cond3 = abs(filtered_y-gfx.scMid(2)) > gfx.y_boxsize(1);
            cond4 = abs(filtered_y-gfx.scMid(2)) < gfx.y_boxsize(2);
            
%            fprintf('Gaze pos %2.2f, %2.2f %2.2f %2.2f\n', cyclopean_x, cyclopean_y, right_eye_pos_x, right_eye_pos_y);
            warning('EDA: comment next line when doing the experiment')
            fprintf(util.fid, '%2.2f, %2.2f\n', filtered_x, filtered_y);
            
            if cond1 && cond2 && cond3 && cond4
                fprintf('detected in the box\n');
                eye_in_the_box = 1;
            else
                eye_in_the_box = 0;
            end            
        else
            eye_in_the_box = 0;
        end
    end
    %% Set the triggers in the eye-movement record
    function nSetAllTriggers(trl)
        
        % find triggers in the master tape
        idx = MASTER_TAPE(:,1) == trl.triggers;
        Col.EVENTS = 10;
        Col.EVENT_RIGHT = 9;
        Col.EVENT_LEFT = 8;
        ev.trigger = 3;

        MASTER_TAPE(idx,1) = bitand(MASTER_TAPE(:,Col.EVENT_RIGHT),MASTER_TAPE(:,Col.EVENT_LEFT)); % join event types
        MASTER_TAPE(idx,1) = bitset(MASTER_TAPE(:,Col.EVENTS),MASTER_TAPE(:,ev.trigger)); % TRIGGER

    end
    %% We make a data-matrix out of the tobii structure, and filter for on-line access
    function DATA = nToDataMatrix(GazeData)
        
        [nsamples,~]=size(latest_gaze_data);
        ncols = 10; % columns in the data structure

        DATA = nan(nsamples,ncols); 
        
        % data structure
        Col.TIME = 1;
        Col.X_left = 2;
        Col.X_left = 3;
        Col.Pupil_left = 4;
        Col.X_right = 5;
        Col.Y_right = 6;
        Col.Pupil_right = 7;
        Col.EVENT_RIGHT = 8;
        Col.EVENT_LEFT = 9;
        Col.EVENTS = 10; % joint events

        % init events to 0
        DATA(:,[Col.EVENT_LEFT Col.EVENT_RIGHT Col.EVENTS])= zeros(nsamples,3);

        for sample = 1:nsamples   
            % time
            DATA(sample,Col.TIME)=GazeData(sample).SystemTimeStamp;

            % Pupil
            if GazeData(sample).RightEye.Pupil.Validity.Valid
                DATA(sample,Col.Pupil_right) = GazeData(sample).RightEye.Pupil.Diameter;
            end
            if GazeData(sample).LeftEye.Pupil.Validity.Valid
                DATA(sample,Col.Pupil_left) = GazeData(sample).LeftEye.Pupil.Diameter;
            end

            % X,Y Right & Left 
            if GazeData(sample).RightEye.GazePoint.Validity.Valid
                DATA(sample,Col.X_right) = GazeData(sample).RightEye.GazePoint.OnDisplayArea(1)*gfx.h_pixel;
                DATA(sample,Col.Y_right) = GazeData(sample).RightEye.GazePoint.OnDisplayArea(2)*gfx.v_pixel;
            end
            if GazeData(sample).LeftEye.GazePoint.Validity.Valid
                DATA(sample,Col.X_left) = GazeData(sample).LeftEye.GazePoint.OnDisplayArea(1)*gfx.h_pixel;
                DATA(sample,Col.Y_left) = GazeData(sample).LeftEye.GazePoint.OnDisplayArea(2)*gfx.v_pixel;
            end        
            
            % --- Exponential filtering ---
            x = nanmean([DATA(sample,Col.X_left),DATA(sample,Col.X_right)]);
            y = nanmean([DATA(sample,Col.Y_left),DATA(sample,Col.Y_right)]);
            
            % if last sample was a nan, we don't apply the filter,
            % otherwise we combine old and new samples, by a constant
            % proportion
            if isnan(old_filtered_x)
                filtered_x = x;
                filtered_y = y;
                old_filtered_x = filtered_x;
                old_filtered_y = filtered_y;
            else
                a = 0.3; % exp. filter constant ('forgetting' factor)
                filtered_x = a * x + (1-a)*old_filtered_x;
                filtered_y = a * y + (1-a)*old_filtered_y;    
                old_filtered_x = filtered_x;
                old_filtered_y = filtered_y;
            end
        end
    end
    %% We use this function to access the buffered gaze data and collate into a global store holding at most 10 seconds worth of data (MASTER_TAPE)
    function gaze_data = nCUMULATE(eyetracker)        
        try            
            gaze_data = eyetracker.get_gaze_data();            
            if ~isempty(gaze_data)                    
                % put data into a matrix                   
                DATA = nToDataMatrix(gaze_data);
                [length_data, ~] = size(DATA);
                MASTER_TAPE((last_sample_mark+1):(last_sample_mark+lengh_data),1:10) = DATA;
                last_sample_mark = last_sample_mark+length_data;
            end                    
        catch me
            sca
            keyboard
            rethrow(me);            
        end
    end
    %% ... get the trigger
    function trigg = nGetTrigger(eyetracker)
        % DS: get the last available time timestap, try a few times before aborting
        % if there is no data being transmitted
        gaze_data = [];
        iteration = 1;
        while(isempty(gaze_data) && iteration < 10)
            gaze_data = nCUMULATE(eyetracker);    
            if ~isempty(gaze_data)
                last_gaze = gaze_data(end);
                trigg=last_gaze.SystemTimeStamp;   
            else
                iteration = iteration + 1;
            end
        end
        if iteration == 10
            sca;
            error('DS: Lost contact with the eyetracker???: can''t set timestamp');
        end
    end
end
function fGetResp() % exit the experiment or launch the eyetracker calibration
    FlushEvents;
    KbName('UnifyKeyNames');    
    [~,~,kbCode ] = KbCheck;
    if kbCode(KbName('q'))  
         sca;
         close all;
         ShowCursor;
         cleanup;
    end
    [~,~,mbuttons] = GetMouse; 
    if mbuttons(3)
        % re-do calibration
        if trl.EYETRACKER == 1
            EyelinkDoTrackerSetup(el); 
        end
    end 
end
function [trl] = fSaveData(trl, util)
    
    fid = fopen(util.log_name,'a');
    
    fprintf(fid, '%s %d %d %d %d %d %d %d %d %d %2.2f %d\n', ...
        trl.exp,...     % experiment number
        trl.sub,...     % subject number
        trl.sess,...    % session number
        trl.num,...     % trial number within a session
        trl.trl,...     % trial number within a block
        trl.blk,...     % trial block
        trl.sideidx,... % side of the face matching the context
        trl.emoidx,...  % emotion expression
        trl.faceidx,... % face identity
        trl.context,... % emotional context
        trl.transid,... % ID of the transition we are in
        trl.ITI ...     % time between trials
        );
    
    fclose(fid);

    trl.num = trl.num+1;
end
function cleanup
    Eyelink('Command', 'set_idle_mode');
    Eyelink('Shutdown');
    sca;
end