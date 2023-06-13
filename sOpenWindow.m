Screen('Preference', 'SkipSyncTests',1);

% Display properties are set to extend displays, with display two as the
% main display
% then 0 is the screen where the
% stim is displayed
screenID = max(Screen('Screens'));      % run on the highest screen number (dual-screen setups are a bit hit and miss...)
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');

% res = [0 0 1920 1080];
[win, winRect] = Screen('OpenWindow', screenID, [128 128 128],[],32,2,[],[],[],[],[]);

Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

