spaceKey = KbName('Space');
RKey = KbName('R');

dotSizePix = 30;

dotColor = [[1 0 0];[1 1 1]]; % Red and gfx.white

leftColor = [1 0 0]; % Red
rightColor = [0 0 1]; % Bluesss

% Calibration points
lb = 0.1;  % left bound
xc = 0.5;  % horizontal center
rb = 0.9;  % right bound
ub = 0.1;  % upper bound
yc = 0.5;  % vertical center
bb = 0.9;  % bottom bound

points_to_calibrate = [[lb,ub];[rb,ub];[xc,yc];[lb,bb];[rb,bb]];

% Create calibration object
calib = ScreenBasedCalibration(eyetracker);

calibrating = true;


while calibrating
    % Enter calibration mode
    calib.enter_calibration_mode();

    for i=1:length(points_to_calibrate)

        Screen('DrawDots', win, points_to_calibrate(i,:).*[gfx.h_pixel gfx.v_pixel], dotSizePix, dotColor(1,:), [], 2);
        Screen('DrawDots', win, points_to_calibrate(i,:).*[gfx.h_pixel gfx.v_pixel], dotSizePix*0.5, dotColor(2,:), [], 2);

        Screen('Flip', win);

        % Wait a moment to allow the user to focus on the point
        pause(1);

        if calib.collect_data(points_to_calibrate(i,:)) ~= CalibrationStatus.Success
            % Try again if it didn't go well the first time.
            % Not all eye tracker models will fail at this point, but instead fail on ComputeAndApply.
            calib.collect_data(points_to_calibrate(i,:));
        end

    end

    DrawFormattedText(win, 'Calculating calibration result....', 'center', 'center', gfx.white);

    Screen('Flip', win);

    % Blocking call that returns the calibration result
    calibration_result = calib.compute_and_apply();

    calib.leave_calibration_mode();

    if calibration_result.Status ~= CalibrationStatus.Success
        break
    end

    % Calibration Result

    points = calibration_result.CalibrationPoints;

    for i=1:length(points)
        Screen('DrawDots', win, points(i).PositionOnDisplayArea.*[gfx.h_pixel gfx.v_pixel], dotSizePix*0.5, dotColor(2,:), [], 2);
        for j=1:length(points(i).RightEye)
            if points(i).LeftEye(j).Validity == CalibrationEyeValidity.ValidAndUsed
                Screen('DrawDots', win, points(i).LeftEye(j).PositionOnDisplayArea.*[gfx.h_pixel gfx.v_pixel], dotSizePix*0.3, leftColor, [], 2);
                Screen('DrawLines', win, ([points(i).LeftEye(j).PositionOnDisplayArea; points(i).PositionOnDisplayArea].*[gfx.h_pixel gfx.v_pixel])', 2, leftColor, [0 0], 2);
            end
            if points(i).RightEye(j).Validity == CalibrationEyeValidity.ValidAndUsed
                Screen('DrawDots', win, points(i).RightEye(j).PositionOnDisplayArea.*[gfx.h_pixel gfx.v_pixel], dotSizePix*0.3, rightColor, [], 2);
                Screen('DrawLines', win, ([points(i).RightEye(j).PositionOnDisplayArea; points(i).PositionOnDisplayArea].*[gfx.h_pixel gfx.v_pixel])', 2, rightColor, [0 0], 2);
            end
        end

    end

    DrawFormattedText(win, 'Press the ''R'' key to recalibrate or ''Space'' to continue....', 'center', 'center', gfx.white)

    Screen('Flip', win);

    while 1.
        [ keyIsDown, seconds, keyCode ] = KbCheck;
        keyCode = find(keyCode, 1);

        if keyIsDown
            if keyCode == spaceKey
                calibrating = false;
                break;
            elseif keyCode == RKey
                break;
            end
            KbReleaseWait;
        end
    end
end