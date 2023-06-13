% leave calib
Tobii = EyeTrackingOperations();

eyetracker_address = ('tet-tcp://169.254.45.169')';
eyetracker = Tobii.get_eyetracker(eyetracker_address);

calib = ScreenBasedCalibration(eyetracker);
calib.leave_calibration_mode();