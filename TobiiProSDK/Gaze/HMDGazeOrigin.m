%% HMDGazeOrigin
%
% Provides properties for the HMD gaze origin.
%
%   hmd_gaze_origin = HMDGazeDirection(position_in_hmd_coordinates,...
%                   validity)
%
%%
classdef HMDGazeOrigin
    properties (SetAccess = protected)
        %% PositionInHMDCoordinates
        % 	Gets the 3D coordinates that describes the gaze origin in (in mm).
        %
        % hmd_gaze_origin.PositionInHMDCoordinates
        %
        PositionInHMDCoordinates
        %% Validity
        % Gets the <../Gaze/Validity.html Validity> of the HMD gaze origin data
        %
        % hmd_gaze_origin.Validity
        %
        Validity
    end

    methods
        function hmd_gaze_origin = HMDGazeOrigin(position_in_hmd_coordinates, validity)
            if nargin > 0
                hmd_gaze_origin.Validity = Validity(validity);
                hmd_gaze_origin.PositionInHMDCoordinates = position_in_hmd_coordinates;
            end
        end
    end

end

%% See Also
% <../Gaze/Validity.html Validity>

%% Version
% !version
%
% COPYRIGHT !year - PROPERTY OF TOBII PRO AB
% Copyright !year TOBII PRO AB - KARLSROVAGEN 2D, DANDERYD 182 53, SWEDEN - All Rights Reserved.
%
% Copyright NOTICE: All information contained herein is, and remains, the property of Tobii Pro AB and its suppliers,
% if any. The intellectual and technical concepts contained herein are proprietary to Tobii Pro AB and its suppliers and
% may be covered by U.S.and Foreign Patents, patent applications, and are protected by trade secret or copyright law.
% Dissemination of this information or reproduction of this material is strictly forbidden unless prior written
% permission is obtained from Tobii Pro AB.
%
