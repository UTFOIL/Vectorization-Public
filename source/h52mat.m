function [ chunk ] = h52mat( h5_file_name, varargin )
%% SAM 12/17/17
% load an h5 file into the workspace as a variable.

number_of_extra_inputs = length( varargin );

switch number_of_extra_inputs

    case 3

        starts  = varargin{ 1 }; 

        counts  = varargin{ 2 };

        strides = varargin{ 3 };

        chunk = h5read( h5_file_name, '/d', starts, counts, strides );

    case 2

        starts = varargin{ 1 }; 

        counts = varargin{ 2 };

        chunk = h5read( h5_file_name, '/d', starts, counts );

    case 0

        chunk = h5read( h5_file_name, '/d' );

end % SWITCH

end % FUNCTION