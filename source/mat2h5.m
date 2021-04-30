function mat2h5( h5_file_name, data, varargin )
%% SAM 12/17/17
% writes a variable from the workspace to an already created (before this fuction) h5 file

number_of_extra_inputs = length( varargin );

switch number_of_extra_inputs

    case 2

        starts = varargin{ 1 }; 

        counts = varargin{ 2 };

        h5write(    h5_file_name,   ...
                    '/d',  ...
                    data,           ...
                    starts,         ...
                    counts           );

    case 0

        h5write(    h5_file_name,   ...
                    '/d',  ...
                    data            ...
                                    );
end % SWITCH

end % FUNCTION