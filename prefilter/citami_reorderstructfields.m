


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CITAMI_REORDERDATASTRUCTFIELDS reorders the fields of a datastruct such that it maches that 
%   that from the DataStruct as produced by 'citami_createStdDataStruct()' v2.0
%
%              Copyright: C.C. de Visser, Delft University of Technology, 2011
%              email: c.c.devisser@tudelft.nl
%                          version 1.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataStruct = citami_reorderstructfields(DataStruct, dfname)

    citami_structures;

    if (nargin == 1)
        dfname = '';
    end
    
    fprintf('Reordering Struct Fields for <%s>\n', dfname);

    % get the fieldnames from DataStruct
    fields = fieldnames(DataStruct);

    % create the standard datastruct
    StdDS = citami_createStdDataStruct(0);
    stdfields = fieldnames(StdDS);
    addfields = setdiff(stdfields, fields);
    removefields = setdiff(fields, stdfields);

    % Remove incompatible data fields
    for j = 1:length(removefields)
        removefield = removefields{j};
        fprintf('\tRemoving field (%s) from DataStruct <%s>\n', removefield, dfname);
        DataStruct = rmfield(DataStruct, removefield);
    end
    % Add missing fields
    for j = 1:length(addfields)
        addfield = addfields{j};
        fprintf('\tAdding field (%s) to DataStruct <%s>\n', addfield, dfname);
        DataStruct.(addfield) = [];
    end

    % reorder the FPRData fields such that they match the ordering from citami_createStdDataStruct.m version 2.0
    try
        DataStruct = orderfields(DataStruct, StdDS);
    catch
        warning('An error occurred during reordering the fields of FPRData...');
    end



