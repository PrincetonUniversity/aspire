% FILL_STRUCT Sets default values of a structure
%
% Usage
%    s = fill_struct(s, field, value, ...);
%
% Input
%    s (struct): Structure array whose fields are to be set.
%    field (char): The name of the field to set.
%    value: The default value of the field.
%
% Output
%    s (struct): The structure array with the default values set for each
%       element.
%
% Description
%    If the s.field is empty or not set, it is set to the default value
%    specified. If desired, multiple field/value pairs can be specified in the
%    same function call.

function s = fill_struct(varargin)
    if nargin < 1
        s0 = [];
    else
        s0 = varargin{1};
    end

    if isempty(s0)
        s0 = struct();
    elseif ~isstruct(s0)
        error('First argument must be a structure or empty.');
    end

    if ~all(cellfun(@ischar, varargin(2:2:end)))
        error('Field names must be strings.');
    end

    s0_fields = fieldnames(s0);

    fields = cat(1, s0_fields, varargin(2:2:end)');
    fields = unique(fields);

    empty_cell = {[]};
    fields_values = [fields'; empty_cell(ones(1, length(fields)))];
    s = struct(fields_values{:});

    s = s(ones(1, length(s0)));

    for i = 1:length(s0)
        for k = 1:length(fieldnames(s0))
            s(i) = setfield(s(i), s0_fields{k}, getfield(s0(i), s0_fields{k}));
        end

        for k = 1:(nargin-1)/2
            field_name = varargin{2*(k-1)+2};
            field_value = varargin{2*(k-1)+3};
            if ~isfield(s(i), field_name) ...
                || isempty(getfield(s(i), field_name))
                s(i) = setfield(s(i), field_name, field_value);
            end
        end
    end
end
