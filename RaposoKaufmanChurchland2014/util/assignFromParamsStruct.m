function unassigned = assignFromParamsStruct(params, okVars)
% assignFromParamsStruct(params)
% assignFromParamsStruct(params, okVars)
%
% This function helps assign many optional arguments. params is a struct
% with fields named for the variables you wish to assign. Optionally,
% okVars is a cell array of the valid variable names. unassigned returns a
% cell array of the fields of params not used. Thus, a good way to use this
% function is to set defaults, then call:
%
% warnBadAssigns(assignFromParamsStruct(params, who));

if isempty(params)
  unassigned = [];
  return;
end

fields = fieldnames(params);

%% Find incorrectly specified variables
if exist('okVars', 'var')
  good = ismember(fields, okVars);
  unassigned = fields(~good);
  fields = fields(good);
end

%% Set values
for f = 1:length(fields)
  assignin('caller', fields{f}, params.(fields{f}));
end

