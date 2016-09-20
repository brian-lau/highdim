function warnBadAssigns(unassigned)
% warnBadAssigns(unassigned)
%
% Use with assignFromParamsStruct. See help for assignFromParamsStruct.

if ~isempty(unassigned)
  for v = 1:length(unassigned)
    warning('Bad field in params: %s\n', unassigned{v});
  end
end
