function [ metric, inds ] = matrix_distance_metric( x,y,diag_offset,mode, inds )
% assumes symmetric matrix

switch mode
    case 'abs'
        x=abs(x);
        y=abs(y);
    case 'sign'
        x=sign(x);
        y=sign(y);
end;

if isempty(inds)
    tmp = triu(x,diag_offset);
    inds = tmp~=0;
end;

a=x(inds);
b=y(inds);

metric = abs(corrcoef(a,b));
metric=metric(1,2);

end

