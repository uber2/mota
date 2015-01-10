%
% motaR.m
% des:    create mota object 'motaR' from mota-output struct.
%         the object can be investigated with the two functions 'print' and 'plot'
% usage:  m=motaR(out)
% author: Stefan Hengl 
% year:   2007
%
function  m=motaR(out)

    m.K=out.X;
    m.S=out.S;
    m.r2=out.rSquared;
    
    m=class(m,'motaR');
end