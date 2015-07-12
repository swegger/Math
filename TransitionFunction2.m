function T = TransitionFunction2(x1,x2,deterministic_function,probabilistic_function,boundry)
%% TransitionFunction
%
%   T = TransitionFunction2(x1,x2,deterministic_function,probabilistic_function,boundry)
%
%   For creating transition matrices.
%
%%

switch boundry.type
    case 'none'
        T = probabilistic_function(deterministic_function(x1,x2),x1,x2);
        
    case 'absorbing'
        T = probabilistic_function(deterministic_function(x1,x2),x1,x2);
        T(x2 < boundry.bound(1)) = 0;
        T(x2 > boundry.bound(2)) = 0;
        
    case 'reflecting'
        T = probabilistic_function(deterministic_function(x1,x2),x1,x2);
        error('reflecting boundary not yet supported!')
end