function check_integrity(s,alpha,beta,epsilon,str)
% =========================================================================
% INPUTS:
%
% str ................ name of function (string)
%
% s, alpha, beta ..... hyperparameters of BACKTRACKING algorithm
%
% epsilon ............ tolerance parameter 
%
% Functions performs a check to ensure no logical mistakes were made during
% initialization of hyperparameters. Correct values of hyperparameters
% ensure sensible value of iterations of the algorithm.
%
% Hyperparameters S,ALPHA,BETA must be set to 'none' if exact 
% line search is used.
%
%
% An upper limit of ul = 0.01 for the epsilon parameter has been chosen
% arbitrarily. If one wants to run the algorithms for error tolerance
% parameter epsilon greater than 0.01, change the variable [ ul ] to the 
% desired upper limit value.
% =========================================================================

ul = 10^-2;
err_txt = "Logical error using function:" + str +" ... Check the parameters!";

if (~strcmp(s,'none'))
    if(s <= 0)
        error(err_txt)
    end
end

if (~strcmp(alpha,'none'))
    if(alpha > 0.5 || alpha < 0)
        error(err_txt)
    end
end

if (~strcmp(beta,'none'))
    if(beta > 1 || beta < 0)
        error(err_txt)
    end
end


if ( epsilon <= 0 || epsilon > ul )
    error(err_txt)
end

end

