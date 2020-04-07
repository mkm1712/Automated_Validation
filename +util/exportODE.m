function [commandLine] = exportODE(speciesNames,paramList,ODElist)
% exportODE
% This function exports the system of differential equations as a
% standalone .m file.
% Input: paramsList = parameter list; speciesNames = list of speciesNames;
% ODElist = list of ODEs
% Output: commandLine = cell that holds everything
% Version 0.08a, 08/31/2011 by JJS

%% Write the ODEfun
commandLine{1} = sprintf('function dydt=ODEfun(t,y,params)');
commandLine{end+1} = sprintf('%% Assign names for parameters');
commandLine{end+1} = sprintf('[rpar,tau,ymax,speciesNames]=params{:};');
for i = 1:length(ODElist)
    commandLine{end+1} = sprintf(ODElist{i});
end
 
%% write utility functions
commandLine{end+1} = sprintf('\n%% utility functions');
commandLine{end+1} = sprintf('function fact = act(x,rpar)');
commandLine{end+1} = sprintf('%% hill activation function with parameters w (weight), n (Hill coeff), EC50');
commandLine{end+1} = sprintf('    w = rpar(1);');
commandLine{end+1} = sprintf('    n = rpar(2);');
commandLine{end+1} = sprintf('    EC50 = rpar(3);');
commandLine{end+1} = sprintf('    beta = (EC50.^n - 1)./(2*EC50.^n - 1);');
commandLine{end+1} = sprintf('    K = (beta - 1).^(1./n);');
commandLine{end+1} = sprintf('    fact = w.*(beta.*x.^n)./(K.^n + x.^n);');
commandLine{end+1} = sprintf('    if fact>w,                 %% cap fact(x)<= 1');
commandLine{end+1} = sprintf('        fact = w;');
commandLine{end+1} = sprintf('    end\n');

commandLine{end+1} = sprintf('function finhib = inhib(x,rpar)');
commandLine{end+1} = sprintf('%% inverse hill function with parameters w (weight), n (Hill coeff), EC50');
commandLine{end+1} = sprintf('    finhib = rpar(1) - act(x,rpar);\n');

commandLine{end+1} = sprintf('function z = OR(x,y)');
commandLine{end+1} = sprintf('%% OR logic gate');
commandLine{end+1} = sprintf('    z = x + y - x*y;\n');

commandLine{end+1} = sprintf('function z = AND(rpar,varargin)');
commandLine{end+1} = sprintf('%% AND logic gate, multiplying all of the reactants together');
commandLine{end+1} = sprintf('    w = rpar(1);');
commandLine{end+1} = sprintf('    if w == 0,');
commandLine{end+1} = sprintf('        z = 0;'); 
commandLine{end+1} = sprintf('    else');
commandLine{end+1} = sprintf('        v = cell2mat(varargin);'); 
commandLine{end+1} = sprintf('        z = prod(v)/w^(nargin-2); % need to divide by w^(#reactants-1) to eliminate the extra ws'); 
commandLine{end+1} = sprintf('    end');

commandLine = commandLine';