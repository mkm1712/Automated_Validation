function [percentMatch, resultChart, BMatch, byClass] = Automated_Validation_V1(Modelfile_path, Validationfile_path, Int_time, Steady_time,Threshold, Model_version)
% Calculates the percent agreement and writes resultChart variable to
% "Validation Results.xlsx".
%   inputs:
%   - Modelfile_path = model file path (.xlsx)
%   - Validationfile_path = validation file path (.xlsx)
%   - Int_time = Initial simulation time
%   - Steady_time = Main simulation time (to steady state)
%   - Threshold = Threshold for determination of increase or decrease (%)
%   - Model_version = Basic Version (1), New Vesrion (2)
%
%   outputs:
%   - percentMatch = percent agreement
%   - resultChart = chart containing results of each individual validation
%     simulation.
%   - BMatch = boolean vector containing the result of each validation (1 =
%     correct, 0 = incorrect)
%   - byClass = percent match by validation class

% Part 1
% Delete the previously formed ODE to make sure it's rewritten
currentfolder= pwd;
if exist([currentfolder '\ODEfun.m'],'file') == 2
    delete('ODEfun.m');
end
%
% Part 2
% Parse out model's name (need by netflux)
namepos = strfind(Modelfile_path,'.xls');
if rem (namepos,1)== 0
    namestr = Modelfile_path(1:namepos-1);
    namestr = cellstr(namestr);
else
    disp ('Error: please insert the path for model file with correct extension');
    return
end
%
% Part 3
% generate ODE from model spreadsheet
[specID,reactionIDs,~,paramList,ODElist,~,~] = util.xls2Netflux(namestr,Modelfile_path);
if Model_version == 1
commandLine = util.exportODE(specID,paramList,ODElist); %Model Version 1
else
commandLine = util.exportODE2(specID,paramList,ODElist); %Model Version 2
end
util.textwrite('ODEfun.m',commandLine);
%
% Part 4
% set up simulation options
tspan_ss = [0 Steady_time]; % run out to ss
tspan_int = [0 Int_time]; % run out to ss
%
% Part 5
% Read the validation sheet & Identification of each main columns
[~, txt, ~] = xlsread(Validationfile_path);
Input1Column = cellfun(@(x)isequal(x,'Input'), txt(1,:))';
Input2Column = cellfun(@(x)isequal(x,'Input 2'), txt(1,:))';
IncodeColumn = cellfun(@(x)isequal(x,'Input Code'), txt(1,:))';
OutputColumn = cellfun(@(x)isequal(x,'Output'), txt(1,:))';
MeasurementColumn = cellfun(@(x)isequal(x,'Measurement'), txt(1,:))';
IDColumn =  cellfun(@(x)isequal(x,'ID'), txt(1,:))';
ValIDColumn = cellfun(@(x)isequal(x,'Validation tag'), txt(1,:))';
controlColumn = cellfun(@(x)isequal(x,'Control'), txt(1,:))';
%
% Part 6
% remove rows without data
noData = cellfun(@(x)isequal(x,'No Data'), txt(1:end, MeasurementColumn));
txt(noData, :) = [];
noData = cellfun(@isempty, txt(1:end, MeasurementColumn));
txt(noData, :) = [];
%
% Part 7
% Having validation inputs in workspace for debugging & Extracting Data from Validation file
assignin('base', 'txt', txt);
validationIDs = txt(2:end, IDColumn);
input1 = txt(2:end, Input1Column);
input2 = txt(2:end, Input2Column);
inputCode = txt(2:end, IncodeColumn);
measurement = txt(2:end,MeasurementColumn);
outputSpec = txt(2:end, OutputColumn);
validationTags = txt(2:end, ValIDColumn);
control = txt(2:end, controlColumn); % when we have predefined control state
UninputCode = unique(inputCode);
%
% Part 8
% convert species and rxn names to integer values to map name and reaction
% ID's to numerical integers, allowing the input code to be evaluated
% directly
for k = 1:length(specID)
    if isempty(specID{k})
        disp (['specID ',num2str(k),' missing']);
    else
        eval([specID{k},' = ',num2str(k),';']);
    end
end
for i = 1:length(reactionIDs)
    if isempty(reactionIDs{i})
        disp (['reactionIDs ',num2str(i),' missing']);
    else
        eval([reactionIDs{i},' = ',num2str(i),';']);
    end
end
for j = 1:length(validationIDs)
    if isempty(validationIDs{j})
        disp (['validationIDs ',num2str(j),' missing']);
    else
        eval([validationIDs{j}, ' = ', num2str(j), ';']);
    end
end
%
% Part 9
% Set validation threshold change
thresh1 = Threshold/100;
thresh2 = 1e-3*thresh1; % Not considering very low changes in abs value (maybe due numerical errors) 
options = odeset('RelTol',0.1*thresh1,'AbsTol',0.1*thresh2); % Specifying computational accuracy based on threshold
inc = {'Increase'};
dec = {'Decrease'};
noc = {'No Change'};
numMatching = 0; % number of predictions consistent with the qualitative literature species behavior
%
% Part 10
% Define the size of some variable changing in the loops and find indices of output species
outputSpeciesIndex = zeros(1, length(measurement));
yStartL = cell(1, length(UninputCode));
yEndL = cell(1, length(UninputCode));
prediction = cell(1, length(inputCode));
predChange = cell(1, length(inputCode));
match = zeros(1, length(measurement));
c =zeros (length(UninputCode),1);
for k = 1:length(outputSpec)
    [~,outputSpeciesIndex(k)] = ismember(outputSpec{k},specID);
end
%
% Part 11
% loop over all validation simulations read from the excel sheet
for i = 1:length(UninputCode)
    disp(['Simulation # ', num2str(i), ' of ',num2str(length(UninputCode))]) % write the simulation number to the command line to track loop progress
    [w,n,EC50,tau,ymax,y0] = paramList{:}; % reset params
    % Initial (control) Simulation
    a= strfind (UninputCode{i},';');
    b= strfind (UninputCode{i},'w(');
    if length(a)> 1.1 && length(b)> 0.1
        eval([UninputCode{i}(1:a(1)), '%', UninputCode{i}(a(1)+1:end)]);
        c(i)=1;
    end
    rpar = [w;n;EC50];
    params = {rpar,tau,ymax,specID};
    [~,y] = ode15s(@ODEfun, tspan_int, y0, options, params);
    yStart = y(end,:)'; % use the "no input" steady state as control
    
    % evaluate validation conditions from excel sheet
    eval(UninputCode{i});
    
    % Main Simulation
    rpar = [w;n;EC50];
    params = {rpar,tau,ymax,specID};
    [~,y] = ode15s(@ODEfun, tspan_ss, y0, options, params);
    yEnd = y(end,:)';
    
    % Determine Change of Species' Activity after Stimulation
    yStartL{i} = yStart;
    yEndL{i} = yEnd;
end
% Part 12
% Determination of activity change in each experiment
for i=1:length(inputCode)
    idx = find(ismember(UninputCode, inputCode{i}));
    %     if isempty(control{idx}) % if control validation defined
    activityChange = (real(yEndL{idx}(outputSpeciesIndex(i)))-real(yStartL{idx}(outputSpeciesIndex(i))));
    activityChange_check = abs(((real(yEndL{idx}(outputSpeciesIndex(i)))-real(yStartL{idx}(outputSpeciesIndex(i))))/real(yStartL{idx}(outputSpeciesIndex(i)))));
    %     else
    %         indexControl = eval(control{idx});
    %         assignin('base', 'indexControl',indexControl);
    %         activityChange = real(yEndL{idx}(outputSpeciesIndex(i))) - real(yStartL{indexControl}(outputSpeciesIndex(i)));
    %     end
    
    % Determine type of Changes
    if activityChange > thresh2 && activityChange_check > thresh1 % increase
        prediction{i} = 'Increase';
        predChange{i} = num2str(activityChange);
        if isequal(inc,measurement(i))
            numMatching = numMatching + 1;
            match(i) = 1; %if the simulation matches the experimental validation put a 1 in the vector
        else
            match(i) = 0; %if the simulation does not match put a 0 in the matrix
        end
    elseif activityChange < -thresh2 && activityChange_check > thresh1 % decrease
        prediction{i} = 'Decrease';
        predChange{i} = num2str(activityChange);
        if isequal(dec,measurement(i))
            numMatching = numMatching + 1;
            match(i) = 1;
        else
            match(i) = 0;
        end
    else % no change
        prediction{i} = 'No Change';
        predChange{i} = num2str(activityChange);
        if isequal(noc,measurement(i))
            numMatching = numMatching + 1;
            match(i) = 1;
        else
            match(i) = 0;
        end
    end
    
end
%
% Part 12
% Produce the experiment-simulation match result and define the class of each experiment
match2 = cell (1,length(match));
for j = 1:length(match)
    if match(j) == 1
        match2{j} = 'yes';
    else
        match2{j} = 'no';
    end
end

BMatch = match';
ind = find(strcmp('Inp-out', validationTags));
byClass.Inpout = sum(BMatch(ind))/length(BMatch(ind))*100;

ind = find(strcmp('Inp-int', validationTags));
byClass.Inpint = sum(BMatch(ind))/length(BMatch(ind))*100;

ind = find(strcmp('Int-over', validationTags));
byClass.Intover = sum(BMatch(ind))/length(BMatch(ind))*100;

ind = find(strcmp('Int-inh', validationTags));
byClass.Intinh = sum(BMatch(ind))/length(BMatch(ind))*100;
%
% Part 13
% produe the outputs
resultChart = {validationIDs, input1, input2, outputSpec, measurement, prediction', predChange', match2', validationTags}; %create a cell array showing input, output, and whether they matched validation
header = {'ID', 'input' ,'input 2', 'output', 'measurement', 'prediction', 'predicted change', 'match', 'tag'};
resultChart = horzcat(resultChart{:});
resultChart = vertcat(header, resultChart);
% Write the results in a file
delete ('Validation_Results.xlsx');
xlswrite('Validation_Results.xlsx', resultChart);
disp(['wrote ', 'Validation_Results.xlsx']);
percentMatch = numMatching/length(measurement)*100;
delete('ODEfun.m');