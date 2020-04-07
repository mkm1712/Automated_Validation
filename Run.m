Modelfile_path ='C:\Users\Ali\Desktop\AutomatedV\Toy_model.xlsx'; % Change the path to the file current path
Validationfile_path ='C:\Users\Ali\Desktop\AutomatedV\Toy_model_exp.xlsx';% Change the path to the file current path
Int_time = 40;
Steady_time = 40;
Threshold = 0.1; % percent of changes
Model_version = 1; % 1= Original 2= Modified
[percentMatch, resultChart, BMatch, byClass] = Automated_Validation_V1(Modelfile_path, Validationfile_path, Int_time, Steady_time,Threshold, Model_version);