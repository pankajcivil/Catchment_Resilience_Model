classdef CatchmentResilienceModel <  handle
% CATCHMENTRESILIENCEMODEL Class definition for building a hillslope vadose zone model having two 'steady states'
%
% Description
%   This class allows the building and time-integration of the coupled
%   hillslope-vadose zone Boussinesq model from Peterson & Western (2012).
%   This model contains a positive feedback resulting from the reduction of
%   Leaf Area Index (LAI) when a saline water table intersects the root
%   zone. This positive feedback can result in two steady state water table
%   elevations throughout the catchment for a single parameter set (henceforth, 
%   the steady states are called 'attractors'). Using this code, the model
%   can be bult and solved using a wide range of climate forcing
%   time-steps. By linking with the companion class definition,
%   'LimitCycleContinuation.m' the depth to water table of the attractors 
%   and threshold (henceforth,referred to as the repellor) can be
%   quantified with a change in a single model parameter.
%
%   The above mentioned positive feedback can be turn on or off by setting
%   the input 'model.modelParameters.aquifer.doPositiveFeedback' to true.
%   If set to false, then only one attractor can emerge.
%
%   For an example of how to build and run a model see the example below.
%
%   For more details of the model see Peterson et al. (2009). For examples
%   of using this model to explore the attractors under daily forcing see
%   Peteron, Western and Argent (2012a). For concepts and methods in
%   hydrological resilience see Peterson, Western and Argent (2012b)
%   
% Example: 
%   The following example solves the model using monthly average forcing:
%   -----------------------------------------------------------------------
%   % Get model parameters used within Peterson and Western 2012:
%   modelParams = CatchmentResilienceModel.getModelParametersStructure();
%
%   % Create the model object:
%   model = CatchmentResilienceModel(modelParams);
% 
%   % Set the initial conditions as defined by a fraction of the maximum
%   % saturated thickness of 0.9 and the soil porosity of 0.4:
%   setInitialConditions(model, [0.9, 0.4],'fractional');
%
%   % Set differential equation solve options to the defaults:
%   setSolverOptions(model);
%
%   % Set the climate forcing data to the monthly average with a spline fit to
%   % between the end of each month. For addition forcing options see the
%   % help for 'setClimateData' below.
%   setClimateData(model, 3, true, climate_daily,[]);
%
%   % Solve the model from year 1890 to 2000 and output 12 evenly spaced 
%   results per year of simulation:
%   solveModel(model, 1890, 2000, 12);
%
%   % Plot a cross section heads from 1890 to 2000 at 10 years increments:
%   plotCrossSections(model, [1890: 10 : 2000], true, false);
%
%   % Get a list of possible fluxes and variables to plot over time.
%   fluxesAndVariables = plotTimeSeries(model);
%
%   % Plot time series of depth to water table and Leaf Area Index (LAI)
%   over time from 1890 to 2000 and return flux data in variable 'fluxData'.
%   fluxData = plotTimeSeries(model, fluxesAndVariables( [19; 27] ), 1890, 2000, [250; 1000])
%
%   The following example solves the model using the up-scaled stochastic forcing
%   from Peterson and Western (2012) at two initial conditions. For the
%   two initial conditions:
%   -----------------------------------------------------------------------
%   % Get model parameters used within Peterson and Western 2012:
%   modelParams = CatchmentResilienceModel.getModelParametersStructure();
%
%   % Create two model objects, one with a deep initial condition and one 
%   % with a shallow initial condition:
%   model_deep = CatchmentResilienceModel(modelParams);
%   model_shal = CatchmentResilienceModel(modelParams);
% 
%   % Set the initial conditions to 20% and 90% of the saturated thickness.
%   setInitialConditions(model_deep, [0.2, 0.4],'fractional');
%   setInitialConditions(model_shal, [0.9, 0.4],'fractional');
%
%   % For both models, set differential equation solve options to the defaults:
%   setSolverOptions(model_deep);
%   setSolverOptions(model_shal);
%
%   % Set the climate forcing data to the up-scaled daily data (as detailed 
%   within Peterson and Western 2012) using 50 percentiles.
%   setClimateData(model_deep, 4, true, climate_daily, 50);
%   setClimateData(model_shal, 4, true, climate_daily, 50);
%
%   % Solve both models from year 1890 to 2000 and output 12 evenly spaced 
%   results per year of simulation:
%   solveModel(model_deep, 1890, 2000, 12);
%   solveModel(model_shal, 1890, 2000, 12);
%
%   % Get a list of possible fluxes and variables to plot over time.
%   fluxesAndVariables = plotTimeSeries(model_deep);
%
%   % Plot time series of depth to water table and soil moisture (theta)
%   over time from 1890 to 2000 and return flux data in variable 'fluxData'.
%   fluxData = plotTimeSeries(model_deep, fluxesAndVariables( [19; 23] ), 1890, 2000, [250; 1000])
%   fluxData = plotTimeSeries(model_shal, fluxesAndVariables( [19; 23] ),
%   1890, 2000, [250; 1000])
%
% See also
%   getModelParametersStructure: default_model_parameters;
%   CatchmentResilienceModel: model_construction;
%   setInitialConditions: set_initial_conditions_for_state_variables;
%   setSolverOptions: set_differential_equation_solve_options;
%   setClimateData: set_the_climate_data;
%   solveModel: time_integration_solution_of_model;
%   plotTimeSeries: plot_time_series_of_solution_outputs;
%   plotCrossSections: plot_cross_sections_of_solutions;
%   plotAttractors: estimate_attractors_and_plot;
%
% References:
%   Peterson, T. J., R. M. Argent, A. W. Western, and F. H. S. Chiew
%   (2009), Multiple stable states in hydrological models: An ecohydrological
%   investigation, Water Resour. Res., 45, W03406, doi:10.1029/2008WR006886. 
%
%   Peterson, T. J. and A. W. Western (2012), On the 
%   existence and emergence of multiple hydrological attractors under 
%   stochastic daily forcing: 1. Identifying the attractors, Water 
%   Resour. Res., (Submitted August 2012)
%
%   Peterson, T. J., A. W. Western, and R. M. Argent (2012a), On the 
%   existence and emergence of multiple hydrological attractors under 
%   stochastic daily forcing: 2: Exploring the switching between
%   attractors, Water Resour. Res., (Submitted August 2012)
%
%   Peterson, T. J., A. W. Western, and R. M. Argent (2012b), Analytical 
%   methods for ecosystem resilience: A hydrological investigation, Water 
%   Resour. Res., (Accepted August 2012)
%   
% Author: 
%   Dr. Tim Peterson
%   The Department of Infrastructure
%   The University of Melbourne
%   Parkville, VIC, AUSTRALIA
%   ph: +613 8344 9950
%   email: timjp@unimelb.edu.au
%
% License:
%    Copyright (C) 2012  Dr. Tim J. Peterson, Prof. Andrew W. Western and Dr. Robert M. Argent.
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>
%
   properties
       
       % Model parameters
       modelParameters;
       
       % climate forcing
       climateData;
       
       % State variable initial conditions
       initialConditions
       
       % PDE solver 
       solverOptions; 
       
       % Time integration results
       results;
       
    end 
    
    properties (SetAccess=protected, GetAccess=protected)
       % Derived constants and state variable dependent variables.
       geometryVariables;
       variables;
       
       % Model state variables
       stateVariables =  [];                     
       
       % Model public fluxes
       fluxes = [];
   end
   
   methods(Static)
%% Default model parameters       
       function modelParameterExample = getModelParametersStructure()
% GETMODELPARAMETERSSTRUCTURE Get the default model parameters
%
% Syntax:
%   modelParameterExample = getModelParametersStructure()
%
% Description:
%   Static function that returns the model parameters (and parameter
%   structure) used within Peterson and Western (2012) and Peterson,
%   Western and Argent (2012).
%
% Inputs:
%   (none)
%
% Outputs:
%   modelParameterExample - structure variable containing the model
%   parameters and model geometry.
%
% Dependencies
%   (none)
%
% Author: 
%   Dr. Tim Peterson, The Department of Infrastructure
%   Engineering, The University of Melbourne.
%
% Date:
%   2 August 2012

           % Model geometry data:
           u_nodes =  (0: 10: 2000)';
           nNodes = size(u_nodes,1);
           modelParameterExample.geometry.u_nodes =u_nodes;
           
           modelParameterExample.geometry.width = 2*250 * exp(0.0004 * u_nodes);
           modelParameterExample.geometry.basementElevation = 1.75e-5 * u_nodes.^2 + 0.003 * u_nodes;
           modelParameterExample.geometry.surfaceElevation = 1e-5 * u_nodes.^2 + 50;
           modelParameterExample.geometry.soilDepth = 0.05 * (modelParameterExample.geometry.surfaceElevation - modelParameterExample.geometry.basementElevation);
           modelParameterExample.geometry.streamDepth = zeros(nNodes,1);

           % Aquifer parameters.
           % The first line assigns all nodes to have unit type number one.
           modelParameterExample.aquifer.lambda_m = 0.1;
           modelParameterExample.aquifer.singularityFraction = 1.0;
           modelParameterExample.aquifer.doPositiveFeedback = true;     % Turn on the positive feedback!
           modelParameterExample.aquifer.UnitCoverage = [modelParameterExample.geometry.u_nodes, ones(nNodes,1)];
           modelParameterExample.aquifer.Units{1,1}.k = 1.0;
           modelParameterExample.aquifer.Units{1,1}.f = 0.05;
           modelParameterExample.aquifer.Units{1,1}.d_khalf = 2;
           modelParameterExample.aquifer.Units{1,1}.tor = 4;

           % Vadose zone parameters
           modelParameterExample.vadose.UnitCoverage = [modelParameterExample.geometry.u_nodes, ones(nNodes,1)];
           modelParameterExample.vadose.Units{1,1}.theta_s = 0.43;
           modelParameterExample.vadose.Units{1,1}.theta_r = 0.109;
           modelParameterExample.vadose.Units{1,1}.BC_phi = 0.168;
           modelParameterExample.vadose.Units{1,1}.BC_psi_a = -291.7;
           modelParameterExample.vadose.Units{1,1}.k_vert = 28.8;
           modelParameterExample.vadose.Units{1,1}.d_evap = 0.5;
           modelParameterExample.vadose.Units{1,1}.Io= 200;

           % Vegetation cover parameters for grasslands
           modelParameterExample.veg.UnitCoverage = [NaN, 0; modelParameterExample.geometry.u_nodes, ones(nNodes,1)];
           modelParameterExample.veg.Units{1,1}.F_canopy = 1.0;
           modelParameterExample.veg.Units{1,1}.k_light = 0.6;
           modelParameterExample.veg.Units{1,1}.wiltingPoint_MPa = -3.0;
           modelParameterExample.veg.Units{1,1}.stomataClosure_MPa = -0.03;
           modelParameterExample.veg.Units{1,1}.d_0p5LAI = 2;
           modelParameterExample.veg.Units{1,1}.alpha = 3;           
           modelParameterExample.veg.Units{1,1}.LAI_avg = {'sinLAI'; 0.8064; 3.205; -0.08111};
           
           % Stream parameters
           modelParameterExample.stream.river_conductScaler = 0;
           
           % Set boundary conditions
           modelParameterExample.boundaryConditions.Lower.Type = 'GHB';
           modelParameterExample.boundaryConditions.Lower.Param1 = 0.2;
           modelParameterExample.boundaryConditions.Lower.Param2 = 48.0;          
           
           modelParameterExample.boundaryConditions.Upper.Type = 'Q';
           modelParameterExample.boundaryConditions.Upper.Param1 = 0;
           modelParameterExample.boundaryConditions.Upper.Param2 = 0;                   
           
       end
   end
   
   methods 

%% Construct the model
       function obj=CatchmentResilienceModel(modelInputs)
% CATCHMENTRESILIENCEMODEL Model construction.
%
% Syntax:
%   obj = CatchmentResilienceModel(modelInputs)
%
% Description:
%   Constructs an instance of the model object. This can contain only the
%   model parameters or, if the input is an object of the same type, then a
%   duplicate model is created to that of the input.
%
%   Once the model object is created, the following steps should be followed
%   to run the model:
%   1. Set the state variable initial conditions.
%   2. Set the differential equation solver options;
%   3. Set the climate forcing;
%   4. Solve the model over a user defined time range;
%
%   Once the model is successfully solved, the results can be plotted as a
%   time series at any model node or as cross sections of groundwater level 
%   and soil moisture. The successfully solved model can also be used to
%   estimate the attractors at over a defined parameter range or to start
%   the continuation analysis. In the example section below, details are
%   given for building and solving a model and plotting the results.
%   
% Inputs:
%   modelInputs - an input structure containing all model parameters or an 
%   object of type 'CatchmentResilienceModel' 
%
% Outputs:
%   obj - an output object of type 'CatchmentResilienceModel'. If the input
%   was also an object of this type, then the output will be a duplicate of
%   it. 
%
% See also:
%   CatchmentResilienceModel: class_description;
%
% Dependencies
%   (none)
%
% Author: 
%   Dr. Tim Peterson, The Department of Infrastructure
%   Engineering, The University of Melbourne.
%
% Date:
%   2 August 2012

           if isa(modelInputs,'CatchmentResilienceModel')
                obj.modelParameters = modelInputs.modelParameters; 
                obj.climateData = modelInputs.climateData; 
                obj.initialConditions = modelInputs.initialConditions;                 
                obj.solverOptions = modelInputs.solverOptions;     
                obj.results = modelInputs.results;                     
           else
               obj.modelParameters = modelInputs;                 
           end            
           
           % Calc. model constants
           setGeometryVariables(obj);
       end
       
%% Set climate forcinf data
       function setClimateData(obj, userOption, isAtDailyTimestep, inputClimateData, npcntiles)
% SETCLIMATEDATA Set the model climate data
%
% Syntax:
%   setClimateData(obj, userOption, isAtDailyTimestep, inputClimateData, npcntiles)
%
% Description:    
%   Sets the climate data to one of four options. The options range from
%   the simulation under stochastic forcing (daily or monthly) to daily
%   data up-scaled to a monthly time step using fixed or time varying probability
%   functions. For details of the climate scaling see Peterson & Western
%   2012. Below are details of the options:
% 
%   1. Simulation to use the input climate data. This data will be
%      converted to units of mm/day (assumed to be in mm over time step). The
%      option 'scalingPercentiles' will be ignored. This option was used to
%      produce the stochastic simulations within Peterson, Western and
%      Argent 2012.
%
%   2. Simulation to use the input climate data scaled up from
%      daily to monthly. For each month of each year the
%      percentiles within 'scalingPercentiles' will be calculated. This user
%      option requires input of daily data. This option can be used to check
%      the up-scaling method against simulations using option 1. 
%      
%   3. Simulations to use smooth sigmoidal repeated functions of
%      average climate forcing. The input 'inputClimateData' must be: 
%      (i) daily data that will be upscaled to monthly and a sigmoidal 
%      function fit to the monthly mean forcing; or (ii) monthly data in which
%      a sigmoidal function will fit to the monthly mean. For both options 
%      the input 'scalingPercentiles' will be ignored. This type of 
%      simulation was included for limit-cycle continuation without up-scaling
%      (as undertaken in Peterson Argent Western and Chiew 2009).
%
%   4. Simulations to use smooth sigmoidal repeated functions of
%      climate forcing BUT at various user defined percentiles.
%      This option requires inputClimateData to of daily forcing.
%      It was included to allow limit-cycle continuation with
%      inclusion of sub-monthly climate variability. This option was used to
%      quantify the attractors and repellor in Peterson & Western 2012.
%
% Inputs:
%   obj - object of type 'CatchmentResilienceModel'.
%
%   userOption - scalar integer from 1 to 4 used to set the type of climate
%   forcing (see description above);
%
%   isAtDailyTimestep -  logical variable defining the time step of the
%   input climate data. A true value indicates daily time step, while false
%   indicates a monthly time step.
%
%   inputClimateData - nx5 matrix of climate data. The columns of data must
%   in the following order [year, month, day, precipitation, potential
%   areal evaporation]. If 'isAtDailyTimestep' is false then the daily
%   column should be be the last day of each month.
%
%   npcntiles - the number of percentiles to use in the up-scaling. The
%   percentile are from zero to one and are distributed as follows:
%   [1-logspace(3,0,npcntiles )./10^3, 1];
%
% Outputs:
%   (model object, obj, is updated to contain the climate forcing.)
%
% See also:
%   CatchmentResilienceModel: class_description;
%
% Dependencies
%   (none)
%
% Author: 
%   Dr. Tim Peterson, The Department of Infrastructure
%   Engineering, The University of Melbourne.
%
% Date:
%   2 August 2012

           
           % Define climate scaling
           scalingPercentiles = [1-logspace(3,0,npcntiles )./10^3, 1];
           
           % Check input format of the data.
           %----------------------
           if ~ isnumeric(inputClimateData); 
               error('Input "inputClimateData" must be data of contain five columns: year, month, day, rainfall from prior time point, FAO56 ET from prior time point.');
           end
           
           % Check format of input matrix   
           if size(inputClimateData,2) ~= 5;
               error('Input climate data must contain five columns: year, month, day, rainfall from prior time point, FAO56 ET from prior time point.');
           elseif max(inputClimateData(:,2))>12 || max(inputClimateData(:,3))>31
               error('Input climate data error. The second column must be the month and third colum the day of the month.')
           elseif (userOption==2 || userOption==4)  && (nargin==3 || isempty(scalingPercentiles))
               error('User options 2 and 4 require input of percentiles for scaling but "scalingPercentiles" is empty or not input.')
           end    

           % If calculating percentiles, check the input percentiles. 
           if userOption==2 || userOption==4 || userOption==5
               if isempty(scalingPercentiles)
                   error('"ScalingPercentiles" must be input for user options 2, 4 and 5.');
               elseif ~isvector(scalingPercentiles) || length(scalingPercentiles) < 3;
                   error('"ScalingPercentiles" must be a vector of at least three percentiles.');
               elseif max(scalingPercentiles) > 1 || min(scalingPercentiles) < 0
                   error('"ScalingPercentiles" must contain percentiles >=0 and <=1 .');
               %elseif max( diff(scalingPercentiles)) ~= max( diff(scalingPercentiles))
               %    error('The percentiles must be of a uniform step size! (limitation of extended Simpsons composite integration method).');                   
               elseif mod(length(scalingPercentiles),2) == 0
                   error('The number of percentiles must be of an uneven number! (limitation of extended Simpsons composite integration method).');                      
               end
               display('... Simulations will be conducted using down-scaled climate data at the user-supplied percentiles');
               npcntiles = length(scalingPercentiles);
           end
               
           % Check if there is any missing data or inconsistent format. 
           for i=2:size(inputClimateData,1)

               % Define expected date.
               if isAtDailyTimestep
                   expectedDate =datenum( inputClimateData(i-1,1), inputClimateData(i-1,2), inputClimateData(i-1,3)) + 1;
                   tolerance = 0.5;     %set tolerance at 1/2 day.
               else
                   if inputClimateData(i-1,2)==12;
                        expectedDate =datenum( inputClimateData(i-1,1)+1, 1, inputClimateData(i-1,3));
                   else
                        expectedDate =datenum( inputClimateData(i-1,1), inputClimateData(i-1,2)+1, inputClimateData(i-1,3));
                        tolerance = 1;     %set tolerance at 1 day.
                   end                                              
               end

               % Check expected data against current data.
               if abs(datenum( inputClimateData(i,1), inputClimateData(i,2), inputClimateData(i,3)) - expectedDate) > tolerance
                   if isAtDailyTimestep
                        error(['The input data does not appear to be of a consistent time step. Please check the data at row ',num2str(i),'.']);
                   else
                        error(['The input data does not appear to be of a consistent time step. Note, for monthly timetsep simulations, the input data must be monthly. Please check the data at row ',num2str(i),'.']); 
                   end
               end
           end
           %----------------------
           
           % Input ORIGINAL variable step climate data.
           obj.climateData.input = inputClimateData;
           
           % If "inputClimateData" is simply to be converted to mm/day and
           % input to object, then do this operation and return.
           %----------------------
           if userOption==1                           
              % Aggregate year, month and day to a single number of year+fraction of year.
              t_data = datenum(inputClimateData(:,1),inputClimateData(:,2),inputClimateData(:,3)) - datenum(inputClimateData(:,1),1,1) + 1;
              days_per_year = eomday(inputClimateData(:,1), 2) + 365-28;
              t_data = inputClimateData(:,1) + t_data./days_per_year;              

              % Convert data series to daily fluxes.
              diff_days=[];
              diff_days(2:size(t_data,1),1) = days_per_year(2:end).*(t_data(2:end,1) - t_data(1:end-1,1));
              diff_days(1,1)=diff_days(2,1);
              if any(diff_days>1+eps) || any(diff_days+eps<1-eps);
                  inputClimateData(:,4) = inputClimateData(:,4)./diff_days;    
                  inputClimateData(:,5) = inputClimateData(:,5)./diff_days;    
              end  
              
              % Climate data in the required model format to the object.
              obj.climateData.type = 'data';
              obj.climateData.doClimateScaling = false;    
              obj.climateData.scalingPercentiles = [];
              obj.climateData.data.P = [t_data, inputClimateData(:,4)];          
              obj.climateData.data.ET = [t_data, inputClimateData(:,5)];          
              obj.climateData.functionParameters_P = [];
              obj.climateData.functionParameters_ET = []; 
              return;               
           end                      
           %----------------------
           
           % Aggregate data to monthly matrix if userOption is to use
           % monthly data. If the data is already monthly then this
           % operation has no effect.
           %----------------------
           
           % Check that all months are complete, ie no partial months.
           nyears = inputClimateData(end,1) - inputClimateData(1,1) +1;
           monthly_countDays = cell( nyears,12);
           for i=1: size(inputClimateData,1)                
                imonth = inputClimateData(i,2);                
                irow = min(nyears, inputClimateData(i,1) - inputClimateData(1,1)+1);                   
                if isAtDailyTimestep
                    diff_days = 1;
                else
                    diff_days = inputClimateData(i,3) - inputClimateData(i-1,3);
                end                 
                monthly_countDays{irow,imonth} = [monthly_countDays{irow,imonth}; diff_days];
           end
           for irow=1:nyears
               for imonth=1:max(inputClimateData(:,2) )
                   % Check that a full month has been input
                   if sum(monthly_countDays{irow,imonth}) ~= eomday(irow - 1 + inputClimateData(1,1), imonth)
                       error('Input climate data cannot contain incomplete months of data!')
                   end    
               end
           end
           
           % Only monthly percentile or monthly average forcing data is required.            
           if userOption~=2 
               nyears = 1;
           end
           
           % Format data into a matrix of years versus months
           monthly_P = cell( nyears,12);
           monthly_ET = monthly_P;
           monthly_countDays = monthly_P;
           monthly_countDays{1 ,inputClimateData(1,2)} = 1;
           monthly_P{1 ,inputClimateData(1,2)} = inputClimateData(1,4);
           monthly_ET{1 ,inputClimateData(1,2)} = inputClimateData(1,5);                                            
           for i=2: size(inputClimateData,1)
                iyear = inputClimateData(i,1);
                imonth = inputClimateData(i,2);                
                irow = min(nyears, iyear - inputClimateData(1,1)+1);
                monthly_P{irow,imonth} = [monthly_P{irow,imonth}; inputClimateData(i,4)];
                monthly_ET{irow,imonth} = [monthly_ET{irow,imonth}; inputClimateData(i,5)];
                
                 if isAtDailyTimestep
                     diff_days = 1;
                 else
                     diff_days = inputClimateData(i,3) - inputClimateData(i-1,3);
                 end
                 
                 monthly_countDays{irow,imonth} = [monthly_countDays{irow,imonth}; diff_days];
           end
           % Convert to mm/day
           for irow=1:size(monthly_P,1);
               for imonth=1:size(monthly_P,2);                                                                           
                   monthly_P{irow, imonth} = monthly_P{irow, imonth}./monthly_countDays{irow,imonth};
                   monthly_ET{irow, imonth} = monthly_ET{irow, imonth}./monthly_countDays{irow,imonth};
               end
           end
           %----------------------
           
           % Process the cell structure of monthly data according to the
           % user options.
           %-------------------------
           switch userOption
               case 2       % Monthly simulation using percentiles of daily forcing rates.
                   % Calculate percentiles of daily forcing.
                   for irow=1:size(monthly_P,1);
                       for imonth=1:size(monthly_P,2);
                             for ipcnt = 1: npcntiles
                                 if size( monthly_P{irow, imonth}, 1) > 0
                                     monthly_P_pcnt(imonth, irow, ipcnt ) =  prctile(monthly_P{irow, imonth}, 100.*scalingPercentiles(ipcnt) );
                                     monthly_ET_pcnt(imonth, irow, ipcnt ) =  prctile(monthly_ET{irow, imonth}, 100-100.*scalingPercentiles(ipcnt) );
                                 else
                                     monthly_P_pcnt(imonth, irow, ipcnt ) = NaN;
                                     monthly_ET_pcnt(imonth, irow, ipcnt ) = NaN;
                                 end
                             end
                       end
                   end  
                   
                   % Derive monthly time points of climate data and convert to vectors.
                   nyears = size(monthly_P,1);
                   nmonths = size(monthly_P,2);
                   years = repmat(  (inputClimateData(1,1):inputClimateData(end,1))', 1, 12);
                   months = repmat(  (1:nmonths), nyears, 1);
                   years = reshape(years', nyears* nmonths, 1);
                   months = reshape(months', nyears* nmonths, 1);
                   days = eomday( years, months);
                   
                   % Aggregate year, month and day to a single number of
                   % year+fraction of year.
                   t_data = datenum(years,months,days) - datenum(years,1,1) + 1;
                   days_per_year = eomday(years, 2) + 365-28;
                   t_data = years + t_data./days_per_year;
                   
                   % Add climate data in the required model format to the object.
                   obj.climateData.type = 'data';
                   obj.climateData.doClimateScaling = true;    
                   obj.climateData.scalingPercentiles = scalingPercentiles;
                   climateData_temp = [t_data, ...
                       reshape(monthly_P_pcnt, nyears* nmonths, npcntiles), ...
                       reshape(monthly_ET_pcnt, nyears* nmonths, npcntiles)];
                   
                   obj.climateData.functionParameters_P = [];
                   obj.climateData.functionParameters_ET = [];                  
                   
                   % Remove rows of NaN climate data from climateData_temp
                   climateData_temp = climateData_temp( ~isnan(climateData_temp(:,2)),:);

                   % Add climate data to two structures within model
                   % object.
                   obj.climateData.data.P = climateData_temp(:,1:2);
                   obj.climateData.data.ET = climateData_temp(:,[1,3]);
                   clear climateData_temp                  
                   
                   return;
                   
               case 3       % Monthly average simulation using a smooth sigmoidal curve.
                   % Calculate monthly mean.
                   for irow=1:size(monthly_P,1);
                       for imonth=1:size(monthly_P,2);
                             monthly_P_avg(imonth,1) =  mean(monthly_P{irow, imonth});
                             monthly_ET_avg(imonth,1) =  mean(monthly_ET{irow, imonth});
                       end
                   end
                   
                   % Add end of year flux to start of year
                   monthly_P_avg = [monthly_P_avg(end); monthly_P_avg];
                   monthly_ET_avg = [monthly_ET_avg(end); monthly_ET_avg];

                   time_points = [0; (cumsum(eomday(1999,1:12))./365)'];                   
                   obj.climateData.functionParameters_P = pchip([ time_points(end-1) - 1; time_points; 1+time_points(2)], [ monthly_P_avg(end-1,:); monthly_P_avg; monthly_P_avg(2,:) ]' );
                   obj.climateData.functionParameters_ET = pchip([ time_points(end-1) - 1; time_points; 1+time_points(2)], [ monthly_ET_avg(end-1,:); monthly_ET_avg; monthly_ET_avg(2,:) ]' );
                   
                   % Input results to model.
                   obj.climateData.type = 'function';
                   obj.climateData.doClimateScaling = false;    
                   obj.climateData.scalingPercentiles = [];
                   obj.climateData.data.P = [];
                   obj.climateData.data.ET = [];
                   
                   % Plot results
                   figure();
                   scatter(time_points, monthly_P_avg,'.b');
                   hold on;
                   scatter(time_points, monthly_ET_avg,'.g');
                   legend('Precip','ET');
                   time_points = [0:0.01:1];

                   monthly_P_smooth = ppval(obj.climateData.functionParameters_P, time_points )';
                   monthly_ET_smooth = ppval(obj.climateData.functionParameters_ET, time_points )';                   
                   plot(time_points,monthly_P_smooth, '-b');
                   plot(time_points,monthly_ET_smooth, '-g');
                   xlabel('Fraction of year');
                   ylabel('Climate flux (mm/day)');
                   hold off;
                   return;
                   
               case {4,5}       % Monthly simulation using a smooth sigmoidal curve for percentiles of daily forcing rates.
                   % Calculate percentiles of daily forcing.
                   nyears = size(monthly_P,1);
                   nmonths = size(monthly_P,2);
                   monthly_P_pcnt = zeros(nmonths, nyears, npcntiles);
                   monthly_ET_pcnt = zeros(nmonths, nyears, npcntiles);
                   for irow=1:nyears
                       for imonth=1:nmonths
                             for ipcnt = 1: npcntiles
                                 if size( monthly_P{irow, imonth}, 1) > 0
                                     if sum(monthly_P{irow, imonth}) > 0                                        
                                        monthly_P_pcnt(imonth, irow, ipcnt ) =  prctile(monthly_P{irow, imonth}, 100.*scalingPercentiles(ipcnt) );
                                     end
                                     if sum(monthly_ET{irow, imonth}) > 0                                        
                                        monthly_ET_pcnt(imonth, irow, ipcnt ) =  prctile(monthly_ET{irow, imonth}, 100-100.*scalingPercentiles(ipcnt) );
                                     end
                                 else
                                     monthly_P_pcnt(imonth, irow, ipcnt ) = NaN;
                                     monthly_ET_pcnt(imonth, irow, ipcnt ) = NaN;
                                 end

                             end
                       end
                   end
 
                   monthly_P_pcnt = reshape(monthly_P_pcnt, nmonths, npcntiles );
                   monthly_ET_pcnt = reshape(monthly_ET_pcnt, nmonths, npcntiles );
                   monthly_P_pcnt = [monthly_P_pcnt(end,:); monthly_P_pcnt];
                   monthly_ET_pcnt = [monthly_ET_pcnt(end,:); monthly_ET_pcnt];
                   
                   % Calculate monthly mean forcing.
                   monthly_P_mean= zeros(12,1);
                   monthly_ET_mean= zeros(12,1);
                   monthly_Pcount= zeros(12,1);
                   monthly_ETcount= zeros(12,1);
                   for imonth=1:12
                   % Sum monthly forcing for latter calc of mean.
                       for irow=1:nyears
                           monthly_P_mean(imonth,1) = monthly_P_mean(imonth,1) + sum(monthly_P{irow, imonth});
                           monthly_Pcount(imonth,1) = monthly_Pcount(imonth,1) + size(monthly_P{irow, imonth},1);
                           
                           monthly_ET_mean(imonth,1) = monthly_ET_mean(imonth,1) + sum(monthly_ET{irow, imonth});
                           monthly_ETcount(imonth,1) = monthly_ETcount(imonth,1) + size(monthly_ET{irow, imonth},1);
                       end                       
                   end
                   
                   monthly_P_mean = monthly_P_mean./monthly_Pcount;
                   monthly_ET_mean = monthly_ET_mean./monthly_ETcount;
                                                        
                   % Fit curve to P and ET percentiles.
                   if userOption==4
                       obj.climateData.type = 'function';
                       obj.climateData.doClimateScaling = true;    
                       obj.climateData.scalingPercentiles = scalingPercentiles;
                       obj.climateData.data.P = [];
                       obj.climateData.data.ET = [];

                       time_points = [0; (cumsum(eomday(1999,1:12))./365)'];                       
                       time_points_ppfn = [0:0.01:1];
                       
                       obj.climateData.functionParameters_P = pchip([ time_points(end-1) - 1; time_points; 1+time_points(2)], [ monthly_P_pcnt(end-1,:); monthly_P_pcnt; monthly_P_pcnt(2,:) ]' );
                       obj.climateData.functionParameters_ET = pchip([ time_points(end-1) - 1; time_points; 1+time_points(2)], [ monthly_ET_pcnt(end-1,:); monthly_ET_pcnt; monthly_ET_pcnt(2,:) ]' );                                              
                       monthly_P_smooth = ppval(obj.climateData.functionParameters_P, time_points_ppfn )';
                       monthly_ET_smooth = ppval(obj.climateData.functionParameters_ET, time_points_ppfn )';
                       
                       figure();
                       %lsqcurvefit_options = optimset('Display','off');
                       for ipcnt = 1: npcntiles
                         
                           % Plot results                          
                           subplot(4,2,1:2);
                           scatter(time_points, monthly_P_pcnt(:,ipcnt),'.b');                       
                           hold on;
                           scatter(time_points, monthly_ET_pcnt(:,ipcnt),'.g');                                                                     
                          
                           plot(time_points_ppfn, monthly_P_smooth(:,ipcnt),'-b');
                           plot(time_points_ppfn, monthly_ET_smooth(:,ipcnt),'-g');                       
                           legend_str(1,4*(ipcnt-1)+1 : 4*(ipcnt-1)+4) = {['Precip. - ',num2str(scalingPercentiles(ipcnt)),' %ile'], ...
                                                                                            ['ET- ',num2str(scalingPercentiles(ipcnt)),' %ile'], ...
                                                                                            ['Fitted precip. - ',num2str(scalingPercentiles(ipcnt)),' %ile'], ...
                                                                                            ['Fitted ET- ',num2str(scalingPercentiles(ipcnt)),' %ile'] };

                       end
                       legend(legend_str);
                       legend off;
                       xlabel('Fraction of year');
                       ylabel('Climate flux (mm/day)');
                       box on;
                       hold off;
                   else
                       % Sigmoidal function is not fit the empirical cdf
                       obj.climateData.type = 'data';
                       obj.climateData.doClimateScaling = true;    
                       obj.climateData.scalingPercentiles = scalingPercentiles;

                       
                       % Derive monthly time points of climate data and convert to vectors.
                       nyears = max(inputClimateData(:,1)) - min(inputClimateData(:,1)) + 1;
                       nmonths = 12;
                       years = repmat(  (inputClimateData(1,1):inputClimateData(end,1))', 1, 12);
                       months = repmat(  1:nmonths, nyears , 1);
                       years = reshape(years', nyears* nmonths, 1);
                       months = reshape(months', nyears* nmonths, 1);
                       days = eomday( years, months);

                       % Aggregate year, month and day to a single number of
                       % year+fraction of year.
                       t_data = datenum(years,months,days) - datenum(years,1,1) + 1;
                       days_per_year = eomday(years, 2) + 365-28;
                       t_data = years + t_data./days_per_year;                      
                      
                       % Extract empirical percentile data for required
                       % months and assign to object.                       
                       obj.climateData.data.P = [t_data, monthly_P_pcnt(months,:)];
                       obj.climateData.data.ET = [t_data, monthly_ET_pcnt(months,:)];          
                       obj.climateData.functionParameters_P = [];
                       obj.climateData.functionParameters_ET = [];                        
                       
                       return
                   end
                   
                   % Add line plot CDFS of climate forcing.
                   month_labels = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
                   h = subplot(4,2,3);                    
                   plot(h, scalingPercentiles, monthly_P_pcnt(1:end-1,:)');
                   legend( month_labels{1:size(monthly_P_pcnt(1:end-1,1))}, 'location','Best');
                   legend off;
                   xlabel('Percentile');
                   ylabel('Precip. (mm/day)');
                   title('Monthly Precipitation from Empirical CDFs');
                   box on;
                   
                   h = subplot(4,2,4);                    
                   plot(h, scalingPercentiles, monthly_ET_pcnt(1:end-1,:)'); 
                   legend( month_labels{1:size(monthly_ET_pcnt(1:end-1,1))}, 'location','Best');
                   legend off;
                   xlabel('Percentile');
                   ylabel('ET. (mm/day)');
                   title('Monthly Precipitation from Empirical CDFs');
                   box on; 
                   
                   % To assess the adequacy of the number of percentile,
                   % a bar plot is created to summarise the difference in
                   % the mean monthly forcing and that from the integral of
                   % the empirical distributions and the fitted function.
                   monthly_Count= zeros(12,1);       
                   monthly_P_smooth = zeros(12, npcntiles);
                   monthly_ET_smooth = zeros(12, npcntiles);
                   monthly_P_smoothInt = zeros(12, 1);
                   monthly_ET_smoothInt = zeros(12, 1);
                   for imonth=1:12
                       % Integrate empirical monthly flux
                       monthly_P_empiricalInt(imonth,1) = trapz( scalingPercentiles, monthly_P_pcnt(imonth+1,:) );
                       monthly_ET_empiricalInt(imonth,1) = trapz(scalingPercentiles,  monthly_ET_pcnt(imonth+1,:) );
                       
                       % Integrate fitted functions for each percentile at the end of each
                       % month.
                       fractionOfYear = sum(eomday(0,1:imonth))/sum(eomday(0,1:12));
                       monthly_P_smooth(imonth, :) = ppval(obj.climateData.functionParameters_P, fractionOfYear )';
                       monthly_ET_smooth(imonth, :) = ppval(obj.climateData.functionParameters_ET, fractionOfYear )';
                       monthly_P_smoothInt(imonth,1) = trapz( scalingPercentiles, max(0,monthly_P_smooth(imonth,:) )' ); 
                       monthly_ET_smoothInt(imonth,1) = trapz(  scalingPercentiles, max(0,monthly_ET_smooth(imonth,:) )' );
                       
                   end

                   % Add line plot CDFS of fitted-function climate forcing.
                   month_labels = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
                   h = subplot(4,2,5);                    
                   plot(h, scalingPercentiles, monthly_P_smooth');
                   legend( month_labels, 'location','Best');
                   legend off;
                   xlabel('Percentile');
                   ylabel('Precip. (mm/day)');
                   title('Monthly Precipitation from Sigmoidal Functions fitted to CDFs');
                   box on;
                   
                   h = subplot(4,2,6);                    
                   plot(h, scalingPercentiles, monthly_ET_smooth'); 
                   legend( month_labels, 'location','Best');
                   legend off;
                   xlabel('Percentile');
                   ylabel('ET. (mm/day)');
                   title('Monthly ETo from Sigmoidal Functions fitted to CDFs');
                   box on;


                   h = subplot(4,2,7);                    
                   bar(h, [ monthly_P_mean, monthly_P_empiricalInt, monthly_P_smoothInt] ); 
                   legend('Observed','Empirical integration','Function integration');
                   set(h, 'XTickLabel', month_labels );
                   xlabel('Month');
                   ylabel('Mean monthly Precip. (mm/day)');
                   box on; 

                   h = subplot(4,2,8);
                   bar(h, [ monthly_ET_mean, monthly_ET_empiricalInt, monthly_ET_smoothInt] ); 
                   legend('Observed','Empirical integration','Function integration');
                   set(h, 'XTickLabel', month_labels );
                   xlabel('Month');
                   ylabel('Mean monthly Precip. (mm/day)');
                   box on; 
                   
                   return;
           end
           %-------------------------
           function result = prctile(data, percent)
               if percent==0
                   result = min(data);
               elseif percent==100
                   result = max(data);    
               else
                   N = length(data);
                   data = sort( data,'ascend');                   
                   data_prcntile = 100 .* (1:N)'./ N;
                   lowerRow = find(data_prcntile <= percent,1, 'last');
                   upperRow = find(data_prcntile > percent,1, 'first');               
                   result = data(lowerRow) + (percent - data_prcntile(lowerRow)) * ( data(upperRow) - data(lowerRow)  )/( data_prcntile(upperRow) - data_prcntile(lowerRow));
               end
           end
       end
       
%% Set state variable initial conditions
       function setInitialConditions(obj, input_initialConditions, ICType)           
% SETINITIALCONDITIONS Set model initial conditions
%
% Syntax:
%   setInitialConditions(obj, input_initialConditions, ICType)   
%
% Description:    
%   The method sets the state variable initial conditions. The state
%   variable set are the groundwater elevation at all model nodes and the
%   soil moisture (v/v) at all model nodes. Three methods are available for 
%   setting the initial conditions and each is defined by the input 'ICType'
%   and is decribed below: 
%
%   1. 'fractional': a 2x1 vector of the aquifer head as fraction of the 
%   maximum saturated thickness and volumetric soil moisture (v/v).
%
%   2. 'node_storage': a matrix mX2 for the groundwater storage [L^2] 
%   and soil moisture storage [L] where 'm' is the number of model nodes.
%   Importantly, the inputs are the model state variable values, not the
%   head (as an elevation) and fractional soil moisture. This option is
%   useful for input of the results from a previous simulation.
%
%   3. 'node_head': : a matrix mX2 for the groundwater elevation [L] 
%   and soil moisture storage (v/v) where 'm' is the number of model nodes.
%   
% Inputs:
%   obj - object of type 'CatchmentResilienceModel'.
%
%   input_initialConditions - input vector or matrix containing the initial
%   conditions.
%
%   ICType - string decribing the type of input condition. Options are:
%   'fractional', 'node_storage', 'node_storage'. See the description
%   section above for details of each.
%
% Outputs:
%   (model object, obj, is updated to contain the initial conditions.)
%
% See also:
%   CatchmentResilienceModel: class_description;
%
% Dependencies
%   (none)
%
% Author: 
%   Dr. Tim Peterson, The Department of Infrastructure
%   Engineering, The University of Melbourne.
%
% Date:
%   2 August 2012

            % The input 'ICType' is added to distinguish between 
            % the node aquifer values being of head, h, or storage, S.
            obj.initialConditions.type = ICType;
            switch ICType
                case 'fractional'  % vector 2x1 aquifer head as fraction of the maxium and volumetric soil mositure.
                    b_vert = input_initialConditions(1) * obj.geometryVariables.basementDepth ;
                    b = b_vert .* obj.geometryVariables.cos_BedAngle;
                    bmax = (obj.geometryVariables.basementDepth - obj.modelParameters.geometry.soilDepth) ...
                        .* obj.geometryVariables.cos_BedAngle;

                    obj.initialConditions.S = obj.modelParameters.geometry.width .* obj.geometryVariables.f ...
                        .* min(b,bmax) + (b > bmax).* (b - bmax) .* obj.modelParameters.geometry.width ...
                        .* obj.geometryVariables.theta_s;                    
                    obj.initialConditions.m = 1e3.*min(obj.modelParameters.geometry.soilDepth, obj.geometryVariables.basementDepth - b_vert).* input_initialConditions(2);                
                    
                case 'node_storage' % matrix nx X 2 of groundwater storage, S, and soil moisture storage, m.
                    obj.initialConditions.S = input_initialConditions(:,1);
                    obj.initialConditions.m = input_initialConditions(:,2);
                    
                case 'node_head'  % matrix nx X 2 of groundwater head elevation, h, and volumetric soil moisture storage, theta.
                    head = min(obj.modelParameters.geometry.surfaceElevation - 0.01 , input_initialConditions(:,1));
                    head = max( obj.modelParameters.geometry.basementElevation + 0.01, head);
                    b_vert = head - bed;

                    bmax_vert = D - obj.modelParameters.geometry.soilDepth;        
                    b = b_vert .* obj.geometryVariables.cos_BedAngle;
                    bmax = bmax_vert .* obj.geometryVariables.cos_BedAngle;

                    obj.initialConditions.S = obj.modelParameters.geometry.width .* obj.geometryVariables.f .* min(b,bmax) + (b > bmax).* (b - bmax) .* obj.modelParameters.geometry.width .* obj.geometryVariables.theta_s;
                    obj.initialConditions.m = 1e3.*min(d_vert, D - b_vert).* input_initialConditions(:,2);  
                otherwise
                    error('Input initial condition data type is not known!');
                    
            end
       end
       
%% Set differential equation solver options.
       function setSolverOptions(obj, optionsODE) 
% SETSOLVEROPTIONS Set differential equation solver options.
%
% Syntax:
%   setSolverOptions(obj)   
%   setSolverOptions(obj, optionsODE)   
%
% Description:    
%   The method sets the options for the differential equation solver. A
%   default structure can be set, or modifications can be made to this and
%   applied to the object.
%
%   Importantly, the results within Peterson, Western and Argent
%   (2012) used an ODE solver that completed simulation at the
%   end of each climate forcing time step. This was achieved
%   using a modified version of ODE15s. However, this code
%   cannot be distributed because it would likely be in
%   violation of The Deployment Addendum of The MathWorks,
%   Inc. Software License Agreement. Hence, in the distributed
%   version of this code, the standard version of ODe15s is
%   used. To ensure that each simulation time-step finishes at
%   the end of each climate time step, the solver is re-called
%   for each clime time step. 
%   
%   To modify the standard solver 'ode15s' to that used within Peterson, 
%   Western and Argent (2012), the following modifications can be made to
%   ode15s.m. Once modified, the file MUST be renamed to 'ode15s_StepLimit'
%   for it to be used AND must be placed into the same directory as
%   'ode15s.m'.
%
%   Replace the following code: 
%       % THE MAIN LOOP
% 
%       done = false;
%       at_hmin = false;
%       while ~done
%   
%           hmin = 16*eps(t);
%           absh = min(hmax, max(hmin, absh));
%           if absh == hmin
%               if at_hmin
%                   absh = abshlast;  % required by stepsize recovery
%               end  
%               at_hmin = true;
%           else
%               at_hmin = false;
%           end  
%           h = tdir * absh;
%
%   With the following modification:
%          % Variable to signal that a solution at end of timestep.
%          hIsCapped=false;
%
%          % THE MAIN LOOP 
%          done = false;
%          at_hmin = false;
%          while ~done
%               [irow,icol] = find( (t<tspan),1,'first');
%               t_endOfStep = tspan(irow,1);
%               hmin = 16*eps*abs(t);
%               if hIsCapped && abs(t - t_endOfStep) < hmin;
%                   absh = min(hmax, max(abshprev,absh) );        
%                   hIsCapped=false;
%               elseif  hIsCapped && t ==  tspan(irow-1,1);         %previous step was end of discrete period (ie month/yr)
%                   absh = min(hmax, min(t_endOfStep -tdir*t, max(abshprev,absh) ) );      
%                   hIsCapped=false;
%               elseif t_endOfStep - t < tdir*absh*1.1 && t_endOfStep - t <hmax; 
%                   abshprev = absh;
%                   h= t_endOfStep -tdir*t;  
% 
%                   % get exact time point 
%                   temp = h + t_endOfStep;
%                   h = temp - t_endOfStep;
% 
%                   absh=abs(h);      
%                   hIsCapped=true;
%               else
%                   abshprev = min(abs(absh),hmax);
%                   hIsCapped=false;
%               end  
% 
%               absh = min(hmax, max(hmin, absh));
%               h = tdir * absh;
%
% Inputs:
%   obj - object of type 'CatchmentResilienceModel'.
%
%   optionsODE - structure variable containing a full set of solver options
%   modified to the user's requirements. To obtain an example of the
%   default structure variable, build a model object and call this method
%   without the input 'optionsODE'. The default structure is given within
%   'obj.solverOptions' where 'obj' is the name of the model.
%
% Outputs:
%   (model object, obj, is updated to contain the solver options.)
%
% See also:
%   CatchmentResilienceModel: class_description;
%
% Dependencies
%   (none)
%
% Author: 
%   Dr. Tim Peterson, The Department of Infrastructure
%   Engineering, The University of Melbourne.
%
% Date:
%   2 August 2012

           % Clear old options
           obj.solverOptions = odeset();
           
           if nargin==1
               display(' The default solver options have been adopted.');
               % Set standard ODE options
               obj.solverOptions = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);                              
               
               % Set solver type
               obj.solverOptions.solverType = 'ode15s';
               
               % Add a new option that specifies the maximum step size as a
               % fraction of the climate forcing step size.
               obj.solverOptions.MaxStep_Scaler = 0.49;                              
               
               % Define if the model should be treated an a differential
               % algebraic equation.
               obj.solverOptions.useDAE = false;
               
               % Show warnings issued by the model
               obj.solverOptions.showWarnings = true;
               
               % Note, the results within Peterson, Western and Argent
               % (2012) used an ODE solver that completed simulation at the
               % end of each climate forcing time step. This was achieved
               % using a modified version of ODE15s. However, this code
               % cannot be distributed because it would likely be in
               % violation of The Deployment Addendum of The MathWorks,
               % Inc. Software License Agreement. Hence, in the distributed
               % version of this code, the standard version of ODe15s is
               % used. To ensure that each simulation time-step finishes at
               % the end of each climate time step, the solver is re-called
               % for each clime time step. 
               %optionsODE.solverType = 'ode15s_StepLimit';
           else
               
               if ~isfield(optionsODE,'solverType')
                   display(' Input solver options must contain a filed called "solverType"');
                   display(' The default solver "ode15s_StepLimit" has been adoped.');
                   optionsODE.solverType = 'ode15s_StepLimit';
               else
                   % Add the user set solver options.
                   obj.solverOptions = optionsODE;                      
               end
           end
           
       end
       
%% Call ODE solver
       function solveModel(obj, ta, tb, nForcingSubSteps)
% SOLVEMODEL Time-integration solution of the model.
%
% Syntax:
%   solveModel(obj, ta, tb, nForcingSubSteps)
%
% Description:              
%   Solves the non-linear differential equation model from time point 'ta'
%   to time point 'tb' using the user input climate data, initial
%   conditions and differential equation solveroptions. Once the model is
%   successfully solved, the following can be undertaken:
%       - plot cross sections of head and soil moisture;
%       - plot time series of a large number of fluxes, variables, state 
%       variables, inputs; and
%       - undertake continuation analysis for estimation of attractors and
%       repellors.
%
%   Importantly, the model is a complex non-linear partial differential
%   equation that can become unstable. This is most likely when the
%   groundwater storage approaches the land surface or the aquifer
%   basement. The stability can be improved by increasing the spatial
%   discretization, but some parameter combinations may still produce
%   unstable solutions. At the aquifer basement, instability can occur if
%   the lateral flow at the boundary is too high, causing the model nodes
%   to become dry or have a negative groundwater storage. To address this
%   problem, the following model parameters can be modified:
%
%       - 'obj.modelParameters.aquifer.Units{i,1}.d_khalf' controls the
%       saturated depth (at geological unit i) at which k_sat is half the 
%       maximum. Increasing this parameter will increase model stability;
%
%       - 'obj.modelParameters.aquifer.Units{i,1}.tor' controls
%       the rate at which the k_sat declines as the groundwater approaches
%       the basement. Reducing this parameter will increase model
%       stability.
%
%   Also, as the water table approaches the land surface the model becomes
%   unstable because of a singularity. Specifically, as the water table
%   approaches the surface, the thickness of the unsaturated zone declines
%   accordingly. With a reduced soil moisture capacity, a smaller rainfall
%   event is required to cause saturation of the vadose zone. Considering
%   that the head by which the water table rises as a result of recharge is
%   a function of the factional soil moisture, a near saturated soil moisture
%   could cause an extremely large rate of change in the groundwater
%   level. Overcoming this instability is challenging. Options to increase
%   the stability include modifying the model configuration (e.g. boundary
%   conditions) or reconfiguring the equations to a differential algebraic
%   equation (DAE). The latter can be achieved by changing the solver
%   option 'obj.solverOptions.useDAE' to 'true'.
%
% Inputs:
%   obj - object of type 'CatchmentResilienceModel'.
%
%   ta - scalar integer of the start year for simulation.
%
%   tb - scalar integer of the end year for simulation.
%
%   nForcingSubSteps - scalar integer >0 defining the number of sub-steps
%   from that of the forcing data to report. For example, if the forcing data
%   is observed monthly data and 'nForcingSubSteps=4' then solutions will
%   be reported approximately every week. If the forcing rate is defined by a
%   smooth function (e.g. a spline or sine function) then 'nForcingSubSteps'
%   is the number of simulations per year. This feature was included to
%   facilitate limit-cycle continuation.
%
% Outputs:
%   (model object, obj, is updated to contain the time-integration solution.)
%
% See also:
%   CatchmentResilienceModel: class_description;
%
% Dependencies
%   (none)
%
% Author: 
%   Dr. Tim Peterson, The Department of Infrastructure
%   Engineering, The University of Melbourne.
%
% Date:
%   2 August 2012

           % Update model constants
           setGeometryVariables(obj);                                  
           
           % Convert initial conditions to a vector    
           y0(1 : 2 : 2 * obj.geometryVariables.nx  ,1) = obj.initialConditions.S;
           y0(2 : 2 : 2 * obj.geometryVariables.nx+1 , 1) = obj.initialConditions.m;
           
           % Check solution time range against forcing data.
           if ta>=tb;
              error('When extracting the climate date for PDE solver, tb  must be greater than ta');
           elseif strcmp(obj.climateData.type ,'data')
               if ta <= obj.climateData.data.P(1,1);
                  error(' The start date, ta, must be greater than the first data within the climate data series');
               elseif tb> obj.climateData.data.P(end,1);
                  error(' The end date, tb, is greater than the end data within the climate data series');
               end 
           end 
                     
           if strcmp(obj.climateData.type ,'data')
               % Get solution time points from forcing data.
               t = obj.climateData.data.P(:,1);
               t_filt  = t >= ta & t <= tb;
               t = t(t_filt);

               % Sub-divide time steps from forcing data in nForcingSubSteps
               % sub-steps.
               t = [t(1:end-1), t(2:end)];
               for i=1: size(t,1); 
                   t(i,1:max(1,nForcingSubSteps)+1) = t(i,1) : (t(i,2) - t(i,1))/max(1,nForcingSubSteps) : t(i,2); 
               end
               t = t(:,1:end-1);
               t = reshape(t, prod(size(t)), 1);
               t = sort(t, 'ascend');               
           else
               if nForcingSubSteps<=1;
                   error('The number of forcing timesteps per year, when using a smooth function and not data, must be greater than 1');
               end
               t = ta: 1/nForcingSubSteps: tb;
           end
           hmax = min(diff(t)).*min(1,obj.solverOptions.MaxStep_Scaler);
           
           % Set Jacobian sparisity pattern
           ny = obj.geometryVariables.nx * 2;
           dFdy = ones(ny,6);
           dFdy(2:2:ny,1) = 0; 
           dFdy(2:2:ny,2) = 0; 
           dFdy(2:2:ny,3) = 0;           
           dFdy(2:2:ny,6) = 0;    
           dFdy = spdiags(dFdy,-3:2,ny,ny);   
           
           % Make temporary copy of solver options
           optionsODE = obj.solverOptions;   
           
           % Clear prior solver options and set new solver options           
           if obj.solverOptions.useDAE
               % Set mass matrix sparisity pattern
               diags = repmat([1;0], obj.geometryVariables.nx,1);
               diags(:,2) = repmat([0;1], obj.geometryVariables.nx,1);
               massPattern = spdiags( diags ,[0, 1],ny,ny);
               
               % Set solver options
               obj.solverOptions = odeset(obj.solverOptions,'Vectorized','on', 'MaxStep',hmax, 'JPattern', dFdy, 'Mass',@getDerivatives_DAEmass,'MassSingular','maybe','MvPattern',massPattern,'MStateDependence','weak');                                     
               
               % Output messaage
               display('... Solving as differential-algebraic equation, ie using mass matrix');
           else
               % Set solver options
               obj.solverOptions = odeset(obj.solverOptions,'Vectorized','on', 'MaxStep',hmax, 'JPattern', dFdy);      
               
               % Output messaage
               display('... Solving as non-differential-algebraic equation, ie not using mass matrix');
           end
           
           % Create options structure and handle the Jacobian.
           currentDir = pwd;
           if isunix || ismac
                cd([matlabroot , '/toolbox/matlab/funfun/private/']);          
           elseif ispc
               cd([matlabroot , '\toolbox\matlab\funfun\private\']);          
           else
               error('The code only recognises Unix. Mac and Windows operating systems!');
           end
           [junk_Jconstant,junk_Jac,junk_Jargs, Joptions] = feval('odejacobian', true,@getDerivatives_DAE, ta, y0, obj.solverOptions,obj);
           cd(currentDir);
           
           % Add original options back into obj.solverOptions
           obj.solverOptions.MaxStep_Scaler = optionsODE.MaxStep_Scaler;
           obj.solverOptions.solverType = optionsODE.solverType;
           obj.solverOptions.useDAE = optionsODE.useDAE;
           obj.solverOptions.showWarnings = optionsODE.showWarnings;
           clear optionsODE;
           
           % Add Jacobian options to obj.solverOptions
           obj.solverOptions.Jacobian_Options = Joptions;
           
           % Solve PDE   . 
           % ------------------------------------------             
           warning on;                            
           if strcmp(obj.solverOptions.solverType, 'ode15s_StepLimit')
               [t, y] = feval(obj.solverOptions.solverType, @getDerivatives_DAE, t, y0, obj.solverOptions, obj, Joptions.g, dFdy );
           else
               % This less efficient solver is used if the modified solver
               % 'ode15s_StepLimit' cannot be legally supplied to the user.
               [tout, y] = feval(obj.solverOptions.solverType, @getDerivatives_DAE, t(1:3), y0, obj.solverOptions, obj, Joptions.g, dFdy );
               for i=4:length(t)
                    [tout_temp, y_temp] = feval(obj.solverOptions.solverType, @getDerivatives_DAE, t(i-1:i), y(i-1,:), obj.solverOptions, obj, Joptions.g, dFdy );
                    tout(i) = tout_temp(end);
                    y(i,:) = y_temp(end,:);
               end
               t=tout;
           end
           % ------------------------------------------ 
           
           % Extract state variables from vector, y.
           npde=2;
           y=y';
           S = y( 1:npde:npde * obj.geometryVariables.nx-1 , :);
           m= y( 2:npde:npde * obj.geometryVariables.nx , :);          
           
           % Add results to object, obj.
           % ------------------------------------------
           obj.results.S = S; 
           obj.results.m = m; 
           obj.results.t = t';
           obj.results.u = obj.modelParameters.geometry.u_nodes;
           obj.results.x = obj.geometryVariables.x;           
           
       end
       
%% Plot cross sections of the model solutions
       function plotCrossSections(obj, plotTimePoints, plotHeads, plotSoilMositure)
% PLOTCROSSSECTIONS Plot cross sections of the model solutions.
%
% Syntax:
%   plotCrossSections(obj, plotTimePoints, plotHeads, plotSoilMositure)
%
% Description:
%   Plots cross section of the model solutions at user defined time points.
%   The user can specify to plot the groundwater head, soil moisture or
%   both. At least three solution time points are required for plotting.
%
% Inputs:
%   obj - object of type 'CatchmentResilienceModel'.
%
%   plotTimePoints - vector of time points to plot. 
%
%   plotHeads - logical scalar for specifying if the groundwater level
%   should be plotted where 'true' to plot and 'false' to not plot.
%
%   plotSoilMositure - logical scalar for specifying if the soil moisture 
%   should be plotted where 'true' to plot and 'false' to not plot.
%
% Outputs:
%   (none)
%
% See also:
%   CatchmentResilienceModel: class_description;
%
% Dependencies
%   (none)
%
% Author: 
%   Dr. Tim Peterson, The Department of Infrastructure
%   Engineering, The University of Melbourne.
%
% Date:
%   2 August 2012

           % re-calculate model geometry data just in case.
           setGeometryVariables(obj);
           
           % Add initial head time point if not already included
           if obj.results.t(1) ~= plotTimePoints(1);
               plotTimePoints = [obj.results.t(1), plotTimePoints];
           end           
           
           j=0;
           for i=1:max(size(plotTimePoints))
              
               % Filter solution time points to those of t_steps.
               [trow, tcol] = find(obj.results.t >= plotTimePoints(i),1,'first');
               filter = false(size(obj.results.t));
               filter(1,tcol)=true;               
               
               % If there exists a time point then assign state variables
               % to the object and calculate the model variables 
               if max(filter)>0
                   j=j+1;
                   obj.stateVariables.S = obj.results.S(:,filter);
                   obj.stateVariables.m = obj.results.m(:,filter);
                   obj.stateVariables.t = obj.results.t(filter);
                   
                   % calculate state variable dependent variables, namely
                   % heads and soil moisture.
                   setVariables(obj);
                   
                   % Extract data for cross sections and add to matrix
                   timePoints(j) = obj.stateVariables.t;
                   heads(:,j) = obj.variables.h_vert;
                   theta(:,j) = obj.variables.theta;
               end
           end           
           %-------------------------
           
           % Check at least three time points were extracted for plotting
           if size(heads,2)<3
               error(['An insufficient number of time points (<3) where available between the input start and end years. The first and last years within the solution are ' , num2str(obj.results.t(1)), ' and ' , num2str(obj.results.t(end))]);
           end
           
           % Plot cross section
           if plotHeads;
              figure();
              
              % Assemble all cross-section elevation data and heads
              elevationData = [obj.modelParameters.geometry.surfaceElevation, ...                       
                       obj.modelParameters.geometry.surfaceElevation - obj.modelParameters.geometry.soilDepth];
              h1 = plot( obj.results.u, elevationData, 'Color', [.6 .6 .6], 'LineWidth', 0.5, 'LineStyle', '--' );
              hold on;
              h2 = area( obj.results.u, obj.modelParameters.geometry.basementElevation);
              h3 = plot( obj.results.u, [heads(:,1:2), heads(:,end), heads(:,3:end-1)], 'Color', [.6 .6 .6], 'LineWidth', 0.5, 'LineStyle', '--' );
              
              
              % format each series in addition to the transient sets
              set(h1(1,1), 'Color','k','LineStyle','-','LineWidth',0.75);  % natural surface
              set(h1(2,1), 'Color','k','LineStyle','-.','LineWidth',0.5);   % soil - aquifer layer
              set(h2, 'FaceColor',[0.8, 0.8, 0.8],'EdgeColor','k');   % bedrock
              set(h3(1,1), 'Color','k','LineStyle','--','LineWidth',0.75);   % initial head
              set(h3(3,1), 'Color','k','LineStyle','-','LineWidth',2);   % final head    
              
              % Legend labels
              legend({'Land surface','Soil-aquifer boundary','Aquifer basement', ...
                  'Initial head','Head at time steps','Final head'}, 'Location','NorthWest'); 
              
              % some plot labels
              xlabel('Distance from outlet (m)');
              ylabel('Elevation (m)');
              box('on');              
              
              % set order of objects
              uistack(h3(1,1),'top');
              uistack(h3(3,1),'top');

              % set axis range to tight
              axis tight;              
              hold off;                                       
           end                    
           
           if plotSoilMositure;
              figure();
              h1 = plot( obj.results.u, [theta(:,1:2), theta(:,end), theta(:,3:end-1)], 'Color', [.6 .6 .6], 'LineWidth', 0.5, 'LineStyle', '--' );
              hold on;
              
              set(h1(1,1), 'Color','k','LineStyle','--','LineWidth',0.75);   % initial theta
              set(h1(3,1), 'Color','k','LineStyle','-','LineWidth',2);   % final theta    
              
              % Legend labels
              legend({'Initial \theta_{v/v}','\theta_{v/v} at time steps','Final \theta_{v/v}'}, 'Location','NorthWest'); 
              
              % some plot labels
              xlabel('Distance from outlet (m)');
              ylabel('Volumetric soil moisture (\theta_{v/v})');
              box('on');              
              
              % set order of objects
              uistack(h1(1,1),'top');
              uistack(h1(3,1),'top');

              % set axis range to tight
              axis tight;              
              hold off;  
              
           end
       end
       
%% Plot time-series of the model solutions
       function fluxesAndVariables = plotTimeSeries(obj, fluxesAndVariablesToPlot, t_start, t_end, nodesToPlot)
% PLOTTIMESERIES Plot time-series of model solutions, variables and fluxes.
%
% Syntax:
%   fluxesAndVariables = plotTimeSeries(obj)
%   fluxesAndVariables = plotTimeSeries(obj, fluxesAndVariablesToPlot, t_start, t_end, nodesToPlot)
%
% Description:
%   Plots time-series of the model solutions, variables or fluxes at user
%   specified model spatial nodes. A list of possible time-series plots can
%   be obtained by calling the method with only the model object input
%   (e.g. 'fluxesAndVariables = plotTimeSeries(obj)'). This will return an
%   array of strings where each string is an available type of plot. For
%   an example of how to use this method, see the example with the
%   documentation for the class description below. Alternatively, enter the
%   following into the MatLab command window: 'doc CatchmentResilienceModel'.
%
%   The method can also plot an approximation of the water balance errors.
%   This is undertaken by temporal interpolation over each time-step. To
%   ensure a reliable interpolation is obtained (and thus a reliable
%   estimate of the water balance error) the model should be solved with a
%   large number of sub-time steps compared to that of the forcing
%   time-steps (see the documentation for 'solveModel' and its input 
%   'nForcingSubSteps').
%
% Inputs:
%   obj - object of type 'CatchmentResilienceModel'.
%
%   fluxesAndVariablesToPlot - cell vector of strings containing the
%   variables, fluxes rtc to plot.
%
%   t_start - scalar double of the start date for simulation.
%
%   t_end - scalar double of the end date for simulation.
%
%   nodesToPlot - vector of the model spatial nodes at which plots are 
%   required. This is input as a distance, not model node number. 
%
% Outputs:
%   fluxesAndVariables - variable type output depending upon the inputs to
%   the method. If only the model object is input, then a cell array of
%   strings is returned. Each string is the name of a variable, flux etc
%   that can be plotted. If all five inputs are specified, then
%   'fluxesAndVariables' is a structure variable containing the time-series
%   data that was plotted.
%
% See also:
%   CatchmentResilienceModel: class_description;
%
% Dependencies
%   (none)
%
% Author: 
%   Dr. Tim Peterson, The Department of Infrastructure
%   Engineering, The University of Melbourne.
%
% Date:
%   2 August 2012

           % If only the model object is inputs, then return a list of all fluxes able to be
           % plotted.
           if nargin==1;
               if isfield(obj.results,'S')
                   fluxes = fieldnames(obj.fluxes); 
                   variables = fieldnames(obj.variables); 
                   fluxesAndVariables = [fluxes; variables; 'massBalance', 'S', 'm'];
                   return;
               else
                   error(' A model solution must first be derived in order to plot he fluxes!');
               end               
           end                     
           
           % Extract time points for plotting
           timePoints_filter = obj.results.t >=t_start & obj.results.t <=t_end;
           timePoints = obj.results.t(timePoints_filter);
           ntimePoints = size(timePoints,2);
           
           % Re-calculate model geometry data just in case.
           setGeometryVariables(obj);          
           
           % Create output structure of fluxes and check that the
           % requested fluxes to be plotted are an output of the model.
           doMassBalance = false;
           for i=1: length(fluxesAndVariablesToPlot)
              if isfield(obj.fluxes, fluxesAndVariablesToPlot(i)) ... 
              || isfield(obj.variables, fluxesAndVariablesToPlot(i)) ...
              || strcmp('S', fluxesAndVariablesToPlot(i))  ...
              || strcmp('m', fluxesAndVariablesToPlot(i))  ...   
              || strcmp('massBalance', fluxesAndVariablesToPlot(i))                
                  fluxesAndVariables.(char(fluxesAndVariablesToPlot{i})) = [];
              else
                  warning(['... The following input flux/variable label is not an output of the model. Please check its spelling. Errorous input name: ', char(fluxesAndVariablesToPlot(i))]);                      
              end
              
              % Also check if mass balance calculations required.
              if strcmp(fluxesAndVariablesToPlot(i), 'massBalance')
                  doMassBalance = true;
              end
              
           end        
           
           % Calculate fluxes for each time point during one call
           % to 'getDerivatives_DAE'.
           %-----------------------------
           
           % Assign state variables
           obj.stateVariables.S = obj.results.S(:,timePoints_filter);
           obj.stateVariables.m = obj.results.m(:,timePoints_filter);
           obj.stateVariables.t = obj.results.t(timePoints_filter);            
           
           % Calculate ALL fluxes by calling getDerivatives_DAE.
           y(1 : 2 : 2 * obj.geometryVariables.nx  ,:) = obj.stateVariables.S;
           y(2 : 2 : 2 * obj.geometryVariables.nx+1 , :) = obj.stateVariables.m;
           getDerivatives_DAE(obj.stateVariables.t, y, obj);                    
           
           % Aggregate fluxes and variables to temp variable
           allFluxesVariables = obj.fluxes;
           temp_variableNames = fieldnames(obj.variables);
           for j=1:length(temp_variableNames);
               allFluxesVariables.(char(temp_variableNames{j})) = obj.variables.(char(temp_variableNames{j}));
           end
           allFluxesVariables.S = obj.stateVariables.S;
           allFluxesVariables.m = obj.stateVariables.m;

           % Convert unit of some variables from years to days
           allFluxesVariables.k = allFluxesVariables.k./365.25;
           allFluxesVariables.boundaryFlow_lower = allFluxesVariables.boundaryFlow_lower./365.25;
           allFluxesVariables.boundaryFlow_upper = allFluxesVariables.boundaryFlow_upper./365.25;
           allFluxesVariables.Q_toStream = allFluxesVariables.Q_toStream./365.25;               

           % Extract fluxes required for plotting 
           for j=1: length(fluxesAndVariablesToPlot)
              if isfield(allFluxesVariables, fluxesAndVariablesToPlot(j))
                  fluxesAndVariables.(char(fluxesAndVariablesToPlot(j))) = allFluxesVariables.(char(fluxesAndVariablesToPlot{j}));   
              end
           end                    
           %-----------------------------
           
           % Calculate mass balance
           %-----------------------------
           if doMassBalance;
               bsxops(1);
               
               % Time step (in units of years)
               dt = repmat(timePoints(2:end)-timePoints(1:end-1), obj.geometryVariables.nx ,1);
               
               % Linearly integrate fluxes over time step.
               if strcmp( obj.climateData.type, 'data')
                   massBalance.P =  365.25 .* dt .* allFluxesVariables.P(:,1:end-1);
               else
                   massBalance.P = 0.5.* 365.25 .* dt .* (allFluxesVariables.P(:,1:end-1) + allFluxesVariables.P(:,2:end) );                   
               end
               massBalance.intercept = -0.5.* 365.25 .* dt .* (allFluxesVariables.intercept(:,1:end-1) + allFluxesVariables.intercept(:,2:end) );
               massBalance.satRunoff = -0.5.* 365.25 .* dt .* (allFluxesVariables.satRunoff(:,1:end-1) + allFluxesVariables.satRunoff(:,2:end) );               
               massBalance.evap = -0.5.* 365.25 .* dt .* (allFluxesVariables.evap(:,1:end-1) + allFluxesVariables.evap(:,2:end) );
               massBalance.transp = -0.5.* 365.25 .* dt .* (allFluxesVariables.transp(:,1:end-1) + allFluxesVariables.transp(:,2:end) );
               massBalance.evap_gw = -0.5.* 365.25 .* dt .* (allFluxesVariables.evap_gw(:,1:end-1) + allFluxesVariables.evap_gw(:,2:end) );
               massBalance.boundaryFlow_lower = -0.5.* 365.25 .* dt(1,:) .* (allFluxesVariables.boundaryFlow_lower(:,1:end-1) + allFluxesVariables.boundaryFlow_lower(:,2:end) );
               massBalance.boundaryFlow_upper = 0.5.* 365.25 .* dt(1,:) .* (allFluxesVariables.boundaryFlow_upper(:,1:end-1) + allFluxesVariables.boundaryFlow_upper(:,2:end) );
               massBalance.Q_toStream = -0.5.* 365.25 .* dt .* (allFluxesVariables.Q_toStream(:,1:end-1) + allFluxesVariables.Q_toStream(:,2:end) );               
               
               % Calculate fluxes for each time point where the climate
               % changes.
               if strcmp( obj.climateData.type, 'data')                    
                    
                   % Get climate data time points
                   climateTime  = obj.climateData.data.P(:,1);
                   climateTime_filter  = climateTime >=t_start & climateTime <=t_end;
                   climateTime = climateTime(climateTime_filter)';
                   
                   % Create filter within 'timePoints' for those at which
                   % the climate step ends.
                   timePoints_climateStepPoints = false(size(timePoints_filter));                   
                   flux_climateStepPoints = false(size(timePoints));                   
                   for i=1:size(climateTime,2)                       
                       [irow, icol] = find( timePoints == climateTime(i), 1, 'first');
                       if icol < size(timePoints,2); 
                           flux_climateStepPoints(irow,icol)=true; 
                           
                           [irow, icol] = find( obj.results.t == climateTime(i), 1, 'first');
                           timePoints_climateStepPoints(irow,icol)=true;
                       end
                   end
                   flux_climateStepPoints_ind = find(flux_climateStepPoints);

                   % Calculate fluxes at the start of each time step using
                   % the climate data from the forthcoming step.
                   min_dt = min(diff(obj.results.t(:,timePoints_climateStepPoints)));
                   climate = getClimate(obj, obj.results.t(:,timePoints_climateStepPoints ) + min_dt/2);
                   setFluxes(obj, obj.results.S(:,timePoints_climateStepPoints ), obj.results.m(:,timePoints_climateStepPoints ), obj.results.t(:,timePoints_climateStepPoints), climate);   
                   allFluxesVariables_ClimateStep = obj.fluxes;

                   % Convert unit of some variables from years to days
                   allFluxesVariables_ClimateStep.boundaryFlow_lower = allFluxesVariables_ClimateStep.boundaryFlow_lower./365.25;
                   allFluxesVariables_ClimateStep.boundaryFlow_upper = allFluxesVariables_ClimateStep.boundaryFlow_upper./365.25;
                   allFluxesVariables_ClimateStep.Q_toStream = allFluxesVariables_ClimateStep.Q_toStream./365.25;                                                                       
                   
                   % Re-derived flux integreations at time points just
                   % after a change in climate forcing.
                   dt_ClimateStep = dt(:,flux_climateStepPoints_ind);
                   massBalance.P(:, flux_climateStepPoints_ind) =  365.25 .* dt_ClimateStep .* allFluxesVariables_ClimateStep.P;
                   massBalance.intercept(:, flux_climateStepPoints_ind) = -0.5.* 365.25 .* dt_ClimateStep .* (allFluxesVariables.intercept(:,flux_climateStepPoints_ind+1) + allFluxesVariables_ClimateStep.intercept );
                   massBalance.satRunoff(:, flux_climateStepPoints_ind) = -0.5.* 365.25 .* dt_ClimateStep .* (allFluxesVariables.satRunoff(:,flux_climateStepPoints_ind+1) + allFluxesVariables_ClimateStep.satRunoff );
                   massBalance.evap(:, flux_climateStepPoints_ind) = -0.5.* 365.25 .* dt_ClimateStep .* (allFluxesVariables.evap(:,flux_climateStepPoints_ind+1) + allFluxesVariables_ClimateStep.evap );
                   massBalance.transp(:, flux_climateStepPoints_ind) = -0.5.* 365.25 .* dt_ClimateStep .* (allFluxesVariables.transp(:,flux_climateStepPoints_ind+1) + allFluxesVariables_ClimateStep.transp );
                   massBalance.evap_gw(:, flux_climateStepPoints_ind) = -0.5.* 365.25 .* dt_ClimateStep .* (allFluxesVariables.evap_gw(:,flux_climateStepPoints_ind+1) + allFluxesVariables_ClimateStep.evap_gw );
                   massBalance.boundaryFlow_lower(:, flux_climateStepPoints_ind) = -0.5.* 365.25 .* dt_ClimateStep(1,:) .* (allFluxesVariables.boundaryFlow_lower(:,flux_climateStepPoints_ind+1) + allFluxesVariables_ClimateStep.boundaryFlow_lower );
                   massBalance.boundaryFlow_upper(:, flux_climateStepPoints_ind) = 0.5.* 365.25 .* dt_ClimateStep(1,:) .* (allFluxesVariables.boundaryFlow_upper(:,flux_climateStepPoints_ind+1) + allFluxesVariables_ClimateStep.boundaryFlow_upper );                                      
                   massBalance.Q_toStream(:, flux_climateStepPoints_ind) = -0.5.* 365.25 .* dt_ClimateStep(1,:) .* (allFluxesVariables.Q_toStream(:,flux_climateStepPoints_ind+1) + allFluxesVariables_ClimateStep.Q_toStream);
                end

                % convert to horizontal plane
                S = obj.results.S(:,timePoints_filter)./obj.geometryVariables.cos_BedAngle;
                m =  obj.results.m(:,timePoints_filter); 

                % Calculate node value change in soil moisture storages, m, and
                % groundwater storage, S.                
                dm = 1e-3.*( m(:,2:end)  - m(:,1:end-1) );
                dS = S(:,2:end) - S(:,1:end-1);

                % Calculate net flux at land surface (converted to metres) for each node.
                massBalance.Surface_out_catchment = 1e-3 .* (massBalance.P + massBalance.intercept + massBalance.satRunoff + massBalance.evap + massBalance.transp + massBalance.evap_gw ); 

                % Integrate storage changes and land surface flux over the
                % catchment length.
                dm = simp_extend( obj.modelParameters.geometry.u_nodes , obj.modelParameters.geometry.width .* dm );
                dS = simp_extend( obj.modelParameters.geometry.u_nodes , dS);    
                surfaceFlux = simp_extend( obj.modelParameters.geometry.u_nodes , obj.modelParameters.geometry.width .* massBalance.Surface_out_catchment );
                baseFlow = simp_extend( obj.modelParameters.geometry.u_nodes , massBalance.Q_toStream );

                % Calculate net boundary condition flux
                boundaryFlux = massBalance.boundaryFlow_lower+ massBalance.boundaryFlow_upper;

                % Calculate catchment wide mass balance and convert to mm
                % per time step size.
                area = simp_extend( obj.modelParameters.geometry.u_nodes , obj.modelParameters.geometry.width);
                fluxesAndVariables.massBalance = 1e3.*(surfaceFlux + boundaryFlux + baseFlow - dm - dS)./area;
                clear area surfaceFlux boundaryFlux allFluxesVariables_ClimateStep dm dS;

                % Convert mass balance to mm/day.
                fluxesAndVariables.massBalance = fluxesAndVariables.massBalance ./( 365.25 .* dt(1,:)); 
                
                bsxops(0);
           end
           clear allFluxesVariables
           
           % Extract fluxes at required nodes. This first requires creation
           % of a filter for the required nodes.
           nodes_filter = false(size(obj.results.u));
           for i=1: length(nodesToPlot)
               [irow, icol] = find(obj.results.u == nodesToPlot(i), 1, 'first');
               nodes_filter(irow, icol) = true;
           end
           nNodes = sum(nodes_filter);
           
           nfluxes = 0;
           for i=1: length(fluxesAndVariablesToPlot)
                if isfield(fluxesAndVariables, fluxesAndVariablesToPlot(i))
                    % Check if field is of only one node.
                    if isfield(fluxesAndVariables, 'P') ...
                    || isfield(fluxesAndVariables, 'boundaryFlow_lower') ...
                    || isfield(fluxesAndVariables, 'boundaryFlow_upper') ...
                    || isfield(fluxesAndVariables, 'massBalance') 
                        fluxesAtNodes.(char(fluxesAndVariablesToPlot{i})) = fluxesAndVariables.(char(fluxesAndVariablesToPlot{i}));
                    else
                        fluxesAtNodes.(char(fluxesAndVariablesToPlot{i})) = fluxesAndVariables.(char(fluxesAndVariablesToPlot{i}))(nodes_filter,:);
                    end
                    nfluxes = nfluxes+1;
                end
           end

           % Plot fluxes over require time range and required nodes
           figure();           
           plotNumber = 0;
           for i=1: nfluxes;
               for j=1: nNodes
                   plotNumber = plotNumber+1;
                   
                   % Check if flux is non-spatial.
                   if strcmp(fluxesAndVariablesToPlot(i),  'P') ...
                   || strcmp(fluxesAndVariablesToPlot(i), 'boundaryFlow_lower') ...
                   || strcmp(fluxesAndVariablesToPlot(i),  'boundaryFlow_upper') 
                       if j == 1
                            h = subplot(nfluxes, nNodes, [plotNumber : plotNumber+nNodes-1]);                        
                            % Plot flux                        
                            plot(h, timePoints, fluxesAtNodes.(char(fluxesAndVariablesToPlot{i})), '-k');                                                                                              
                       else
                           continue
                       end                                           
                       
                   elseif strcmp(fluxesAndVariablesToPlot{i}, 'massBalance')
                       if j == 1;
                           h = subplot(nfluxes, nNodes, [plotNumber : plotNumber+nNodes-1]);    
                           plot(h, timePoints(2:end), fluxesAtNodes.(char(fluxesAndVariablesToPlot{i})), '-k');                        
                           hold on;
                       else
                           continue                            
                       end                              
                   elseif isfield(fluxesAndVariables, fluxesAndVariablesToPlot(i))
                        % Create sub-plot                        
                        h = subplot(nfluxes, nNodes, plotNumber);
                        
                        % Plot flux  
                        plot(h, timePoints, fluxesAtNodes.(char(fluxesAndVariablesToPlot{i}))(j,:), '-k');                        
                        hold on;
                   else
                        continue;
                   end
                   
                    % Get y-axis label
                    ylabel_str='(No Label)';
                    switch char(fluxesAndVariablesToPlot{i})
                        case 'P'; ylabel_str = 'Precipitation (mm/day)';
                        case 'intercept'; ylabel_str = 'Interception (mm/day)';
                        case 'infilt'; ylabel_str = 'Infiltration (mm/day)';
                        case 'satRunoff'; ylabel_str = 'Runoff (mm/day)';
                        case 'evap'; ylabel_str = 'Soil moisture evap. (mm/day)';
                        case 'transp'; ylabel_str = 'Transpiration (mm/day)';
                        case 'recharge'; ylabel_str = 'Recharge (mm/day)';
                        case 'headChange_masstoSoil'; ylabel_str = 'Head derived soil moisture input (mm/day)';
                        case 'evap_gw'; ylabel_str = 'Groundwater evap. (mm/day)';
                        case 'Q_toStream'; ylabel_str = 'Baseflow (m^2/day)';    
                        case 'boundaryFlow_lower'; ylabel_str = 'Lower boundary sat. flow (m^3/day)';    
                        case 'boundaryFlow_upper'; ylabel_str = 'Upper boundary sat. flow (m^3/day)';                                
                        case 'mu'; ylabel_str = 'Bulk specific yield, \mu';                                
                        case 'b'; ylabel_str = 'Sat. thickness perp. to basement (m)';
                        case 'b_vert'; ylabel_str = 'Sat. thickness (m)';
                        case 'h'; ylabel_str = 'Head perp. to basement (m)';
                        case 'h_vert'; ylabel_str = 'Head (m)';
                        case 'DBNS'; ylabel_str = 'Depth to water table (m)';
                        case 'k'; ylabel_str = 'Sat. conductivity (m/day)';
                        case 'unsatDepth'; ylabel_str = 'Unsat. thickness (m)';                            
                        case 'theta'; ylabel_str = 'Vol. soil moisture (\theta_{v/v})';
                        case 'theta_frac'; ylabel_str = 'Fractional soil moisture';
                        case 'theta_frac_veg'; ylabel_str = 'Fractional veg. soil moisture';                            
                        case 'LAI'; ylabel_str = 'LAI';                                
                        case 'LAIdecay'; ylabel_str = 'LAI (as fraction of potential)'; 
                        case 'k_unsat'; ylabel_str = 'Unsat. vert. conductivity (mm/day)';
                        case 'psi_matric'; ylabel_str = 'Matrix potential (mm)';                             
                        case 'isAquiferInSoil'; ylabel_str = 'Index for aquifer within soil';                             
                        case 'dpsidz'; ylabel_str = 'd \psi dz';
                        case 'massBalance'; ylabel_str = 'Catchment mass balance (mm/day)';    
                        case 'S'; ylabel_str = 'Groundwater storage (m^2)';    
                        case 'm'; ylabel_str = 'Soil moisture (mm)';    
                    end    

                    % Label sub-plot
                    box on;
                    xlabel('Time (years)');
                    ylabel(ylabel_str);
                    title(['Flux at ',num2str( nodesToPlot(j)), ' metres']);
                                           
                    % reverse y-axis if DBNS
                    if strcmp(char(fluxesAndVariablesToPlot{i}),'DBNS');
                        set(gca,'YDir','reverse');
                    end
               end
           end
           
           function I = simp_extend( x, y )
           %SIMP_EXTEND Extended simpson's integration
           %   Taken from Numerical Recipes in fortran, section p.1 page 128, 1992
    
               dx = x(2:end,:) - x(1:end-1,:);
               nx = size(x,1);
               h=min(dx);

               if abs(min(dx) - max(dx))>1e-6*dx;
                   error('Simpsons integration requires a uniform grid spacing');
               else
                   dx = min(dx);
                   % Do simpson's extended composite integration
                   if mod(nx,2) == 0;

                       % Doing simpson's extended integration for all points but last
                       % interval. For the last interval simple trapazidal integration
                       % is done
                       nx =size(y,1)-1;
                       I = dx./3*(y(1,:) + 4.*sum(y(2:2:nx-1,:),1) + 2.*sum(y(3:2:nx-2,:),1) + y(end,:));    
                       I = I + trapz( x(end-1:end,:), y(end-1:end,:), 1);

                   else
                   % Doing simpson's extended integration for all points
                   I = dx./3*(y(1,:) + 4.*sum(y(2:2:nx-1,:),1) + 2.*sum(y(3:2:nx-2,:),1) + y(end,:));    
                   end
               end
           end
       end              

%% Estimation of attractors via time-integration.
       function [DBNS_shallowIC, Params_shallowIC, DBNS_deepIC, Params_deepIC]= plotAttractors( obj, ta, tb, nForcingSubSteps, ParameterName, Parameter_values, node_locations_to_plot )
% PLOTATTRACTORS Calculates and plots estimates of model attractors.
%
% Syntax:
%   [DBNS_shallowIC, Params_shallowIC, DBNS_deepIC, Params_deepIC]=
%       plotAttractors( obj, ta, tb, nForcingSubSteps, ParameterName, Parameter_values, node_locations_to_plot )
%
% Description:
%   This method undertakes time-integration solutions of the input model
%   from a very deep and shallow groundwater initial condition at all values
%   of a user specified model parameter. This provides a means for
%   identifying multiple attractors and their ground water level.
%   However, estimation of a repellor between multiple attractors can only
%   be undertaken by continuation analysis. 
%
%   Importantly, this method requires the climate forcing and solver
%   options to be already defined within the input model object, obj. The
%   climate data CANNOT be stochastic. Also, the solution duration
%   must be sufficient to ensure the model has converged to a steady
%   state.
%
%   Additionally, the solution time for this method can be extremely long.
%   If your computer has multiple cores, this can be reduced by turning on 
%   MatLab parallel processing as follows, where 'nCPUS' is the number
%   of CPUs in your computer: 'matlabpool nCPUS' 
%   fluxesAndVariablesToPlot - cell vector of strings containing the
%   variables, fluxes rtc to plot.
%
% Inputs:
%   obj - object of type 'CatchmentResilienceModel'.
%
%   ta - scalar integer of the start year for simulation.
%
%   tb - scalar integer of the end year for simulation.
%
%   nForcingSubSteps - scalar integer >0 defining the number of sub-steps
%   from that of the forcing data to report. For example, if the forcing data
%   is observed monthly data and 'nForcingSubSteps=4' then solutions will
%   be reported approximately every week. If the forcing rate is defined by a
%   smooth function (e.g. a spline or sine function) then 'nForcingSubSteps'
%   is the number of simulations per year. This feature was included to
%   facilitate limit-cycle continuation.
%
%   ParameterName - string of the model parameter to be explore for
%   multiple. It must be specified as the fully. For example 
%   'aquifer.Units{1,1}.k' for the aquifer saturated lateral hydraulic 
%   conductivity within the first hydrogeological unit.
%   
%   Parameter_values - scalar vector of the model parameter values at which
%   model solutions should be derived.
%
%   nodesToPlot - vector of the model spatial nodes at which plots are 
%   required. This is input as a distance, not model node number. 
%
% Outputs:
%   DBNS_shallowIC - nNodes x nTimePoints x nParameterValues matrix of
%   depth to water table model solutions from the shallow water table
%   initial condition.
%
%   Params_shallowIC - vector of parameter values at which model solutions
%   from the shallow initial condition were successfully derived.
%
%   DBNS_deepIC - nNodes x nTimePoints x nParameterValues matrix of
%   depth to water table model solutions from the deep water table
%   initial condition.
%
%   DBNS_deepIC - vector of parameter values at which model solutions
%   from the deep initial condition were successfully derived.
%
% See also:
%   CatchmentResilienceModel: class_description;
%
% Dependencies
%   (none)
%
% Author: 
%   Dr. Tim Peterson, The Department of Infrastructure
%   Engineering, The University of Melbourne.
%
% Date:
%   2 August 2012        
        
           % Get number of parameter values to investigate
           nParameterValues = length(Parameter_values);        
                      
           % Get original parameter value
           originalParameterValue = getParameterValue(obj,ParameterName);
                                 
           % Check model nodes exist at requested plot locations
           u = obj.modelParameters.geometry.u_nodes;
           for i=1: length(node_locations_to_plot)
           
                % Create cell node filter to extract location.
                u_filter = u==node_locations_to_plot(i);
                if ~sum(u_filter)
                    error(['Model node does not exist at ', num2str(node_locations_to_plot(i))]);
                end
           end
           
           % Test run to get the number of time points within the last year
           % of simulation.
           %------------
           setInitialConditions(obj, [0.5, 0.4],'fractional');
               
           % Run model
           solveModel(obj, max(ta, tb-2), tb, nForcingSubSteps);
           
           % Get filter for the last year of simulation data
           t_filter = find(floor(obj.results.t) == floor(obj.results.t(end)-sqrt(eps())));
           
           % Add last point 
           if t_filter(end) < length( obj.results.t )
               t_filter = [t_filter, length( obj.results.t )];
           end         
           t_filter = length( obj.results.t ) - t_filter;
           %------------
           
           % Initalise matrix.
           DBNS_shallowIC = nan(size(u,1), length(t_filter), nParameterValues);
           DBNS_deepIC = nan(size(u,1), length(t_filter), nParameterValues);
           shallowIC_modelFailure = false( 1, nParameterValues);
           deepIC_modelFailure = false( 1, nParameterValues);
           display('... Starting simulations for detection of attractors: ');
           parfor i=1: nParameterValues

                % Set model parameter
                setParameterValue(obj,ParameterName, Parameter_values(i));                 

                % Set shallow initial condition
                setInitialConditions(obj, [0.9, 0.4],'fractional');

                try
                    % Run model
                    solveModel(obj, ta, tb, nForcingSubSteps);
                    
                    % Turn on BSXOPS
                    bsxops(1);          
                    
                    % Extract final solutions and calculate depth to water table by
                    % calling setVariables().  
                    t_filter_temp = length(obj.results.t) - t_filter;
                    setVariables(obj, obj.results.S(:, t_filter_temp), obj.results.m(:, t_filter_temp), obj.results.t(t_filter_temp));
                    DBNS_shallowIC(:,:,i) =  getVariables(obj, 'DBNS');            
                catch
                     shallowIC_modelFailure(i)=true;
                end

                % Set deep initial condition
                setInitialConditions(obj, [0.2, 0.4],'fractional');

                % Run model
                try
                    solveModel(obj, ta, tb, nForcingSubSteps);
                    
                    % Turn on BSXOPS
                    bsxops(1);
                    
                    % Extract final solutions and calculate depth to water table by
                    % calling setVariables().                    
                    t_filter_temp = length(obj.results.t) - t_filter;
                    setVariables(obj, obj.results.S(:, t_filter_temp), obj.results.m(:, t_filter_temp), obj.results.t(t_filter_temp));
                    DBNS_deepIC(:,:,i) =  getVariables(obj, 'DBNS');        
                catch
                    deepIC_modelFailure(i)=true;
                end
                
                display(['... Finished iteration ',num2str(i), ' of ' , num2str(nParameterValues)]);
           end
           
           % Add original patameter value back into the model.
           setParameterValue(obj,ParameterName, originalParameterValue );   
           
           % Turn off BSXOPS
           bsxops(0);

           % Remove simulations that resulted in errors
           DBNS_shallowIC = DBNS_shallowIC(:, :, ~shallowIC_modelFailure);
           Params_shallowIC = Parameter_values(  ~shallowIC_modelFailure );
           
           DBNS_deepIC = DBNS_deepIC(:, :, ~deepIC_modelFailure);
           Params_deepIC = Parameter_values( ~deepIC_modelFailure );
           
           % Add results to object
           obj.results.attractors.parameter_Name = ParameterName;
           obj.results.attractors.shallowIC_parameterValues = Params_shallowIC;
           obj.results.attractors.shallowIC_DBNS = DBNS_shallowIC;
           obj.results.attractors.deepIC_parameterValues = Params_deepIC;
           obj.results.attractors.deepIC_DBNS = DBNS_deepIC;           
           
           % plot results
           nplots = length(node_locations_to_plot);
           ncols=2;
           nrows = ceil(nplots/2);
           figure();                   
           for j=1:nplots
                % Create cell node filter to extract location.
                u_filter = u==node_locations_to_plot(j);

                % re-shape data to 2 dimensions.
                DBNS_shallowIC_temp = reshape( DBNS_shallowIC( u_filter, :, :), size(DBNS_shallowIC,2), size(DBNS_shallowIC,3));
                DBNS_deepIC_temp = reshape( DBNS_deepIC( u_filter, :, :), size(DBNS_deepIC,2), size(DBNS_deepIC,3));
                
                subplot( nrows,ncols,j); hold on;
                scatter( reshape(repmat(Params_shallowIC, size(DBNS_shallowIC_temp,1), 1),[],1), reshape(DBNS_shallowIC_temp,[],1) , 'xk'); 
                scatter( reshape(repmat(Params_deepIC, size(DBNS_deepIC_temp,1), 1),[],1), reshape(DBNS_deepIC_temp,[],1), 'ok'); 
                ylabel('Depth to water table (m)');
                xlabel(ParameterName);
                title([' Time integration attractor estimates at: ', num2str(node_locations_to_plot(j))]);
                hold off;            
           end

       end       

%% The following methods are for internal use only.  
%  Very minimal documentation is provided.
%--------------------------------------------------------------------------       
       % Calculate and return state variable rates of change
       function dydt = getDerivatives_DAE(t,y,varargin)
           
           % Set preliminary variables and fluxes
           % ----------------------------
           obj = varargin{1};
           
           % Extract state variables from vector, u, and set to object, obj.
           npde=2;
           obj.stateVariables.S = y( 1:npde:npde * obj.geometryVariables.nx-1 , :);
           obj.stateVariables.m = y( 2:npde:npde * obj.geometryVariables.nx , :);
           obj.stateVariables.t = t;
           
           if size(obj.stateVariables.S,2) >1 || obj.climateData.doClimateScaling;
               bsxops(1);
           else
               bsxops(0);
           end         
           
           % Set fluxes and, within method, set variables
           setFluxes(obj);  
           
           % Check if sat thickness is <0 or above land surface
           if obj.solverOptions.showWarnings
               if any(any(y(1:2:end,:) < 0))
                   display(['... WARNING: Aquifer storage variable(s) is <0 at time =', num2str(t)]);                
               elseif any(any(obj.variables.DBNS <= 0))
                   display(['... WARNING: Depth to water table is <0 at time =', num2str(t)]);            
               elseif any(any(y(2:2:end,:) < 0))
                   display(['... WARNING: Soil moisture storage variable(s) is <0 at time =', num2str(t)]);         
               elseif any(any(isnan(y(1:2:end,:))))
                   display(['... WARNING: Aquifer storage variables(s) equals NaN at time =', num2str(t)]);                    
               elseif any(any(isnan(y(2:2:end,:))))
                   display(['... WARNING: Soil moisture storage variables(s) equals NaN at time =', num2str(t)]);                        
               end
           end
           %----------------------------                      
           
           % Calculate saturated storage rate of change, dS/dt
           %----------------------------
           dbdx = diff( obj.variables.b)./diff(obj.geometryVariables.x);   
           
           branchConductance = 0.5.*(obj.modelParameters.geometry.width(1:end-1,:) + obj.modelParameters.geometry.width(2:end,:))...
                          .*0.5.*(obj.variables.k(2:end,:) .* obj.variables.b(2:end,:) + obj.variables.k(1:end-1,:) .* obj.variables.b(1:end-1,:));
           S_lat = branchConductance.*( obj.geometryVariables.cos_i_avg .* dbdx +  obj.geometryVariables.sin_i_avg);   

           % Est lateral flux component of dSdt at non boundary nodes
           S_index = 1:2:2*obj.geometryVariables.nx-1;
           dydt =  zeros(npde.*obj.geometryVariables.nx , size(y,2)); 
           dydt( S_index(2:end-1) ,:) = (S_lat(2:end,:) - S_lat(1:end-1,:))./obj.geometryVariables.x_midNodes;

           % Est lateral flux component of dSdt at boundary nodes
           dydt(1,:) = (S_lat(1,:) - obj.fluxes.boundaryFlow_lower)./ (0.5.*(obj.geometryVariables.x(2) - obj.geometryVariables.x(1)));
           dydt(end-1,:) = (obj.fluxes.boundaryFlow_upper - S_lat(end,:))./ (0.5.*(obj.geometryVariables.x(end) - obj.geometryVariables.x(end-1)));
          
           % Calc source term at x (365.25e-3 because convertted from days to years and mm to metres)  
           N = 365.25e-3 .* obj.modelParameters.geometry.width ... 
               .*( obj.geometryVariables.cos_BedAngle .* obj.fluxes.recharge - obj.fluxes.evap_gw) ...
                + obj.fluxes.Q_toStream; % - obj.fluxes.satDischarge;

           % Calc denominator for S at boundary nodes (inverted and called SMscaler)                                 
           if obj.solverOptions.useDAE                          
               SMscaler = 1./(1 -  obj.modelParameters.aquifer.singularityFraction .* obj.variables.isAquiferInSoil .* obj.variables.theta ./ obj.geometryVariables.theta_s );
               dydt( S_index,:) = dydt( S_index ,:) + N;                                             
           else
               SMscaler = 1./(1 -  obj.modelParameters.aquifer.singularityFraction .* obj.variables.isAquiferInSoil .* obj.variables.theta ./ obj.geometryVariables.theta_s );
               % Add source term and multiple sum by scaler for accelerated
               % head rise within the vadose layer.            
               dydt( S_index,:) = SMscaler.*(dydt( S_index ,:) + N);               
           end                      
           
           %  Calculate change in soil water due to head change.           
           obj.fluxes.headChange_masstoSoil = -1e3./(365.25.*obj.modelParameters.geometry.width) ...
                       .* obj.variables.isAquiferInSoil .* obj.variables.theta ./ obj.geometryVariables.theta_s .* dydt( S_index ,:);
           if obj.solverOptions.useDAE           
               obj.fluxes.headChange_masstoSoil = SMscaler .* obj.fluxes.headChange_masstoSoil;
           end
           
           % Calculate soil moisture rate of change, dm/dt
           dydt(2:2:2*obj.geometryVariables.nx,:) = 365.25.*( obj.fluxes.infilt - obj.fluxes.evap - obj.fluxes.transp - obj.fluxes.recharge ...
               + obj.fluxes.headChange_masstoSoil);               
           
           if obj.solverOptions.showWarnings
               if any(any(isnan(dydt(1:2:end,:))))
                   display(['... WARNING: d Aquifer storage / dt equals NaN at time =', num2str(t)]);                    
               elseif any(any(isnan(dydt(2:2:end,:))))
                   display(['... WARNING: d Soil moisture storage / dt equals NaN at time =', num2str(t)]);                        
               end
           end
           if size(obj.stateVariables.S,2) >1;
               bsxops(0);
           end
       end
       
       % Get the mass matrix for DAE model solvers.
       function massMatrix = getDerivatives_DAEmass(t,y,varargin)
           
           % Set preliminary variables and fluxes
           % ----------------------------
           obj = varargin{1};
           
           % Extract state variables from vector, u, and set to object, obj.
           npde=2;
           obj.stateVariables.S = y( 1:npde:npde * obj.geometryVariables.nx-1 , :);
           obj.stateVariables.m = y( 2:npde:npde * obj.geometryVariables.nx , :);
           obj.stateVariables.t = t;
           
           if size(obj.stateVariables.S,2) >1;
               bsxops(1);
           else
               bsxops(0);
           end
           
           if obj.solverOptions.showWarnings
               if any(any(y(1:2:end,:) < 0))
                   display(['... WARNING: Aquifer storage variable(s) is <0 at time =', num2str(t)]);                
               elseif any(any(y(2:2:end,:) < 0))
                   display(['... WARNING: Soil moisture storage variable(s) is <0 at time =', num2str(t)]);         
               elseif any(any(isnan(y(1:2:end,:))))
                   display(['... WARNING: Aquifer storage variables(s) equals NaN at time =', num2str(t)]);                    
               elseif any(any(isnan(y(2:2:end,:))))
                   display(['... WARNING: Soil moisture storage variables(s) equals NaN at time =', num2str(t)]);                        
               end
           end
           % Set fluxes and, within method, set variables
           setVariables(obj);  
           
           soilMoistureScaler = 1 -  obj.modelParameters.aquifer.singularityFraction .* obj.variables.isAquiferInSoil .* obj.variables.theta ./ obj.geometryVariables.theta_s;           
           
           massMatrix = ones(size(y));
           massMatrix(1:2:end,:) = soilMoistureScaler;
           massMatrix = spdiags(massMatrix,0, size(massMatrix,1), size(massMatrix,1));
       end       
       
       % Set the user specified model parameters.
       function setParameterValue(obj,ParameterName, ParameterValue)
       
           % Divide ParameterName into structure names
           fieldDots = strfind(ParameterName,'.');
           nfield = length(fieldDots)+1;                      
           for i=1: nfield
              if i==1 
                fieldName{i} = ParameterName(1:fieldDots(i)-1);
              elseif i==nfield
                  fieldName{i} = ParameterName(fieldDots(i-1)+1:end);
              else
                fieldName{i} = ParameterName(fieldDots(i-1)+1:fieldDots(i)-1);                  
              end
           end
           
           % Check if second field is a cell structure
           hasUnitField = false;
           if nfield==3
               if strcmp(fieldName{2}(1:5),'Units')
                   hasUnitField=true;
                   cellLBracket = strfind(fieldName{2},'{');
                   cellRBracket = strfind(fieldName{2},'}');
                   cellComma = strfind(fieldName{2},',');
                   if ~isempty(cellComma)
                       unitNumber1 = str2double(char(fieldName{2}(cellLBracket+1 : cellComma-1)));
                       unitNumber2 = str2double(char(fieldName{2}(cellComma +1: cellRBracket-1 )));
                       if unitNumber1 ~=1 && unitNumber2~=1 
                           error(' "ParameterName" sent to the function "setParameterValue" contains an incorrect index for "Units" field!');
                       else
                           unitNumber = max(unitNumber1, unitNumber2);
                       end
                   else
                       unitNumber = str2double(char(fieldName{2}(cellLBracket+1 : cellRBracket-1)));
                   end
                   fieldName{2} = fieldName{2}(1:5);
               end
           end
           
           % Assign field value.
           try
               switch nfield
                   case 1
                       obj.modelParameters.(char(fieldName{1})) = ParameterValue;
                   case 2
                       obj.modelParameters.(char(fieldName{1})).(char(fieldName{2})) = ParameterValue;
                   case 3
                       if hasUnitField
                           obj.modelParameters.(char(fieldName{1})).(char(fieldName{2})){unitNumber}.(char(fieldName{3})) = ParameterValue;
                       else    
                           obj.modelParameters.(char(fieldName{1})).(char(fieldName{2})).(char(fieldName{3})) = ParameterValue;
                       end
                   otherwise
                       error(' "ParameterName" sent to the function "setParameterValue" should not contain more than three field names.');
               end
           catch
               error(' "ParameterName" sent to the function "setParameterValue" contains a field name or unit number that is incorrect.');
           end                              
       end
       
       % Get values for the user specified parameter(s).
       function ParameterValue = getParameterValue(obj,ParameterName)
                      % Divide ParameterName into structure names
           fieldDots = strfind(ParameterName,'.');
           nfield = length(fieldDots)+1;                      
           for i=1: nfield
              if i==1 
                fieldName{i} = ParameterName(1:fieldDots(i)-1);
              elseif i==nfield
                  fieldName{i} = ParameterName(fieldDots(i-1)+1:end);
              else
                fieldName{i} = ParameterName(fieldDots(i-1)+1:fieldDots(i)-1);                  
              end
           end
           
           % Check if second field is a cell structure
           hasUnitField = false;
           if nfield==3
               if strcmp(fieldName{2}(1:5),'Units')
                   hasUnitField=true;
                   cellLBracket = strfind(fieldName{2},'{');
                   cellRBracket = strfind(fieldName{2},'}');
                   cellComma = strfind(fieldName{2},',');
                   if ~isempty(cellComma)
                       unitNumber1 = str2double(char(fieldName{2}(cellLBracket+1 : cellComma-1)));
                       unitNumber2 = str2double(char(fieldName{2}(cellComma +1: cellRBracket-1 )));
                       if unitNumber1 ~=1 && unitNumber2~=1 
                           error(' "ParameterName" sent to the function "getParameterValue" contains an incorrect index for "Units" field!');
                       else
                           unitNumber = max(unitNumber1, unitNumber2);
                       end
                   else
                       unitNumber = str2double(char(fieldName{2}(cellLBracket+1 : cellRBracket-1)));
                   end
                   fieldName{2} = fieldName{2}(1:5);
               end
           end
           
           % Assign field value.
           try
               switch nfield
                   case 1
                       ParameterValue = obj.modelParameters.(char(fieldName{1}));
                   case 2
                       ParameterValue = obj.modelParameters.(char(fieldName{1})).(char(fieldName{2}));
                   case 3
                       if hasUnitField
                           ParameterValue = obj.modelParameters.(char(fieldName{1})).(char(fieldName{2})){unitNumber}.(char(fieldName{3}));
                       else    
                           ParameterValue = obj.modelParameters.(char(fieldName{1})).(char(fieldName{2})).(char(fieldName{3}));
                       end
                   otherwise
                       error(' "ParameterName" sent to the function "getParameterValue" should not contain more than three field names.');
               end
           catch
               error(' "ParameterName" sent to the function "getParameterValue" contains a field name or unit number that is incorrect.');
           end                              
       end
       
       % Set the geometry variables for the model.
       function setGeometryVariables(obj)
           
           % Convert horizontal plane nodes,u, into model nodes along the
           % aquifer basement. This is done by integrating to find x-plane nodes corresponding to u-plane 
           % (horizontal) nodes.
           %--------------------------
           u = obj.modelParameters.geometry.u_nodes;
           obj.geometryVariables.x = zeros(size(u));
           obj.geometryVariables.x(1) = u(1);
           for i=2:size(u,1)
               obj.geometryVariables.x(i) = obj.geometryVariables.x(i-1) + sqrt( (u(i)-u(i-1))^2 ...
                   + (obj.modelParameters.geometry.basementElevation(i) ....
                   - obj.modelParameters.geometry.basementElevation(i-1))^2);
           end
           
           % number of 'x' nodes
           obj.geometryVariables.nx = size(obj.geometryVariables.x,1);
           
           % Additional mid-node points.
           obj.geometryVariables.x_midNodes = diff(obj.geometryVariables.x(1:end-1,:) ...
               + 0.5.*diff(obj.geometryVariables.x));           
           %--------------------------
           
           % Assign class of aquifer, vadose zone and vegetation landcover
           % types to model nodes. Note, the vegetation cover can chnage
           % with time and is as such of a slightly different structure.
           % See below for details.
           %--------------------------
           nu = size(u,1);
           
           % Aquifer types
           if size(obj.modelParameters.aquifer.Units,1)==1;
                   obj.geometryVariables.k = obj.modelParameters.aquifer.Units{1,1}.k;
                   obj.geometryVariables.f = obj.modelParameters.aquifer.Units{1,1}.f;
                   obj.geometryVariables.d_khalf = obj.modelParameters.aquifer.Units{1,1}.d_khalf;
                   obj.geometryVariables.tor = obj.modelParameters.aquifer.Units{1,1}.tor;               
           else
               obj.geometryVariables.k= zeros(nu,1);
               obj.geometryVariables.f= zeros(nu,1);
               obj.geometryVariables.d_khalf= zeros(nu,1);
               obj.geometryVariables.tor= zeros(nu,1);

               for i=1:size(obj.modelParameters.aquifer.Units,1);
                   obj.geometryVariables.k(obj.modelParameters.aquifer.UnitCoverage(:,2)==i,1) ...
                       = obj.modelParameters.aquifer.Units{i,1}.k;
                   obj.geometryVariables.f(obj.modelParameters.aquifer.UnitCoverage(:,2)==i,1) ...
                       = obj.modelParameters.aquifer.Units{i,1}.f;
                   obj.geometryVariables.d_khalf(obj.modelParameters.aquifer.UnitCoverage(:,2)==i,1)  ...
                       = obj.modelParameters.aquifer.Units{i,1}.d_khalf;
                   obj.geometryVariables.tor(obj.modelParameters.aquifer.UnitCoverage(:,2)==i,1) ...
                       = obj.modelParameters.aquifer.Units{i,1}.tor;
               end
           end
           
           % Vadose zone
           if size(obj.modelParameters.vadose.Units,1)==1
               obj.geometryVariables.theta_s = obj.modelParameters.vadose.Units{1,1}.theta_s;
               obj.geometryVariables.theta_r = obj.modelParameters.vadose.Units{1,1}.theta_r;
               obj.geometryVariables.BC_phi = obj.modelParameters.vadose.Units{1,1}.BC_phi;
               obj.geometryVariables.BC_psi_a = obj.modelParameters.vadose.Units{1,1}.BC_psi_a;
               obj.geometryVariables.k_vert = obj.modelParameters.vadose.Units{1,1}.k_vert;
               obj.geometryVariables.d_evap = obj.modelParameters.vadose.Units{1,1}.d_evap;               
               obj.geometryVariables.Io = obj.modelParameters.vadose.Units{1,1}.Io;               
           else
               obj.geometryVariables.theta_s = zeros(nu,1);
               obj.geometryVariables.theta_r = zeros(nu,1);
               obj.geometryVariables.BC_phi = zeros(nu,1);
               obj.geometryVariables.BC_psi_a = zeros(nu,1);
               obj.geometryVariables.k_vert = zeros(nu,1);
               obj.geometryVariables.d_evap = zeros(nu,1);        
               obj.geometryVariables.Io = zeros(nu,1);        
               for i=1:size(obj.modelParameters.vadose.Units,1);
                   obj.geometryVariables.theta_s(obj.modelParameters.vadose.UnitCoverage(:,2)==i,1) = obj.modelParameters.vadose.Units{i,1}.theta_s;
                   obj.geometryVariables.theta_r(obj.modelParameters.vadose.UnitCoverage(:,2)==i,1) = obj.modelParameters.vadose.Units{i,1}.theta_r;
                   obj.geometryVariables.BC_phi(obj.modelParameters.vadose.UnitCoverage(:,2)==i,1) = obj.modelParameters.vadose.Units{i,1}.BC_phi;
                   obj.geometryVariables.BC_psi_a(obj.modelParameters.vadose.UnitCoverage(:,2)==i,1) = obj.modelParameters.vadose.Units{i,1}.BC_psi_a;
                   obj.geometryVariables.k_vert(obj.modelParameters.vadose.UnitCoverage(:,2)==i,1) = obj.modelParameters.vadose.Units{i,1}.k_vert;
                   obj.geometryVariables.d_evap(obj.modelParameters.vadose.UnitCoverage(:,2)==i,1) = obj.modelParameters.vadose.Units{i,1}.d_evap;
                   obj.geometryVariables.Io(obj.modelParameters.vadose.UnitCoverage(:,2)==i,1) = obj.modelParameters.vadose.Units{i,1}.Io;
               end 
           end
           
           % Vegetation.
           % Note, the vegetation cover can chnage with time. As such
           % obj.modelParameters.veg.UnitCoverage is of a different
           % structure to that for vadose zone and aquifer detailed above.
           % The matrix 'veg.UnitCoverage' has the following structure: 
           %
           % - first row is the start time (in years) of the vegetation cover. 
           % - rows 2 to n nodes+1 are the vegetation cover types for each
           %   node. 
           % - first column is the node locatons in metres from the
           %   catchment outlet. 
           % - columns 2 to nVegTimePoints+1 are the vegetation cover for
           %   each time point detailed in the first row.
           if size(obj.modelParameters.veg.Units,1)==1;
               obj.geometryVariables.F_canopy = obj.modelParameters.veg.Units{1,1}.F_canopy; 
               obj.geometryVariables.k_light = obj.modelParameters.veg.Units{1,1}.k_light; 
               obj.geometryVariables.wiltingPoint_MPa = obj.modelParameters.veg.Units{1,1}.wiltingPoint_MPa; 
               obj.geometryVariables.stomataClosure_MPa = obj.modelParameters.veg.Units{1,1}.stomataClosure_MPa; 
               obj.geometryVariables.d_0p5LAI = obj.modelParameters.veg.Units{1,1}.d_0p5LAI;
               obj.geometryVariables.alpha = obj.modelParameters.veg.Units{1,1}.alpha;
               
           else
               nVegTimePoints = size(obj.modelParameters.vadose.UnitCoverage,2)-1;
               obj.geometryVariables.F_canopy = zeros(nu,nVegTimePoints);
               obj.geometryVariables.k_light = zeros(nu,nVegTimePoints);           
               obj.geometryVariables.wiltingPoint_MPa = zeros(nu,nVegTimePoints);
               obj.geometryVariables.stomataClosure_MPa = zeros(nu,nVegTimePoints);               
               obj.geometryVariables.d_0p5LAI = zeros(nu,nVegTimePoints);
               obj.geometryVariables.alpha = zeros(nu,nVegTimePoints);               
               obj.geometryVariables.LAI_avg = zeros(nu, 12, nVegTimePoints);

               for i=1:size(obj.modelParameters.veg.Units,1);
                   for j=1:nVegTimePoints;
                       obj.geometryVariables.F_canopy(obj.modelParameters.veg.UnitCoverage(2:end, j+1)==i,j) = obj.modelParameters.veg.Units{i,1}.F_canopy; 
                       obj.geometryVariables.k_light(obj.modelParameters.veg.UnitCoverage(2:end, j+1)==i,j) = obj.modelParameters.veg.Units{i,1}.k_light; 
                       obj.geometryVariables.wiltingPoint_MPa(obj.modelParameters.veg.UnitCoverage(2:end, j+1)==i,j) = obj.modelParameters.veg.Units{i,1}.wiltingPoint_MPa; 
                       obj.geometryVariables.stomataClosure_MPa(obj.modelParameters.veg.UnitCoverage(2:end, j+1)==i,j) = obj.modelParameters.veg.Units{i,1}.stomataClosure_MPa; 
                       obj.geometryVariables.d_0p5LAI(obj.modelParameters.veg.UnitCoverage(2:end, j+1)==i,j) = obj.modelParameters.veg.Units{i,1}.d_0p5LAI;
                       obj.geometryVariables.alpha(obj.modelParameters.veg.UnitCoverage(2:end, j+1)==i,j) = obj.modelParameters.veg.Units{i,1}.alpha;
                       if ~iscell(obj.modelParameters.veg.Units{i,1}.LAI_avg);
                           for k=1:12;
                              obj.geometryVariables.LAI_avg(obj.modelParameters.veg.UnitCoverage(2:end, j+1)==i,k,j) = obj.modelParameters.veg.Units{i,1}.LAI_avg(k);                               
                           end
                       end
                   end
               end
               if ~iscell(obj.modelParameters.veg.Units{i,1}.LAI_avg);
                   obj.geometryVariables.LAI_avg = [obj.geometryVariables.LAI_avg, obj.geometryVariables.LAI_avg(:,1)];
               end
           end
           %--------------------------
           
           % Calculate derived constants 
           %--------------------------           
           
           % Derived depths.
           obj.geometryVariables.streamElevation = obj.modelParameters.geometry.surfaceElevation - obj.modelParameters.geometry.streamDepth;
           obj.geometryVariables.basementDepth = obj.modelParameters.geometry.surfaceElevation - obj.modelParameters.geometry.basementElevation;            

           % Bed angle.
           dBeddu = (obj.modelParameters.geometry.basementElevation(3:end,:) - obj.modelParameters.geometry.basementElevation(1:end-2,:) ) ...
               ./( u(3:end,:) - u(1:end-2,:) );
           dBeddu = [ (obj.modelParameters.geometry.basementElevation(2,:) - obj.modelParameters.geometry.basementElevation(1,:) ) ...
               ./(u(2) - u(1)) ; dBeddu;  ...
               (obj.modelParameters.geometry.basementElevation(end) - obj.modelParameters.geometry.basementElevation(end-1) ) ...
               ./(u(end) - u(end-1))];
           obj.geometryVariables.bedAngle = atan(dBeddu);
           obj.geometryVariables.cos_BedAngle = cos(obj.geometryVariables.bedAngle);
           obj.geometryVariables.sin_BedAngle = sin(obj.geometryVariables.bedAngle);
                      
           % Bed angle between model nodes.
           obj.geometryVariables.cos_i_avg = 0.5.*(obj.geometryVariables.cos_BedAngle(2:end,:) ...
               + obj.geometryVariables.cos_BedAngle(1:end-1,:));
           obj.geometryVariables.sin_i_avg = 0.5.*(obj.geometryVariables.sin_BedAngle(2:end,:) ...
               + obj.geometryVariables.sin_BedAngle(1:end-1,:));                                        
           
           % Aquifer maxium capacity below soil layer, S_aq.
           obj.geometryVariables.aquiferMaxStorage = obj.modelParameters.geometry.width .* obj.geometryVariables.f ...
               .* (obj.geometryVariables.basementDepth - obj.modelParameters.geometry.soilDepth) .* obj.geometryVariables.cos_BedAngle;                      
           %--------------------------
                      
           % Constants for conversion of dates.
           obj.geometryVariables.daysPerMonth = [31 31	28.25	31	30	31	30	31	31	30	31	30	31];
           obj.geometryVariables.cumDaysPerMonth=[0 31	59.25	90.25	120.25	151.25	181.25	212.25	243.25	273.25	304.25	334.25	365.25]';

           % Convert wilting point and stomata closure point from PSI in
           % MPa to mm H2O ia Van Genuchten retention equations.
           %--------------------------
           obj.geometryVariables.theta_wp = theta_from_psi(obj.geometryVariables.wiltingPoint_MPa);
           obj.geometryVariables.theta_star = theta_from_psi(obj.geometryVariables.stomataClosure_MPa);
           
           function theta_result = theta_from_psi(psi_inMPa)
               psi = psi_inMPa * 4014.63075975562*25.4;
               n = obj.geometryVariables.BC_phi +1;
               m = 1 - 1./n;        
               theta_result = (1./(1 + (psi./obj.geometryVariables.BC_psi_a).^n)).^m.*(obj.geometryVariables.theta_s - obj.geometryVariables.theta_r) + obj.geometryVariables.theta_r;
           end
           %--------------------------
       end
       
       % Get a user specified model variable.
       function var= getVariables(obj, variableName)           
           var = obj.variables.(variableName);           
       end
       
       % Get a user specified model flux.
       function var= getFluxes(obj, fluxName)           
           var = obj.fluxes.(fluxName);           
       end
       
       % Calculate and set transient model variables;
       function setVariables(obj, varargin)
           
           % Update state variables if not already within object.
           if nargin==4
               obj.stateVariables.S = varargin{1};
               obj.stateVariables.m= varargin{2};
               obj.stateVariables.t = varargin{3};               
           end
           
           % Extract state variables.
           S=obj.stateVariables.S;
           m=obj.stateVariables.m;
           t=obj.stateVariables.t;
           
           % Saturated zone variables
           %--------------------------
                     
           % Calculate aquifer bulk specific yield, mu
           F = obj.modelParameters.geometry.width .* obj.modelParameters.aquifer.lambda_m .* log( exp( (S - obj.geometryVariables.aquiferMaxStorage)./(obj.modelParameters.geometry.width.*obj.modelParameters.aquifer.lambda_m)) + 1);
           obj.variables.mu = S .* obj.geometryVariables.f./(S - F.*(obj.geometryVariables.theta_s - obj.geometryVariables.f)./obj.geometryVariables.theta_s);
                      
           % Convert aquifer storage to saturated thickness
           obj.variables.b = S./(obj.modelParameters.geometry.width .* obj.variables.mu);
           
           % Convert saturated thickness to vertical saturated thickness
           obj.variables.b_vert = obj.variables.b ./ obj.geometryVariables.cos_BedAngle;
           
           % Convert vertical saturated thickness to vertical head
           obj.variables.h_vert = obj.modelParameters.geometry.basementElevation + obj.variables.b_vert;
           
           % Convert head to depth to water table, DBNS;
           obj.variables.DBNS = obj.modelParameters.geometry.surfaceElevation - obj.variables.h_vert;
           
           % Saturated lateral conductivity decays exponentially from a large sat
           % thickness to zeros: k = k_sat_max*(1 - exp(-(b-d_khalf)/tor ))
           % Integrating over sat thickness and assuming k=0 at basement:
           % note the 'b' outer demoninator is needed has in calculating the
           % conductance this eqn is multiplied by b (ie within dudt.m)                                 
           C = log(exp(obj.geometryVariables.tor.*obj.geometryVariables.d_khalf)) ...
               - log(1+exp(obj.geometryVariables.tor.*obj.geometryVariables.d_khalf));
           
           obj.variables.k = 365.25 .* obj.geometryVariables.k ./(obj.geometryVariables.tor.*obj.variables.b) .* ...
               (log(1+exp((obj.geometryVariables.d_khalf -obj.variables.b) .* obj.geometryVariables.tor )) -  ...
               log(exp((obj.geometryVariables.d_khalf -obj.variables.b) .* obj.geometryVariables.tor))+C);                      
           %--------------------------

           % Un-saturated zone variables
           %--------------------------
           % Calculate unsaturated thickess 
           obj.variables.unsatDepth = obj.variables.DBNS  - (obj.modelParameters.geometry.width .* obj.modelParameters.aquifer.lambda_m) ...
               .*log(1+exp(-(S - obj.geometryVariables.aquiferMaxStorage)./ (obj.modelParameters.geometry.width .* obj.modelParameters.aquifer.lambda_m))) ...
               ./(obj.modelParameters.geometry.width .* obj.geometryVariables.f);
                      
           % Calculate if the water table is within the unsaturated layer.
           obj.variables.isAquiferInSoil = 1./(1 + exp(-(S - obj.geometryVariables.aquiferMaxStorage) ./ (obj.modelParameters.geometry.width .* obj.modelParameters.aquifer.lambda_m) ));                                
           
           % It is assume that the vadose zone is saturated below the air
           % entry point above the water table. 
           %obj.variables.unsatDepth = obj.variables.unsatDepth - 1e-3 .* obj.variables.isAquiferInSoil .* obj.geometryVariables.BC_psi_a;
           
           % Convert soil moisture depth to volumetric soil moisture (with
           % smoothing at thresholds at theta_r and theta_star.           
           obj.variables.theta = 1e-3.*m./obj.variables.unsatDepth;              
           k = 1./(obj.geometryVariables.theta_s - obj.geometryVariables.theta_r);
           
           lambda_theta = 0.0005;        
           obj.variables.theta = obj.geometryVariables.theta_r + (obj.geometryVariables.theta_s - obj.geometryVariables.theta_r) ...
               .*min(1 ,k.*lambda_theta.*log( (1 + exp((obj.variables.theta - obj.geometryVariables.theta_r) ...
               ./lambda_theta))./(1 + exp((obj.variables.theta - obj.geometryVariables.theta_s)./lambda_theta)))); 
           
           obj.variables.theta_frac_max = obj.geometryVariables.theta_r + (obj.geometryVariables.theta_s - obj.geometryVariables.theta_r) ...
               .*min(1 ,k.*lambda_theta.*log( (1 + exp((obj.geometryVariables.theta_s - obj.geometryVariables.theta_r) ...
               ./lambda_theta))./(1 + exp((obj.geometryVariables.theta_s - obj.geometryVariables.theta_s)./lambda_theta))));        
           obj.variables.theta_frac_max = (obj.variables.theta_frac_max - obj.geometryVariables.theta_r)./(obj.geometryVariables.theta_s - obj.geometryVariables.theta_r);
           
           % Convert theta_i to normalised values (NOTE: WRR 2009 uses
           % smoothing function of theta_veg).
           obj.variables.theta_frac = min(obj.variables.theta_frac_max, (obj.variables.theta - obj.geometryVariables.theta_r)./(obj.geometryVariables.theta_s - obj.geometryVariables.theta_r));
           k = 1./(obj.geometryVariables.theta_star - obj.geometryVariables.theta_wp); 
           obj.variables.theta_frac_veg = min(1, k.*lambda_theta.*log( (1 + exp((obj.variables.theta - obj.geometryVariables.theta_wp)./lambda_theta))./(1 + exp((obj.variables.theta  - obj.geometryVariables.theta_star )./lambda_theta))));
           
           % Get maximum LAI at month of current time-step  
           if iscell(obj.modelParameters.veg.Units{1,1}.LAI_avg)
               obj.variables.LAI = feval(obj.modelParameters.veg.Units{1,1}.LAI_avg{1,1}, obj.modelParameters.veg.Units{1,1}.LAI_avg{2:end},t)';
           else
               % Assumes geometryVariables.LAI_avg is a matix of nx rows
               % and 13 columns. It finds the column for the current month
               dateFilter = (obj.geometryVariables.cumDaysPerMonth/365.25 > t - floor(t) - sqrt(eps(t)) );
               
               obj.variables.LAI = zeros(size(t,2));
               for i=1:size(dateFilter,2)
                   [irow,icol]= find(dateFilter(i), 1, 'first');

                   dateFilter_temp = false(13,1);
                   dateFilter_temp(irow,1)=1;
                   obj.variables.LAI(1,i) = obj.geometryVariables.LAI_avg(:,dateFilter_temp);
               end
           end                          
           
           % If doPositiveFeedback=true then the model will have a positive
           % feedback. Importantly, this feedback results from a reduction
           % in LAI as a greater % of the root zone is intersected by a
           % saline water table.
           if obj.modelParameters.aquifer.doPositiveFeedback
               % Fraction of roots that are within the unsaturated zone
               % (depth converted from m to cm).
               obj.variables.rootFractionDry = max(0 , obj.variables.DBNS .^ obj.geometryVariables.alpha ...
                   ./( obj.variables.DBNS .^ obj.geometryVariables.alpha ...
                   + obj.geometryVariables.d_0p5LAI.^ obj.geometryVariables.alpha));
           else
                obj.variables.rootFractionDry = ones( size(obj.variables.DBNS));
           end
           
           % Scale max LAI by decay factor
           obj.variables.LAI = obj.variables.LAI .* obj.variables.rootFractionDry;
           
           % Calculate  unsaturated vertical conductivity from Van Genuchten.         
           n = obj.geometryVariables.BC_phi +1;
           m = 1 - 1./n;        
           obj.variables.k_unsat = obj.geometryVariables.k_vert .* obj.variables.theta_frac.^0.5.*(1-(1 - obj.variables.theta_frac.^(1./m)).^m).^2;
           
           % Calculate matrix suction at centre of unsaturated layer. 
           obj.variables.psi_matric = obj.geometryVariables.BC_psi_a .* (obj.variables.theta_frac.^(-1./m) -1 ).^(1./n);                    
         
           % d matrix potential/ dz. If watertable is below soil layer then dpsidz=0;
           obj.variables.dpsidz= obj.variables.isAquiferInSoil .* obj.variables.psi_matric./(0.5e3 .* obj.variables.DBNS);    
           %--------------------------           
       end
       
       % Calculate and set fluxes    
       function setFluxes(obj, varargin)
           
           % Clear Fluxes
           obj.fluxes=[];
           
           % Update state variables if not already within object.
           if nargin==5
               obj.stateVariables.S = varargin{1};
               obj.stateVariables.m= varargin{2};
               obj.stateVariables.t = varargin{3};
               climateData_timePoint = varargin{4};
           else
               % Get climate data for current time point, t
               climateData_timePoint = getClimate(obj, obj.stateVariables.t);               
           end
                      
           % Set variables based upon state variables
           setVariables(obj);
           
           % Un-saturated fluxes
           %------------------------------
           % Calc. intercept maximum.
           interceptPotential = obj.geometryVariables.F_canopy./10 .* obj.variables.LAI;    
           
           % Calc. infiltration capacity.
           SM_capacity = 1e3 .* obj.geometryVariables.theta_s ...
                        .* obj.variables.unsatDepth .* (obj.variables.theta_frac_max - obj.variables.theta_frac);                                          
                                 
           % Calc. some variables for soil evap.
           f = 0.05;
           a = -log(f)./obj.geometryVariables.d_evap;
           relDemand_unsat = max(0,(1 - exp(-a.*obj.variables.unsatDepth))./a);           
           relDemand_gw = exp(-a.*max(0,obj.variables.DBNS + 1e-3 .* obj.geometryVariables.BC_psi_a));
                   
           % Set infilt constant
           %lambda_p = 0.01;               
           lambda_p = 0.01;               
           
           % Calc fluxes dependent upon precipitation and ET.
           if obj.climateData.doClimateScaling
               P = nan( [1, size(climateData_timePoint,1), length(obj.climateData.scalingPercentiles)]);
               PET = P;               
               intercept = nan( [size(obj.stateVariables.S), length(obj.climateData.scalingPercentiles)]);
               Peff   = intercept; infilt  = intercept; satRunoff = intercept; 
               evap = intercept;   transp  = intercept; evap_gw   = intercept;
               
               npcntiles = length(obj.climateData.scalingPercentiles); 
               nTimePoints = size(climateData_timePoint,1);
               SM_capacity = repmat( SM_capacity, [1, 1 ,npcntiles] );
                              
               for i=1:nTimePoints 
                   P(1, i, :) = climateData_timePoint(i,1:npcntiles )';
                   PET(1, i, :) = climateData_timePoint(i,npcntiles + 1:end)';
                   
                   if nTimePoints==1 && size( obj.stateVariables.S,2)>1
                       colIndex = 1:size( obj.stateVariables.S,2);
                   else
                       colIndex  = i;
                   end
                   
                   % calculate interception, 'Intercept'                          
                   warning('OFF','MATLAB:divideByZero');
                   intercept_temp = interceptPotential(:, colIndex ,:) .*(1 - exp(-P(1, i  , :)./interceptPotential(:, colIndex ,:) ));
                   warning('ON','MATLAB:divideByZero');
                   intercept_temp( repmat(isnan(obj.variables.LAI(:, colIndex)) | obj.variables.LAI(:, colIndex)<=0,[1 1  npcntiles]))= 0;
                   intercept(:, colIndex ,:) = intercept_temp;
                   
                   % Calculate effective precipitation, 'Peff'
                   Peff(:, colIndex ,:) = P(1, i , :) - intercept(:, colIndex  ,:);
                  
                   % Smooth function to est max P for infiltration
                   % NOTE: The smoothening fn below caused a lot of problems. If lambda_p is
                   % > approx 1 the Pinfiltration can be greater than P.
                   % Also because of limitation in floating point precision,
                   % exp((Peff-SM_capacity)/lambda_p) can be >> abs(1e308)                    
                   infilt(:, colIndex ,:) = zeros(size(Peff(:, colIndex ,:) ));
                   if any( any( any( Peff(:, colIndex ,:)>0 )));
                     
                      P_infilt = zeros(size(Peff(:, colIndex ,:)));
                      P_infilt_filt = Peff(:,colIndex ,:)>0 & SM_capacity(:,colIndex ,:) >lambda_p;
                      Peff_temp = Peff(:,colIndex,:);                      
                      P_infilt(P_infilt_filt) = lambda_p.* Peff_temp(P_infilt_filt).*log( (1+exp((Peff_temp(P_infilt_filt) - lambda_p )./(Peff_temp(P_infilt_filt).*lambda_p) )) ./ (1+exp((Peff_temp(P_infilt_filt) - SM_capacity(P_infilt_filt) )./(Peff_temp(P_infilt_filt).*lambda_p) )) );    

                      % The above smooth infiltration potential function
                      % can result in negative estimates smaller than
                      % -eps(1) when precipitation is approx. zero.
                      P_infilt( P_infilt <0) = 0;
                      
                      if any(any(any(  isinf(P_infilt) | isnan(P_infilt) )))
                          error(['The lambda_precip infiltration smoothing parameter is too small!', char(13), 'The potential infiltration equals infinity when precipitation is below the soil moisture deficit.']);                          
                      end                                            
                      
                      Ipot = obj.geometryVariables.Io.*exp( obj.variables.theta_frac(:,colIndex) );
                      infilt(:,colIndex ,:) = Ipot.*(1-exp(-P_infilt./Ipot));                        
                   end 

                   % Calculate saturation excess runoff
                   satRunoff(:, colIndex , :) = Peff(:, colIndex , :) - infilt(:, colIndex , :);  
                                      
                   % Calculate soil evaporation, evap.
                   evap(:, colIndex , :) = PET(:, i , :).*exp(-obj.geometryVariables.k_light .* obj.variables.LAI(:, colIndex , :) ) ... 
                       .* obj.variables.theta_frac(:, colIndex , :) .* relDemand_unsat(:, colIndex ,:);

                   % Calculate transpiration, transp.
                   transp(:, colIndex , :)  = PET(:, i , :).*(1-exp(-obj.geometryVariables.k_light .* obj.variables.LAI(:, colIndex , :) )) ...
                       .* obj.variables.theta_frac_veg(:, colIndex , :) ;   

                   % Groundwater evap.                   
                   evap_gw(:, colIndex , :) =  ( PET(:, i , :) - transp(:, colIndex , :) - evap(:, colIndex , :) ).* exp(-obj.geometryVariables.k_light .* obj.variables.LAI(:, colIndex , :) ) ...
                       .* relDemand_gw(:, colIndex ,:);                        
               end
               
               % Now that climate dependent fluxes have been calculated at
               % each percentile, numerically integrate the fluxes over the
               % distribution. The integration method used is Simpson's
               % extended method. It requires an ODD number of equally
               % spaced points.
               obj.fluxes.P = trapz(obj.climateData.scalingPercentiles, P, 3);
               obj.fluxes.PET = trapz(obj.climateData.scalingPercentiles, PET, 3);
               obj.fluxes.intercept = trapz(obj.climateData.scalingPercentiles, 3);
               obj.fluxes.Peff = trapz(obj.climateData.scalingPercentiles, Peff, 3);
               obj.fluxes.infilt = trapz(obj.climateData.scalingPercentiles, infilt, 3);
               obj.fluxes.satRunoff = trapz(obj.climateData.scalingPercentiles, satRunoff, 3);
               obj.fluxes.evap = trapz(obj.climateData.scalingPercentiles, evap,  3);
               obj.fluxes.transp = trapz(obj.climateData.scalingPercentiles, transp, 3);
               obj.fluxes.evap_gw = trapz(obj.climateData.scalingPercentiles, evap_gw, 3);
               
           else
               obj.fluxes.P = climateData_timePoint(:,1)';
               obj.fluxes.PET = climateData_timePoint(:,2)';
               
               % calculate interception, 'Intercept'
               interceptPotential = obj.geometryVariables.F_canopy./10 .* obj.variables.LAI;           
               warning('OFF','MATLAB:divideByZero');
               obj.fluxes.intercept = interceptPotential.*(1 - exp(-obj.fluxes.P./interceptPotential));
               warning('ON','MATLAB:divideByZero');
               obj.fluxes.intercept(obj.variables.LAI<=0) = 0;

               % Calculate effective precipitation, 'Peff'
               obj.fluxes.Peff = obj.fluxes.P - obj.fluxes.intercept;
              
               % Smooth function to est max P for infiltration
               % NOTE: The smoothening fn below caused a lot of problems. If lambda is
               % > approx 1 the Pinfiltration can be greater than P.
               % Also because of limitation in floating point precision,
               % exp((Peff-SM_capacity)/lambda) can be >> abs(1e308)                           
               obj.fluxes.infilt = zeros(size(obj.fluxes.Peff));
               if any(any(obj.fluxes.Peff>0));
                  P_infilt = obj.fluxes.Peff - lambda_p*log( exp((obj.fluxes.Peff - SM_capacity)/lambda_p)+1);                    
                  P_infilt(  isinf(P_infilt) ) = SM_capacity( isinf(P_infilt) );
                  P_infilt(  P_infilt<0 ) = 0;
                  Ipot = obj.geometryVariables.Io.*exp( obj.variables.theta_frac );    
                  obj.fluxes.infilt = Ipot.*(1-exp(-P_infilt./Ipot));                        
               end  
                                             
               % Calculate saturation excess runoff
               obj.fluxes.satRunoff = obj.fluxes.Peff - obj.fluxes.infilt;
               
               % Calculate soil evaporation, evap.
               obj.fluxes.evap = obj.fluxes.PET.*exp(-obj.geometryVariables.k_light .* obj.variables.LAI) ... 
                   .* obj.variables.theta_frac .* relDemand_unsat;

               % Calculate transpiration, transp.
               obj.fluxes.transp  = obj.fluxes.PET.*(1-exp(-obj.geometryVariables.k_light .* obj.variables.LAI)) ...
                   .* obj.variables.theta_frac_veg;   
               
               % Groundwater evap.               
               obj.fluxes.evap_gw =  (obj.fluxes.PET  - obj.fluxes.transp - obj.fluxes.evap).* exp(-obj.geometryVariables.k_light .* obj.variables.LAI ) ...
                   .* relDemand_gw; 
           end 
           
           % Calculate throughflow, recharge
          obj.fluxes.recharge = obj.variables.k_unsat .* ( 1 + obj.variables.dpsidz );  
          obj.fluxes.recharge(isinf(obj.variables.dpsidz) | obj.variables.DBNS<=0 )=0;          
           
           % Saturated fluxes (excludes groundwater evap)
           %------------------------------   
           %  Set change in soil water due to head change to null. It is set within setDerivatives(t,u,varargin).
           obj.fluxes.headChange_masstoSoil = [];             
           
           % Baseflow
           lambda_baseflow=0.01;            
           obj.fluxes.Q_toStream = -365.25 .* obj.modelParameters.stream.river_conductScaler .* obj.geometryVariables.k ...
               .* lambda_baseflow .* log( 1 + exp( ( obj.variables.h_vert - obj.geometryVariables.streamElevation) ./lambda_baseflow));
           if any(any(isinf(obj.fluxes.Q_toStream)))
               filt = isinf(obj.fluxes.Q_toStream);
               Q_toStream_temp = -365.25 .* obj.modelParameters.stream.river_conductScaler .* obj.geometryVariables.k .* ...
               ( obj.variables.h_vert - obj.geometryVariables.streamElevation );
               obj.fluxes.Q_toStream(filt) = Q_toStream_temp(filt);
           end       
           % Lower catchment boundary saturated flux           
           Param1 = obj.modelParameters.boundaryConditions.Lower.Param1;
           Param2 = obj.modelParameters.boundaryConditions.Lower.Param2;
           switch  obj.modelParameters.boundaryConditions.Lower.Type
               case 'dhdx'
                   obj.fluxes.boundaryFlow_lower = obj.variables.k(1,:).*obj.stateVariables.S(1,:) ...
                       ./obj.variables.mu(1,:) .* Param1;                   
               case 'Q'
                   obj.fluxes.boundaryFlow_lower = Param1;
                   
               case 'GHB'
                   obj.fluxes.boundaryFlow_lower = 365.25.*Param1 .* obj.modelParameters.geometry.width(1,:) ...
                       .* (obj.variables.h_vert(1,:) - Param2);
                   
               case 'GHB_2'
                   % GHB like boundary but uses k and Bsat from x=0 and parameters
                   % for head at x<0 (ie x_boundary).
                   obj.fluxes.boundaryFlow_lower = obj.variables.k(1,:) .* obj.modelParameters.geometry.width(1,:)  ...
                       .* obj.variables.b(1,:) .* ( obj.variables.h_vert(1,:) - Param2)./Param1;
           end
           
           % Upper catchment boundary saturated flux           
           Param1 = obj.modelParameters.boundaryConditions.Upper.Param1;
           Param2 = obj.modelParameters.boundaryConditions.Upper.Param2;
           switch  obj.modelParameters.boundaryConditions.Upper.Type
               case 'dhdx'
                   obj.fluxes.boundaryFlow_upper = obj.variables.k(end,:).*obj.stateVariables.S(end,:) ...
                       ./obj.variables.mu(end,:) .* Param1;                  
               case 'Q'
                   obj.fluxes.boundaryFlow_upper = Param1;

               case 'GHB'
                   obj.fluxes.boundaryFlow_upper = 365.25.*Param1 .* obj.modelParameters.geometry.width(end,:) ...
                       .* (obj.variables.h_vert(end,:) - Param2);
               case 'GHB_2'
                   % GHB like boundary but uses k and Bsat from x=0 and parameters
                   % for head at x<0 (ie x_boundary). 
                   obj.fluxes.boundaryFlow_upper = obj.variables.k(end,:) .* obj.modelParameters.geometry.width(end,:)  ...                       
                       .* obj.variables.b(end,:) .* ( obj.variables.h_vert(end,:) - Param2)./Param1;                   
           end           
           %------------------------------      
       end       
       
       % Get the climate forcing rate at a user specified time.
       function climateData = getClimate(obj, t)
            % check t
            if size(t,1) ~= 1;
                error('The input time data must be of only one row (but can be of n columns)');
            elseif size(t,2)>1 && t(1)>=t(end);
               error('When extracting the climate date for PDE solver, t(end)  must be greater than t(1)');
            end

            switch obj.climateData.type
                case 'data'
                    % Subtract eps from t to ensure that when t is very close to the end of
                    % a climate time step that the the climate for the end of the step is
                    % retrieved and not the next months.
                    teps = t - sqrt(eps(t));

                    % Get climate data from matrix.                    
                    for i=1:size(t,2)
                        if obj.climateData.doClimateScaling
                            [irow,icol]= find( obj.climateData.data.P(:,1) > teps(1,i) , 1, 'first');                                                                                                    
                            climateData(i,:) = [obj.climateData.data.P(irow,2:end), obj.climateData.data.ET(irow,2:end)];
                        else
                            [irow,icol]= find( obj.climateData.data.P(:,1) > teps(1,i) , 1, 'first');                                                                                                    
                            climateData(i,:) = [obj.climateData.data.P(irow,2), obj.climateData.data.ET(irow,2)];
                        end
                    end        
                case 'function'                    
                    if obj.climateData.doClimateScaling
                        nprcentiles =  length(obj.climateData.scalingPercentiles);
                        climateData = zeros(length(t),nprcentiles*2);
                        for i=1:size(t,2)                        
                            climateData(i,1:nprcentiles) = ppval( obj.climateData.functionParameters_P, t(i)-floor(t(i)) );
                            climateData(i,1+nprcentiles:2*nprcentiles) = ppval( obj.climateData.functionParameters_ET, t(i)-floor(t(i)) );
                        end                        
                    else
                        climateData(:,1) = ppval( obj.climateData.functionParameters_P, t-floor(t) );
                        climateData(:,2) = ppval( obj.climateData.functionParameters_ET, t-floor(t) );
                    end
                    climateData = max(0, climateData);
                otherwise
                    error('The input climate data type is not known!');
            end                
       end
   end

end 
