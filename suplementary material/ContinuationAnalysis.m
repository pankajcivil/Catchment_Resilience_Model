classdef ContinuationAnalysis < CatchmentResilienceModel
%CONTINUATIONANALYSIS Class definition for continuation analysis of the 'CatchmentResilienceModel'
%
% Description
%   This class allows the number and state-space location of attractors
%   and repellors in the 'CatchmentResilienceModel' to be quantified with
%   a change in a single model parameter.
%
%   This procedure is called codim 0 continuation analysis, but is also 
%   known as bifurcation analysis. It can quantify the number and location
%   of attractors and repellors when they are a within-year seasonal cycle
%   (formally, a limit-cycle continuation) or a non-seasonal steady state
%   (formally, equilibrium continuation). The actually numerical
%   continuation analysis is undertaken using a version of MatCont adapted
%   for limit-cycle continuation adapted resulting from climate forcing,
%   and not internal equation dynamics (see Dhooge et al. 2003).
%   Importantly, the calculation of the Jacobian matrix for the
%   limit-cycle continuation is undertaken using 'C-MEX', not MatLab code.
%   To run the code on your PC you may need to compile this 'C-MEX code. To
%   compile it, change the current path within MatLab to 
%   '.../continuation/modified_matcont/LimitCycle' (where ... is the
%   location at which you have installed this software) and enter the
%   following command into the MatLab command window: mex *.mex
%
%   For an example of how to build and run a model see the documentation
%   for the file 'CatchmentResilienceModel.m'. 
%
%   One a 'CatchmentResilienceModel' has been built a time-integration
%   solution has been derived, the continuation analysis can be commenced
%   from the model parameters within the solved model. For an example of how
%   to undertake continuation analysis see the example below.
%
%   For more details of the model see Peterson et al. (2009). For examples
%   of using this model to explore the attractors under daily forcing see
%   Peterson and Western (2012). For concepts and methods in hydrological 
%   resilience see Peterson, Western and Argent (2012).
%   
% Example: 
%   % Create the continuation object. This input most be an object of type
%   % 'CatchmentResilienceModel' that has been successfully solved.
%   % Experience with continuation suggests that this model should be solved
%   % from a deep water table initial condition. In the example below,
%   % 'model' is an object of type 'CatchmentResilienceModel'.
%   model_cont = ContinuationAnalysis(model);
%
%   % Set the continuation solver options to the default setting.
%   setContinuationOptions(model_cont)
%
%   % Set initial conditions. Note, if the parameter for continuation analysis
%   % contains the terms 'climateData.data' then equilibrium continuation
%   % will be undertaken for the requested climate variable. Importantly, the
%   % other climate variable must be a constant value. For equilibrium
%   % continuation of a daily climate rate, the 'CatchmentResilienceModel'
%   % model must have the climate forcing set to option '1' with each day of
%   % simulation having an identical precipitation and evapotranspiration rate.
%   %
%   % For limit-cycle continuation of a model parameter, in this case
%   % K_sat, enter the following:
%   setContinuationInitialConditions(model_cont, 'aquifer.Units{1,1}.k', 12, 5)
%
%   % Alternatively, for equilibrium continuation of daily precipitation
%   % enter the following: 
%   setContinuationInitialConditions(model_cont, 'climateData.data.P')
%
%   % Alternatively, for equilibrium continuation of daily potential
%   % evapotranspiration enter the following:
%   setContinuationInitialConditions(model_cont, 'climateData.data.ET')
%
%   % If undertaking limit-cycle continuation then it is good practice to
%   % first check that the within-year cycles have converged to a stable cycle.
%   % To assess if the initial condition cycles are stable, they can be
%   % plotted as follows:
%   xnodes = [100; 250; 1000; 1500; 1750];
%   plotInitialSolutionConvergence(model_cont, 10, xnodes )
%
%   % If the initial conditions are stable, then continuation analysis can
%   % be undertaken.
%   doContinuation(model_cont)
%   
%   % If the continuation analysis successfully finished, the results can
%   % be plotted at user input model nodes.
%   xnodes = [100; 250; 1000; 1500; 1750];
%   plotContinuation(model_cont, xnodes)
%
% See also:
%   ContinuationAnalysis: model_construction;
%   setContinuationOptions: set_continuation_solver_settings;
%   setContinuationInitialConditions: set_initial_conditions;
%   plotInitialSolutionConvergence: check_initial_conditions;
%   doContinuation: do_the_continuation;
%   plotContinuation: plot_continuation_results;
% 
% Dependencies (modified from MatCont 2.3.3):
%   init_LC_LC.m
%   init_EP_EP.m
%   contset.m
%   plotcycle.m
%   cpl.m
%   cont.m
%   newtcorr.m
%   contidx.m
%   equilibrium.m
%   ejac.m
%   ejacp.m
%   limitcycle.m

% References:
%   Dhooge, A., W. Govaerts, and Y. A. Kuznetsov (2003), MATCONT: a MATLAB
%   package for numerical bifurcation analysis of ODEs, ACM transactions
%   on mathematical software, 29 (2), 141–164.
%
%   Peterson, T. J., R. M. Argent, A. W. Western, and F. H. S. Chiew
%   (2009), Multiple stable states in hydrological models: An ecohydrological
%   investigation, Water Resour. Res., 45, W03406, doi:10.1029/2008WR006886. 
%
%   Peterson, T. J. and A. W. Western (2012), On the 
%   existence and emergence of multiple hydrological attractors under 
%   stochastic daily forcing: 1. Identifying the attractors, Water 
%   Resour. Res., (Submitted August 2012)
%
%   Peterson, T. J., A. W. Western, and R. M. Argent (2012), Analytical 
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

       % Type of continuation analysis. True: do limit-cycle; False:do equilibrium continuation 
       % (this will be adopted for climate forcing continuations)
       continuationDoLimitCycle;
       
       % Continuation parameter name and initial value.       
       continuationParameterName; 
       continuationInitialConditions;
              
       % Continuation solve options
       continuationOptions; 
       
       % Time integration results
       continuationResults;
       
   end

   methods
%% Create instance of ContinuationAnalysis object
       function obj_CA = ContinuationAnalysis(model)
% CONTINUATIONANALYSIS Model construction.
%
% Syntax:
%   obj_CA = ContinuationAnalysis(model)
%
% Description:
%   Constructs an instance of the continuation model object. The input can
%   be an object of type 'CatchmentResilienceModel' or
%   'ContinuationAnalysis'. If the latter, then the settings for
%   continuation will be copied from the input. 
%
%   Once the model object is created, the following steps should be followed
%   to run the model:
%   to run the model:
%   1. Set the continuation initial conditions.
%   2. Set the continuation solver options;
%   3. Set the climate forcing;
%   4. Solve the model over a user defined time range;
%
% Inputs:
%   modelInputs - an object of type 'CatchmentResilienceModel' OR
%   'ContinuationAnalysis'.
%
% Outputs:
%   obj_CA - an output object of type 'ContinuationAnalysis'.
%
% See also:
%   ContinuationAnalysis: class_description;
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
%

          if ~isa(model,'ContinuationAnalysis') &&  ~isa(model,'CatchmentResilienceModel')
              error([' Incorrect input! The input must be either an instance of the CatchmentResilienceModel', char(13), ...
                     ' object (ie a model solution) or a ContinuationAnalysis object.']);
          elseif isfield(model.climateData, 'type') && ~strcmp(model.climateData.type,'function')
              display(['WARNING: If not undertaking equilibrium continuation for climate forcing, the', char(13), ...
                       'climate forcing must be a smooth repeating function and not discrete and or', char(13), ...
                       'stochastic data!']);
          elseif ~isfield(model.results, 'S') || ~isfield(model.results, 't')
              % Check time integration solution exists.
              error(' Time integration solutions must be derived prior to starting continuation analysis!');              
          end
           
          %Create instance of CatchmentResilienceModel class and assign
          % input model properties.
          obj_CA = obj_CA@CatchmentResilienceModel(model);     
           
          if isa(model,'ContinuationAnalysis')
               % Copy over properties from existing
               % ContinuationAnalysis object. 
                obj_CA.continuationParameterName = model.continuationParameterName;
                obj_CA.continuationInitialConditions = model.continuationInitialConditions;
                obj_CA.continuationOptions = model.continuationOptions;
                obj_CA.continuationResults = model.continuationResults;          
          end
          
          % Modify file paths to ensure the modified MatCont functions are
          % called.
          %----------------------------------------------------------------          
          % Get the path of the current function.
          functionPath = mfilename('fullpath');
          
          % Remove the function name from the file path
          if ispc
              filt = find(functionPath=='\', 1, 'last');
          elseif ismac || isunix
              filt = find(functionPath=='/', 1, 'last');
          else
              error('The continuation algorithm can only be ran on Linux, PC or Mac machines.');
          end
          functionPath = functionPath(1:filt);     
         
          % Add the file paths to the modified MatCont functions
          if ispc
            addpath([functionPath,'continuation\modified_matcont\Continuer']);
            addpath([functionPath,'continuation\modified_matcont\Equilibrium']);
            addpath([functionPath,'continuation\modified_matcont\LimitCycle']);
          elseif ismac || isunix
            addpath([functionPath,'continuation/modified_matcont/Continuer']);
            addpath([functionPath,'continuation/modified_matcont/Equilibrium']);
            addpath([functionPath,'continuation/modified_matcont/LimitCycle']);
          end
          %----------------------------------------------------------------
       end
       
%% Set continuation options.
       function setContinuationOptions(obj_CA, userContinuationOptions)
% SETCONTINUATIONOPTIONS Set the solver options for continuation analysis.
%
% Syntax:
%   setContinuationOptions(obj)   
%   setContinuationOptions(obj, userContinuationOptions)   
%
% Description:
%   The method sets the options for the continuation solver. A default 
%   structure can be set, or modifications can be made to this and
%   applied to the object. The method amends the MatCont solver settings 
%   (from contset.m) and adds additional settings for the parameter range
%   and the 'CatchmentResilienceModel' differential equation solver.
%   Specifically: (i) the differential equations are not solved as a 
%   differential algebraic equation; (ii) differential equation warnings
%   are suppressed; and (iii) the potential for a singularity emerging
%   within the differential equations is reduced by the setting a scaler
%   to 0.5 within the outer terms of the model differential equation for 
%   lateral groundwater flow. For more details of the
%   singularity see the documentation for 'solveModel' within the
%   documentation for 'CatchmentResilienceModel'.
%
% Inputs:
%   obj_CA - an object of type 'ContinuationAnalysis'.
%
%   userContinuationOptions - structure variable containing a full set of 
%   solver options modified to the user's requirements. To obtain an
%   example of the default structure variable, build a model object and
%   call this method without the input 'userContinuationOptions'. The 
%   default structure is given within 'obj_CA.continuationOptions' where 
%   'obj_CA' is the name of the model.
%
% Outputs:
%   (model object, obj_CA, updated to contain the solver options.)
%
% See also:
%   ContinuationAnalysis: class_description;
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

              obj_CA.solverOptions.useDAE = false;
              obj_CA.modelParameters.aquifer.singularityFraction = 0.5;
              obj_CA.solverOptions.showWarnings = false;

          if nargin==1;
              display(' .. Continuation options set to default values for this model.');
              % Get structure of user options.
              obj_CA.continuationOptions = contset;
              
              % Set options to those found to be important for continuation
              % of Ksat paramerter within the CatchmentResilienceModel.
              obj_CA.continuationOptions.InitStepsize = 250;
              obj_CA.continuationOptions.MinStepsize =  0.001;
              obj_CA.continuationOptions.MaxStepsize = 1000;              
              obj_CA.continuationOptions.MaxCorrIters = 12;                            
              obj_CA.continuationOptions.MaxNewtonIters = 4;              
              obj_CA.continuationOptions.Increment = 0 + 1e-6i;              
              obj_CA.continuationOptions.FunTolerance = 1e-4;
              obj_CA.continuationOptions.VarTolerance = 1e-6;              
              
              obj_CA.continuationOptions.MaxNumPoints = 10000;
              obj_CA.continuationOptions.Backward = false;
              obj_CA.continuationOptions.ParamMin = 0.05;              
              obj_CA.continuationOptions.ParamMax = 16;
              obj_CA.continuationOptions.odeset = obj_CA.solverOptions;
              obj_CA.continuationOptions.MaxSecantIters = 1;
              obj_CA.continuationOptions.MoorePenrose = true;
                            
              obj_CA.continuationOptions.SymDerivative = false;
              obj_CA.continuationOptions.SymDerivativeP = false;
              
              obj_CA.solverOptions.useDAE = false;
              obj_CA.modelParameters.aquifer.singularityFraction = 0.5;
              obj_CA.solverOptions.showWarnings = false;
          else
              obj_CA.continuationOptions = userContinuationOptions;
          end
       end
       
%% Set limit-cycle continuation initial conditions
       function setContinuationInitialConditions(obj_CA,  ParameterName, ntst_user, ncol_user )  
% SETCONTINUATIONINITIALCONDITIONS Set limit-cycle continuation initial conditions.
%
% Syntax:
%   setContinuationInitialConditions(obj_CA,  ParameterName )  
%   setContinuationInitialConditions(obj_CA,  ParameterName, ntst_user, ncol_user )  
%
% Description:
%   The method sets the model parameter that will be varied during the
%   continuation analysis and additional settings required for MatCont to
%   undertake the continuation. If required, time-integrations simulations
%   are also undertaken to produce the required state-variable initial 
%   conditions. This is most relevant to limit-cycle continuation as it
%   requires ntst_user * ncol_user time points per year. If equilibrium
%   continuation is to be undertaken then constant daily precipitation and
%   potential evapotranspiration rates are applied and time-integration
%   undertaken. 
%
%   For all time-integration simulations, the simulation duration is 50 
%   years.  For complex problems this may be insufficient. If so, then the
%   model 'CatchmentResilienceModel' should be re-run for a longer duration
%   and assessed, using time-series plots, to have converged to a stable
%   solution.
%
%   Finally, numerical limit-cycle continuation of the nonlinear partial
%   differential equation model from 'CatchmentResilienceModel' is a
%   numerically complex undertaking and very dependent upon the model
%   parameters, climate forcing and continuation settings. The stability of
%   the continuation can be sensitive to the initial depth to water table 
%   and parameter values. The authors have found that it is often best to
%   start the continuation from a deep initial condition and at a  parameter
%   value having only the deep water tabe attractor. However, if a stable
%   limit-cycle continuation solution cannot be achieved it may be
%   nessesary to increase the number of time points per year. This can be
%   achieved by increasing the inputs 'ntst_user' and 'ncol_user'.
%
% Inputs:
%   obj_CA - an object of type 'ContinuationAnalysis'.
%
%   ParameterName - string for the parameter to be investigated in the
%   continuation analysis. The input can be a parameter defined within the
%   'CatchmentResilienceModel' object or one of the following climate
%   inputs: 'climateData.data.P' or 'climateData.data.ET'. If the latter,
%   then equilibrium continuation is undertaken.
%
%   ntst_user - scalar integer greater than zero for limit-cycle
%   continuation controlling the number of test intervals per annum (see
%   MatCan manual for detail). Default value is 9.
%
%   ncol_user - scalar integer greater than zero for limit-cycle
%   continuation controlling the number of collocation points per annum (see
%   MatCan manual for detail). Default value is 4.
%
% Outputs:
%   (model object, obj_CA, updated to contain the initial condition settings.)
%
% See also:
%   ContinuationAnalysis: class_description;
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

            global cds lds eds

            if strcmp(ParameterName,'climateData.data.P') || strcmp(ParameterName,'climateData.data.ET')
                % Doing equilb. continuation for daily climate data.
                obj_CA.continuationDoLimitCycle = false;
                                
                t = obj_CA.results.t(end);
            else
                % Doing limit cycle continuation for a model parameter 
                obj_CA.continuationDoLimitCycle = true;
                
                
                % Set number of limit-cycle co-location and node points.
                %---------------
                if nargin==2;
                    display('.. Number of collocation points, ntst and ncol, per year set to 9 and 4 respectively.');
                    ntst = 9;
                    ncol = 4;
                else                
                    ntst = ntst_user;
                    ncol = ncol_user;                
                end
                
                % Find number of solution points in the final year of
                % simulation and assume last point is an end of year.
                %---------------
                t = obj_CA.results.t;
                nt = size(t,2);
                tend = t(end);

                % Check tend is an end of year point
                if abs(mod(tend,1)) >sqrt(eps(1));
                    error('The time integration solution must finish at the end of a year, ie the final time point must be an integer.');
                end

            end
          
            
            % Set initial parameter name and value.
            %---------------
            doTimeIntegration=false;
            if isempty(obj_CA.continuationParameterName) ...
            || ~strcmp( obj_CA.continuationParameterName, ParameterName) ...
            || (strcmp( obj_CA.continuationParameterName, ParameterName) ...
            && isfield(obj_CA.continuationInitialConditions) ... 
            && isempty(obj_CA.continuationInitialConditions.parameterValue))
        
                if isempty(obj_CA.continuationParameterName) ...
                || ~strcmp( obj_CA.continuationParameterName, ParameterName) 
                    obj_CA.continuationParameterName = ParameterName;
                end
                
                if obj_CA.continuationDoLimitCycle
                    obj_CA.continuationInitialConditions.parameterValue = getParameterValue(obj_CA, ParameterName);
                else
                    climateData = getClimate(obj_CA, t);
                    if strcmp(ParameterName(end), 'P')
                        obj_CA.continuationInitialConditions.parameterValue = climateData(1);
                    elseif strcmp(ParameterName(end), 'ET')
                        obj_CA.continuationInitialConditions.parameterValue = climateData(2);
                    else
                        error('The continuation analysis for the climate data was expecting an input of the form "climateData.data.P" or "climateData.data.ET".')
                    end
                end
                doTimeIntegration = true;
                                
            elseif ~isfield(obj_CA.continuationInitialConditions, 'matcont') ...
            || ~obj_CA.continuationDoLimitCycle ...
            || obj_CA.continuationInitialConditions.matcont.s0.data.ntst ~= ntst ...
            || obj_CA.continuationInitialConditions.matcont.s0.data.ncol ~= ncol        
        
                doTimeIntegration = true;            
            end    
            %---------------
            
            % Assign initial parameter value back as initial condition
            nyears = 50;
            if obj_CA.continuationDoLimitCycle
                setParameterValue(obj_CA, obj_CA.continuationParameterName, obj_CA.continuationInitialConditions.parameterValue) 
            else
                % Create climate series of a constant rate equal to the
                % continuation value.
                time_points = (1:(nyears+1)*365);
                time_points = datevec(time_points);
                time_points = time_points(:,1:3);
                nt = size(time_points,1);
                climateData = getClimate(obj_CA, t);
                setClimateData(obj_CA, 1, true, [ time_points, repmat(climateData, nt,1)],[]);
            end
                        
            % If required, conduct time-integration simulaton.
            %---------------
            if doTimeIntegration
               display('.. The time-integration simulation does not contain the required number of collocation points.');
               display('   A short time-integration simulation will be performed at the required number of points.');
               
               % Extract state variables at end of the user-supplied
               % simulation. 
               if isfield(obj_CA.continuationInitialConditions, 'S0')
                   S0 = obj_CA.continuationInitialConditions.S0;
                   m0 = obj_CA.continuationInitialConditions.m0;
                   t0 = obj_CA.continuationInitialConditions.t0;
               else
                   S0 = real(obj_CA.results.S(:,end));
                   m0 = real(obj_CA.results.m(:,end));               
                   t0 = real(obj_CA.results.t(end));               
               end
            
               % Input state variables as new initial condition.
               setInitialConditions(obj_CA, [S0, m0], 'node_storage');               
               
               % Define time points required for simulation.
               tstart = 1;
               tend = nyears;
               if obj_CA.continuationDoLimitCycle
                    nTimePoints = ntst * ncol;
               else
                    nTimePoints=1;
               end
               
               % Conduct simulations.
               solveModel(obj_CA, tstart, tend, nTimePoints );
                
               % Re-assign prior set time step variables
               t = obj_CA.results.t;
               if obj_CA.continuationDoLimitCycle
                    ntspan = ntst * ncol + 1;    
               else
                   ntspan = 1;
               end
            
               % Extract state variables and initial parameter value.
               %---------------
               u0 = obj_CA.results.u;
               y0=zeros(2*size(u0,1),ntspan);
               y0(1:2:end-1,:) = real(obj_CA.results.S(:,end-ntspan+1:end));
               y0(2:2:end,:) = real(obj_CA.results.m(:,end-ntspan+1:end));
               ParamValue = obj_CA.continuationInitialConditions.parameterValue;
               
               display('.. Getting initial equilibrium point');                              
               % If doing limit-cycle cont., set number of simulation time 
               % points per year and expand to finer mesh.
               %---------------                        
               if obj_CA.continuationDoLimitCycle
                   tfinemesh = t(size(t,2) - ntspan+1:end);        
                   tfinemesh = tfinemesh - tfinemesh(1);
                   tmesh = tfinemesh(1:(ntspan-1)/ntst:ntspan);
                   
                   s0.index=1;
                   s0.data.T=1;
                   s0.data.ntst= ntst;
                   s0.data.ncol=ncol;            
                   s0.data.timemesh= tmesh - tmesh(1);
                   s0.data.parametervalues = [1; ParamValue];
                   y0 = [reshape(y0,1,[])'; s0.data.parametervalues];               
               else
                   y0 = reshape(y0,1,[])';
               end
               
            else
               y0 = obj_CA.continuationInitialConditions.matcont.y0;
               s0 = obj_CA.continuationInitialConditions.matcont.s0;
               ParamValue = obj_CA.continuationInitialConditions.parameterValue;
            end

            % Call MatCont initial condition setup function
            %---------------
            if obj_CA.continuationDoLimitCycle
                [y0, v0] = init_LC_LC(@obj_CA.getDerivatives_Continuation, y0, [], s0, 2, s0.data.ntst, s0.data.ncol);  
            else
                [y0, v0] = init_EP_EP(@obj_CA.getDerivatives_Continuation, y0, ParamValue, 1);  
                cds.options = contset(cds.options, 'SymDerivative', 0);
                cds.options = contset(cds.options, 'SymDerivativeP', 0);
            end
            
            % Add odeset to cds.options
            cds.options.odeset = odeset;
            %---------------
            
            % Assign initial inputs to objectfile
            obj_CA.continuationInitialConditions.parameterValue = ParamValue;
            if ~isfield(obj_CA.continuationInitialConditions, 'S0')
                obj_CA.continuationInitialConditions.S0 = S0;
                obj_CA.continuationInitialConditions.m0 = m0;
                obj_CA.continuationInitialConditions.t0 = t0;
            end
            obj_CA.continuationInitialConditions.matcont.y0 = y0;
            obj_CA.continuationInitialConditions.matcont.v0 = v0;
            
            if obj_CA.continuationDoLimitCycle
                obj_CA.continuationInitialConditions.matcont.s0 = s0;
                obj_CA.continuationInitialConditions.matcont.lds=lds;
            else
                obj_CA.continuationInitialConditions.matcont.eds=eds;                
            end
            obj_CA.continuationInitialConditions.matcont.cds=cds;
            
            clear global;
       end
       
%% Do limit-cycle continuation 
       function doContinuation(obj_CA, priorLCC_solutionIndxes )
% DOCONTINUATION Undertake the actual numerical continuation analysis.
%
% Syntax:
%   doContinuation(obj_CA )
%   doContinuation(obj_CA, priorLCC_solutionIndxes )
%
% Description:
%   This method undertakens the actual codim 1 (i.e. limit-cycle or 
%   equilibrium) continuation analysis. The method first sets up the model
%   object for continuation analysis and then calls the appropriate MatCont
%   function to actually undertake the continuation. 
%
%   As discussed within the documentation for the method
%   'setContinuationInitialConditions', numerical continuation can be a
%   very complex numerical undertaking. The setting defined throughout this
%   code was found to work for Peterson and Western (2012). However, a
%   chnage of parameter or climate forcing may result in problems during
%   continuation analysis. If problems are encountered, please read the
%   documentation for 'setContinuationInitialConditions', Peterson and 
%   Western (2012), Dhooge et al. (2003) and the MatCont manual. 
% 
% Inputs:
%   obj_CA - an object of type 'ContinuationAnalysis' with the continuation
%   initial conditions and solver options set.
%
%   priorLCC_solutionIndxes - a scalar integer for an EXPERIMENTAL feature
%   allowing the continuation to start from a prior calculated continuation
%   step number. That is, if a previous continuation analysis took 100
%   parameter steps, then this feature can be used to re-start the
%   continuation from any one of these 100 prior solution steps.
%
% Outputs:
%   (model object, obj_CA, updated to contain the continuation results.)
%
%
% See also:
%   ContinuationAnalysis: class_description;
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
%
% References: 
%   Dhooge, A., W. Govaerts, and Y. A. Kuznetsov (2003), MATCONT: a MATLAB
%   package for numerical bifurcation analysis of ODEs, ACM transactions
%   on mathematical software, 29 (2), 141–164.
%
%   Peterson, T. J. and A. W. Western (2012), On the 
%   existence and emergence of multiple hydrological attractors under 
%   stochastic daily forcing: 1. Identifying the attractors, Water 
%   Resour. Res., (Submitted August 2012)
%
           
            % Check the the model is setup. 
            if isempty(obj_CA.continuationInitialConditions);
                error(['The initial conditions for continuation analysis must first be defined.',char(13),'Below is an example for setting the initial conditions:',char(13),'setContinuationInitialConditions(model_LCC,  "aquifer.Units{1,1}.k", 12, 4);']);
            elseif isempty(obj_CA.continuationOptions);
                error(['The solver options for continuation analysis must first be defined.',char(13),'Below is an example for setting the option to the default:',char(13),'setContinuationOptions(model_LCC);']);
            end
           
            % if priorLCC_solutionIndxes is input, then check it is <= the
            % existing prior number of solution nodes.
            if nargin>1;
                if priorLCC_solutionIndxes <1
                    error('The input "priorLCC_solutionIndxes" allows continuation from a prior solution. The number must be >=1!');
                end

                if isfield(obj_CA.continuationResults,'y')
                    if priorLCC_solutionIndxes > size(obj_CA.continuationResults.y,2)
                        error('The input "priorLCC_solutionIndxes" allows continuation from a prior solution. The number must be <= the number of prior solution nodes!');
                    end
                else
                    error('The input "priorLCC_solutionIndxes" allows continuation from a prior solution. No prior solutions yet exist!');
                end 
            else
                priorLCC_solutionIndxes=[];
            end
            
            global cds eds lds objGlobal_LC;
            objGlobal_LC = obj_CA;

%            % TEMP: equilibcont  dev
%            obj_CA.continuationOptions.Increment = 1e-6;   
%            obj_CA.continuationOptions.Backward = true;
%            obj_CA.continuationOptions.Adapt = 0;
%            obj_CA.continuationOptions.ParamMin = 0.05;
%            obj_CA.continuationOptions.ParamMax = 15;
%            obj_CA.continuationOptions.MaxStepsize = 250;
%            obj_CA.continuationOptions.FunTolerance = 1e-7;
%            obj_CA.continuationOptions.VarTolerance = 1e-8;
%            obj_CA.continuationOptions.MaxNumPoints = 1000;
%            obj_CA.continuationOptions.MaxCorrIters=24;
           % Check delta is sufficiently large to avoid rounding errors
           %--------------------------------
           % Round initial user estimate           
           del = obj_CA.continuationOptions.Increment;
           y_typical = 1000;
           if isreal(del);
               temp = y_typical + del;
               del = temp - y_typical;      
           else     
               temp = y_typical*1i + del;
               del = temp - y_typical*1i;
           end
           
           % Reproduce reductions in  delta resulting from Ridders method.
           ntab = 10;
           con = 1.4; 
           for ii=1:ntab
                del = del/con;                        
           end
           
           % Check if rounding occurs with final value of del.
           del_final = del;
           if isreal(del_final);               
               temp = y_typical + del_final;
               del_final = temp - y_typical;      
           else     
               temp = y_typical*1i + del_final;
               del_final = temp - y_typical*1i;
           end
           if abs(del - del_final) ~=0
               display([char(13), 'WARNING: The user set increment size for finite difference gradient calculation may produce rounding errors!', ...
                   char(13), 'The error in delta after Ridders method is approx: ', num2str(abs(del - del_final)), ...
                   char(13), 'Consider increasing the size of continuationOptions.Increment.', char(13)]);
           end
           %--------------------------------
            
            % Set limit-cycle continuation algorithm to initial conditions.
            if obj_CA.continuationDoLimitCycle
                if isempty(priorLCC_solutionIndxes)

                   init_LC_LC(@obj_CA.getDerivatives_Continuation, obj_CA.continuationInitialConditions.matcont.y0, ...
                         [], obj_CA.continuationInitialConditions.matcont.s0, ...
                         2, obj_CA.continuationInitialConditions.matcont.s0.data.ntst, ...
                         obj_CA.continuationInitialConditions.matcont.s0.data.ncol);                 
                else
                   s0 = obj_CA.continuationResults.s(1,1);
                   s0.index=1;
                   s0.data.T=1;
                   s0.data.parametervalues = [1; obj_CA.continuationResults.y(end,priorLCC_solutionIndxes)];
                   s0.data = rmfield(s0.data,'phi');

                   init_LC_LC(@obj_CA.getDerivatives_Continuation, obj_CA.continuationResults.y(:,priorLCC_solutionIndxes) , ...
                         []  ,s0 , ...
                         2, obj_CA.continuationInitialConditions.matcont.s0.data.ntst, ...
                         obj_CA.continuationInitialConditions.matcont.s0.data.ncol);                      
     
                end
                lds = obj_CA.continuationInitialConditions.matcont.lds;                
            else
                if isempty(priorLCC_solutionIndxes)
                    init_EP_EP(@obj_CA.getDerivatives_Continuation, obj_CA.continuationInitialConditions.matcont.y0, ...
                        obj_CA.continuationInitialConditions.parameterValue, 1);  
                else
                   init_EP_EP(@obj_CA.getDerivatives_Continuation, obj_CA.continuationResults.y(:,priorLCC_solutionIndxes) , ...
                         obj_CA.continuationResults.y(end,priorLCC_solutionIndxes), 1);                                          
                end
            end
            cds = obj_CA.continuationInitialConditions.matcont.cds;
            
            % do continuation analysis.
            display('.. Starting continuation');      
            if obj_CA.continuationDoLimitCycle
                if isempty(priorLCC_solutionIndxes)
                    [obj_CA.continuationResults.y, obj_CA.continuationResults.v, obj_CA.continuationResults.s, obj_CA.continuationResults.h] ...
                            = cont(@limitcycle, obj_CA.continuationInitialConditions.matcont.y0, obj_CA.continuationInitialConditions.matcont.v0, obj_CA.continuationOptions);            
                else
                    [obj_CA.continuationResults.y, obj_CA.continuationResults.v, obj_CA.continuationResults.s, obj_CA.continuationResults.h] ...
                            = cont(@limitcycle, obj_CA.continuationResults.y(:,priorLCC_solutionIndxes), [], obj_CA.continuationOptions);                            
                end
            else
                if isempty(priorLCC_solutionIndxes)
                    obj_CA.continuationInitialConditions.matcont.v0 = [];
                    [obj_CA.continuationResults.y, obj_CA.continuationResults.v, obj_CA.continuationResults.s, obj_CA.continuationResults.h] ...
                            = cont(@equilibrium, obj_CA.continuationInitialConditions.matcont.y0, obj_CA.continuationInitialConditions.matcont.v0, obj_CA.continuationOptions);            
                else
                    [obj_CA.continuationResults.y, obj_CA.continuationResults.v, obj_CA.continuationResults.s, obj_CA.continuationResults.h] ...
                            = cont(@equilibrium, obj_CA.continuationResults.y(:,priorLCC_solutionIndxes), [], obj_CA.continuationOptions);                            
                end
                
            end
            % input to LC options to obj_CA
            %obj_CA.continuationInitialConditions.matcont.lds=lds;
            %obj_CA.continuationInitialConditions.matcont.cds=cds;
            
            clear global;
       end      
       
%% Create phase space plot to check time-integreation solution is stable.
       function plotInitialSolutionConvergence(obj_CA, nyears, node_locations )
% PLOTINITIALSOLUTIONCONVERGENCE Plot phase space cycles to check initial conditions.
%
% Syntax:
%   plotInitialSolutionConvergence(obj_CA, nyears, node_locations )
%
% Description:
%   This method uses the results from the method
%   'setContinuationInitialConditions' to check that a stable initial
%   solution has been derived prior to undertaking the actual continuation
%   analysis. The method plots the depth to water table and soil moisture 
%   over a user specified number of years of time-integration solution at 
%   user sepcified model spatial nodes. If the time-integration simulation
%   has converged to a stable solution (i.e. an attractor) at all spatial
%   locations then these plots should not show an shifting of the cycles
%   (or points if undertaking equilibrium continuation) over time. 
%
% Inputs:
%   obj_CA - an object of type 'ContinuationAnalysis'.
%
%   nyears - scalar integer greater than zero defining the number of the
%   end years of time-integration to plot.
%
%   node_locations - vector of the model spatial nodes at which plots are 
%   required. This is input as a distance, not model node number. 
%
% Outputs:
%   (none.)
%
% See also:
%   ContinuationAnalysis: class_description;
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

            % Check if model solutions exist.            
            if isfield (obj_CA.results,'S');
                if isempty(obj_CA.results.S);
                    error('Time-integration results must first be derived before plotting the phase cycles!');
                end
            else
                error('Time-integration results must first be derived before plotting the phase cycles!');
            end                        
           
            % Get state variables.            
            u = obj_CA.results.u;
            t = obj_CA.results.t;              
            t_end = t(end);
            t_start = max(t(1), t_end - nyears);
            t_filter = t>=t_start & t<=t_end;
            t_filter_cols = find(t_filter);
            
            % Get calculated depth to water table, DBNS, and soil moisture, theta.
            for i=1:length(t_filter_cols);
                obj_CA.stateVariables.S = obj_CA.results.S(:, t_filter_cols(i));
                obj_CA.stateVariables.m = obj_CA.results.m(:, t_filter_cols(i));
                obj_CA.stateVariables.t = obj_CA.results.t(t_filter_cols(i));
                
                % Calculate depth to water table by calling setVariables().
                setVariables(obj_CA);
                
                DTWT(:,i) = obj_CA.variables.DBNS;
                theta(:,i) = obj_CA.variables.theta;                                                
            end
            
            % Extract Limit-Cycle initial phase cycle
            plotLCC_IC=false;
            if isfield(obj_CA.continuationInitialConditions,'matcont')
                if isfield(obj_CA.continuationInitialConditions.matcont,'y0')
                    plotLCC_IC = true;
                    
                    %reshape state variable vector, y0, to matrix
                    y0 = reshape( obj_CA.continuationInitialConditions.matcont.y0(1:end-2,1), 2*length(u),[]);
                    
                    % Get time points for phase cycle
                    t0 = obj_CA.continuationInitialConditions.matcont.lds.finemsh;
                    
                    for i=1: length(t0)
                        % Set state varaibles
                        obj_CA.stateVariables.S = y0(1:2:end,i);
                        obj_CA.stateVariables.m = y0(2:2:end,i);
                        obj_CA.stateVariables.t = t0(i);

                        % Calculate depth to water table by calling setVariables().
                        setVariables(obj_CA);

                        DTWT_LCC_IC(:,i) = obj_CA.variables.DBNS;
                        theta_LCC_IC(:,i) = obj_CA.variables.theta;                                                                        
                    end
                    
                end
            end
            
            % extract solutions at each x location
            nplots = size(node_locations,1);
            ncols=2;
            nrows = ceil(nplots/2);
            figure();        
            for j=1:nplots
                % Create cell node filter to extract location.
                u_filter = u==node_locations(j,1);
                
                % Filter DBNS data
                DTWT_u = DTWT( u_filter, :);
                theta_u = theta( u_filter, :);                                

                % plot data
                subplot( nrows,ncols,j); hold on;
                plot(DTWT_u, theta_u,'-k');
                scatter(DTWT_u(1,1), theta_u(1,1),'or');                
                xlabel('Depth to water table (m)');
                ylabel('Soil moisture fraction');
                
                % PLot limit cycle initial phase
                if plotLCC_IC
                    plot(DTWT_LCC_IC( u_filter, :), theta_LCC_IC( u_filter, :),'-b');
                end
                
                % Add legend
                if plotLCC_IC && j==1
                    legend('Time-integration', 'Time-integration final', 'Limit-cycle initial phase cycle');
                elseif j==1;
                    legend('Time-integration', 'Time-integration final');
                end
                
                title([' Limit cycle at: ', num2str(node_locations(j,1))]);
                hold off;            
            end                         
       end
       
%% Create 3-D plot of continuation results
       function [DTWT_theta, u, parameterValues] = plotContinuation(obj_CA, nodesToPlot, ParamRangeToPlot)
% PLOTCONTINUATION 3-D plots of continuation results.
%
% Syntax:
%   [DTWT_theta, u, parameterValues] = plotContinuation(obj_CA, nodesToPlot)
%   [DTWT_theta, u, parameterValues] = plotContinuation(obj_CA, nodesToPlot, ParamRangeToPlot)
%
% Description:
%   This method creates a 3-D plot of the continuation results at user
%   specified model spatial locations. The plots are produced as the
%   parameter value investigated versus depth to water table versus soil
%   moisture. If equilibrium continuation of the climate forcing rate was
%   undertaken then the plots will be produced as the climatic rate instead
%   of a parameter value. If limit-cycle continuation was undertaken, then
%   at each parameter value the depth to water table and soil moisture
%   should form a closed elliptical-like structure. 
%
%   To assess if multiple attractors exist for a given parameter value (or
%   forcing rate) use the MatLab figure menu to rotate the figure. If
%   multiple attractors exist then, over a parameter range, three values of
%   depth to water table should exist. The middle depth is most likely to
%   be the repellor (i.e the threshold) between the attractors.
%
%   Importantly, the emergence of multiple attractors (i.e. which ones 
%   stochastic forcing can switch the system too) can significantly differ 
%   from those found to exist using continuation analysis. Please read
%   Peterson et al. (2012b, 2012,b).
%
% Inputs:
%   obj_CA - an object of type 'ContinuationAnalysis'.
%
%   nodesToPlot - vector of the model spatial nodes at which plots are 
%   required. This is input as a distance, not model node number. 
%
%   ParamRangeToPlot - 1x2 vector of parameter range to plot. If not input
%   then the entire range from the continuation will be plotted.
%
% Outputs:
%   DTWT_theta - 2nx x n_param_vals matrix (where nx is the number of model
%   nodes) of the depth to water table and soil moisture at each parameter
%   value from the continuation. The rows of the output alternate from
%   depth to water table in row 1 to soil moisture in row 2, repeating nx 
%   times. All model nodes are output.
%
%   u - nx x 1 vector of the spatial distance at each model node.
%
%   parameterValues - 1 x n_param_vals vector of the parameter values from
%   the continuation.
%
% See also:
%   ContinuationAnalysis: class_description;
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
%
% References:
%   Peterson, T. J., A. W. Western and Argent R. M. (2012a), On the 
%   existence and emergence of multiple hydrological attractors under 
%   stochastic daily forcing: 2. Can multiple attractors emerge?, Water 
%   Resour. Res., (Submitted August 2012)
%
%   Peterson, T. J., A. W. Western, and R. M. Argent (2012b), Analytical 
%   methods for ecosystem resilience: A hydrological investigation, Water 
%   Resour. Res., (Accepted August 2012)
%
            global lds cds
            
            if nargin<2
               error('This method required at least two inputs of limit-cycle object and the spatial nodes to plot.'); 
            elseif nargin<3
               ParamRangeToPlot = [-inf, inf];
            end
            
            % Set various constants
            if obj_CA.continuationDoLimitCycle
                lds = obj_CA.continuationInitialConditions.matcont.lds;
                ntst = obj_CA.continuationResults.s(1,1).data.ntst;
                ncol = obj_CA.continuationResults.s(1,1).data.ncol;
            else
                lds = [];
                ntst = 1;
                ncol = 1;
            end
            cds = obj_CA.continuationInitialConditions.matcont.cds;
            
            nu = size(obj_CA.results.u,1);
            npts = size(obj_CA.continuationResults.y,2);             
            nLC_rows = size(obj_CA.continuationResults.y,1);
           
            % Extract state variables for each limit cycle and calculate depth to water table (DTWT) and theta.            
            for k=1:npts
                for i=1:ntst;
                    for j=1:ncol;
                        range0 = 2*nu*( (i-1)*ncol +j-1)+1 : 2 : 2*nu*( (i-1)*ncol +j)-1; 
                        
                        % Set model state variables
                        obj_CA.stateVariables.S = obj_CA.continuationResults.y(range0,k);
                        obj_CA.stateVariables.m = obj_CA.continuationResults.y(range0+1,k);
                        obj_CA.stateVariables.t = 1;

                        % Calculate depth to water table (DTWT) and theta
                        setVariables(obj_CA);
                         
                        % Extract depth to water table and theta.
                        DTWT_theta(range0,k) = obj_CA.variables.DBNS;
                        DTWT_theta(range0+1,k) = obj_CA.variables.theta;                                                 
                    end
                end  
                
                % LC boundary condition points
                if obj_CA.continuationDoLimitCycle
                    range0 = nLC_rows-2*nu-1: 2 : nLC_rows-3; 
                    obj_CA.stateVariables.S = obj_CA.continuationResults.y(range0,k);
                    obj_CA.stateVariables.m = obj_CA.continuationResults.y(range0+1,k);
                    obj_CA.stateVariables.t = 1;                    
                    setVariables(obj_CA);
                    
                    DTWT_theta(range0,k) = obj_CA.variables.DBNS;
                    DTWT_theta(range0+1,k) = obj_CA.variables.theta;
                end
                            
            end
           
            % Plot limit cycle
            nnodesToPlot = length(nodesToPlot);
            for i=1: nu
                u( 2*(i-1)+1,1) = obj_CA.results.u(i);
                u( 2*(i-1)+2,1) = obj_CA.results.u(i);
            end
            if obj_CA.continuationDoLimitCycle
                u = repmat( u, ntst*ncol+1, 1);
            end
            
            for i=1:nnodesToPlot
                % Extract data for limit cycle at node location
                range0 = find(u == nodesToPlot(i)); 
                
                % Check if a solution node exists at this point.
                if sum(range0 )==0
                    error(['No solution nodes were found at the following user input location: ', num2str(nodesToPlot(i)), ' metres']);
                end
                
                DTWT_theta_atU = DTWT_theta(range0, :);
                
                if obj_CA.continuationDoLimitCycle
                    DTWT_theta_atU(2*(ntst*ncol+1)+1,:) = obj_CA.continuationResults.y( end-1,:);
                    DTWT_theta_atU(2*(ntst*ncol+1)+2,:) = obj_CA.continuationResults.y( end,:);                     
                else
                    DTWT_theta_atU(3,:) = obj_CA.continuationResults.y( end,:); 
                end
                
                v_atU = obj_CA.continuationResults.v(range0,:);
                if obj_CA.continuationDoLimitCycle
                    v_atU(2*(ntst*ncol+1)+1,:) = obj_CA.continuationResults.v( end,:);   
                    lds.nphase=2;
                else
                    v_atU(3,:) = obj_CA.continuationResults.v( end,:);   
                end
                
                % Filter data for parameters within range
                range0 =   obj_CA.continuationResults.y( end,:) >= ParamRangeToPlot(1) ...
                         & obj_CA.continuationResults.y( end,:) <= ParamRangeToPlot(2);
                DTWT_theta_atU = DTWT_theta_atU(:,range0);
                v_atU = v_atU(:,range0);
               
                %Plot data using MatCont routine 'plotcycle'
                figure();    
                hold on;                           
                if obj_CA.continuationDoLimitCycle
                    plotcycle(DTWT_theta_atU ,v_atU , obj_CA.continuationResults.s, [size(DTWT_theta_atU,1), 2, 1]);
                else
                    cpl([DTWT_theta_atU(3,:); DTWT_theta_atU(2,:); DTWT_theta_atU(1,:);],v_atU , obj_CA.continuationResults.s)
                end
                title(['u = ',num2str(nodesToPlot(i)), ' metres from outlet']);
                xlabel(['Parameter: ', obj_CA.continuationParameterName]);
                zlabel('Depth to watertable (m)');
                ylabel('Soil moisture ({\theta})');
                set(gca,'ZDir','reverse');      
                xlim(ParamRangeToPlot);
                box on;
 
                %rotate figure
                view(-32.0,30);

                hold off; 
                
            end
            
            parameterValues = DTWT_theta_atU(end,:);
            clear global;
       end

%% The following methods are for internal use only.  
%  Very minimal documentation is provided.
%--------------------------------------------------------------------------       

% Get the threshold at a user input parameter value.
       function [DTWT_thresholds, SM_thresholds]= getThreholds(obj_CA, nodesToExtract, param_target)
% GETTHREHOLDSEXPERIMENTRAL Experimental function to get the threshold at a user input parameter value.          

            % Get indexes for first parameter value greater than and less
            % than 'param_val'. This is achieved by tracing along limit
            % cycle parameter points to find the fold points.
            params_LCC = obj_CA.continuationResults.y( end,:);
            if params_LCC(2) < params_LCC(1)
                dir =-1;
            elseif params_LCC(2) > params_LCC(1)    
                dir=1;
            else
                error('The direction of continuation analysis cannot be determined as the first two points atre of equal vale');
            end
            fold_indexes=[];
            for j=2: size(params_LCC,2)
               
               if  dir==1 && params_LCC(j) < params_LCC(j-1)
                  fold_indexes =  [fold_indexes, j-1];
                  dir= -1;
               elseif  dir== -1 && params_LCC(j) > params_LCC(j-1)
                  fold_indexes =  [fold_indexes, j-1];
                  dir= 1;
               end
            end

            for j=min(fold_indexes)+1:max(fold_indexes)
                if params_LCC(j-1) < param_target && params_LCC(j) >= param_target
                    threshold_indexes=[j-1, j];
                end
            end
            
            % Set various constants
            nu = size(obj_CA.results.u,1);
            npts = size(obj_CA.continuationResults.y,2);             
            ntst = obj_CA.continuationResults.s(1,1).data.ntst;
            ncol = obj_CA.continuationResults.s(1,1).data.ncol;
            nLC_rows = size(obj_CA.continuationResults.y,1);
                           
            % Extract state variables for each limit cycle and calculate depth to water table (DTWT) and theta.            
            for k=threshold_indexes
                for i=1:ntst;
                    for j=1:ncol;
                        range0 = 2*nu*( (i-1)*ncol +j-1)+1 : 2 : 2*nu*( (i-1)*ncol +j)-1; 
                        
                        % Set model state variables
                        obj_CA.stateVariables.S = obj_CA.continuationResults.y(range0,k);
                        obj_CA.stateVariables.m = obj_CA.continuationResults.y(range0+1,k);
                        obj_CA.stateVariables.t = 1;

                        % Calculate depth to water table (DTWT) and theta
                        setVariables(obj_CA);
                         
                        % Extract depth to water table and theta.
                        DTWT_theta(range0,k-min(threshold_indexes)+1) = obj_CA.variables.DBNS;
                        DTWT_theta(range0+1,k-min(threshold_indexes)+1) = obj_CA.variables.theta;                                                 
                    end
                end  
            end
            
            % Extract data for limit cycle at node location
            range0=[];
            for i=1: nu
                u( 2*(i-1)+1,1) = obj_CA.results.u(i);
                u( 2*(i-1)+2,1) = obj_CA.results.u(i);
            end
            u = repmat( u, ntst*ncol, 1);
            
            for j=1:length(nodesToExtract)
                range0 = find(u == nodesToExtract(j)); 
                
                DTWT_theta_atU = DTWT_theta(range0, :);
                
                DTWT = DTWT_theta_atU(1:2:end,:);
                theta= DTWT_theta_atU(2:2:end,:);
                DTWT = DTWT(:,1) + (DTWT(:,2)-DTWT(:,1)) ...
                    ./(params_LCC(threshold_indexes(2)) - params_LCC(threshold_indexes(1))) ...
                    .* (param_target - params_LCC(threshold_indexes(1)));

                theta = theta(:,1) + (theta(:,2)-theta(:,1)) ...
                    ./(params_LCC(threshold_indexes(2)) - params_LCC(threshold_indexes(1))) ...
                    .* (param_target - params_LCC(threshold_indexes(1)));
                
                DTWT_thresholds(j,:) = [nodesToExtract(j), min(DTWT), max(DTWT)];
                SM_thresholds(j,:)= [nodesToExtract(j), min(theta), max(theta)];
            end
       end
       
       function out = getDerivatives_Continuation(obj_CA)
% GETDERIVATIVES_CONTINUATION  Nested set of functions to calculate the required derivatives and Jacobian matrix.         

            out{1} = [];
            out{2} = @fun_eval;
            out{3} = @getJacobian;             
            out{4} = [];
            out{5} = [];
            out{6} = [];
            out{7} = [];
            out{8} = [];
            out{9} = [];
            % set out options
            %   for starters ignore bvp jacobian and paramm jacobian
            % --------------------------------------------------------------------------
            function dydt = fun_eval(t,y, varargin)
                global objGlobal_LC;
                if isempty(t)
                    t=10;
                else
                    t = t + 10;
                end
                
                if nargin == 3
                    ParamValue1  = 1;
                    ParamValue2 = varargin{1};
                elseif nargin == 4
                    ParamValue1 = varargin{1};
                    ParamValue2 = varargin{2};
                else
                     error('The number of variable arguments for the parameter values was expected to be 1 or 2.');
                end
                
                % Assign parameter value
                if obj_CA.continuationDoLimitCycle
                    setParameterValue(objGlobal_LC, objGlobal_LC.continuationParameterName, ParamValue2);
                else
                    if strcmp(objGlobal_LC.continuationParameterName(end),'P')
                        obj_CA.climateData.data.P(:,2) = ParamValue2;

                    elseif strcmp(objGlobal_LC.continuationParameterName(end-1:end),'ET')
                        obj_CA.climateData.data.ET(:,2) = ParamValue2;
                    end
                    
                end


                % Update model constants
                setGeometryVariables(objGlobal_LC);   
                
                dydt = getDerivatives_DAE(t, y, objGlobal_LC);
            end

           % Calculate model Jacobian at a given time point
           function dFdy = getJacobian(t,y, varargin)

                % dFdu Modified odenumjac.m so that the delta function is constant
                % Delta function was required to be constant so that problems in the limit
                % cycle divering from prior convergences could be investigated.
%                global objGlobal_LC;
                global eds lds cds;           
                
                if nargin == 5
                % Doing 1-parameter equilbrium continuation    
                    ParamValue1 =1;
                    ParamValue2 = varargin{1};
                    g = varargin{2};
                    S = varargin{3};                    
                elseif nargin == 6
                % Doing 2-parameter equilbrium continuation or 1 parameter
                % limit cycle continuation.
                    ParamValue1 = varargin{1};
                    ParamValue2 = varargin{2};
                    g = varargin{3};
                    S = varargin{4};                                    
                else
                    error('The number of variable arguments was expected to be 3 or 4.');
                end
                
                ny = size(y,1);  
                ng = max(g);
                nF=ny;
                one2ny = (1:ny)';

                nzcols = find(g > 0); % g==0 for all-zero columns in sparsity pattern  
                JacIndex = (g(nzcols)-1)*ny + nzcols;

                [JacIndex_ii JacIndex_jj] = find(S);
                if ~isempty(JacIndex_ii)
                    JacIndex_ii = JacIndex_ii(:); % ensure that JacIndex is a column vector (S could be a row vector)  
                end

                % get del into an extact number
                %del = 1e-10i;
                del = cds.options.Increment;        
                y_typical = 1000;
                if isreal(del);
                    del=max(eps(1), del);
                    temp = y_typical + del;
                    del = temp - y_typical;      
                else     
                    del=max(1i*eps(y_typical), del);
                    temp = y_typical*1i + del;
                    del = temp - y_typical*1i;
                end

                if del==0;
                    warning('Jacobian del==0. May be due to machine precision rounding');
                end
                
                doJacRichardExtrapolation=false;
                if doJacRichardExtrapolation
                   
                    % Ridder's method for polynomial extrapolation of centre weighted finite
                    % diff derivatives.
                    % From Press et al , 2007, Numerical Recipes, p 230 secion 5.7, editon 3.
                    ntab = 10;
                    con = 1.4;
                    con2 = con^2;
                    safe=2;
                    big = 1e29;
                    a = cell(ntab, ntab);
                    h=del;

                    hh=h;

                    % get est for h                  
                    a{1,1} = getFiniteDiffJac(y, hh);

                    err=big;
                    ii=2;
                    %maxii=5;
                    while ii<=ntab
                    %for ii=iiInitial:maxii
                        hh = hh/con;

                        % try smaller step of h
                        a{1,ii} = getFiniteDiffJac(y, hh);
                        
                        if norm(full(a{1,ii}))==0;
                           dFdy = a{1,ii};
                           return
                        end
                        fac=con2;

                        % compute extrapolations
                        jj=2;
                        while jj<=ii
                           a{jj,ii}  = (a{jj-1,ii}.*fac - a{jj-1,ii-1}) ./ (fac - 1);
                           fac = con2*fac;
                           errt = norm( full(max( abs(a{jj,ii} - a{jj-1,ii}), abs(a{jj,ii} - a{jj-1,ii-1}))) );

                           if errt < err;
                               err=errt;
                               dFdy = a{jj,ii};

                               if any(any(imag(dFdy)~=0));
                                   display('dFdy term is imag within Ridders Algorithm');
                                   error(' ');
                               end
                           end
                           jj=jj+1;
                        end

                        if norm( full(abs(a{ii,ii} - a{ii-1,ii-1})) ) > safe*err;
                            return;
                        end
                        ii=ii+1;
                    end                       
                    
               else    
                    % Calc a single centre weighted finite difference approx
                    dFdy = getFiniteDiffJac(y, del);
               end
               function dfdy = getFiniteDiffJac(y, del)
                    ydel = y(:,ones(1,ng));
                    ydel_pls = ydel;
                    ydel_neg = ydel;  
                    ydel_pls(JacIndex) = ydel_pls(JacIndex) + del;
                    ydel_neg(JacIndex) = ydel_neg(JacIndex) - del;     

                    if ~isempty(lds)
                        Fdel = feval(lds.func, t,[ydel_neg, ydel_pls], ParamValue1, ParamValue2);        
                    elseif ~isempty(eds)
                        Fdel = feval(eds.func, t,[ydel_neg, ydel_pls], ParamValue1, ParamValue2);        
                    else
                        error('Was expected a model function name in global variables eds (for equil. cont.) or lds (limit-clcye cont.)');
                    end
                    
                    Fdiff = Fdel(:,ng+1:end) - Fdel(:,1:ng);
                    ndel=2;

                    Fdiff = sparse(JacIndex_ii, JacIndex_jj ,Fdiff((g(JacIndex_jj)-1)*nF + JacIndex_ii),nF,ny);

                    dfdy = real(Fdiff * sparse(one2ny, one2ny, 1./(ndel*del), ny, ny));                                      
               end
           end

       end
                  
    
   end
end 