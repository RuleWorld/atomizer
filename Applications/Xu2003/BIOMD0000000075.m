% This file works with OCTAVE and is automatically generated with 
% the System Biology Format Converter (http://sbfc.sourceforge.net/)
% from an SBML file.
% To run this file with Matlab you must edit the comments providing
% the definition of the ode solver and the signature for the 
% xdot function.
%
% The conversion system has the following limitations:
%  - You may have to re order some reactions and Assignment Rules definition
%  - Delays are not taken into account
%  - You should change the lsode parameters (start, end, steps) to get better results
%

%
% Model name = Xu2003 - Phosphoinositide turnover
%
% is http://identifiers.org/biomodels.db/MODEL3095606944
% is http://identifiers.org/biomodels.db/BIOMD0000000075
% isDescribedBy http://identifiers.org/pubmed/12771127
% isDerivedFrom http://identifiers.org/pubmed/10579714
% isDerivedFrom http://identifiers.org/pubmed/10866945
%


function main()
%Initial conditions vector
	x0=zeros(13,1);
	x0(1) = 0.0;
	x0(2) = 0.0;
	x0(3) = 142857.0;
	x0(4) = 1.0;
	x0(5) = 0.0;
	x0(6) = 4000.0;
	x0(7) = 2857.0;
	x0(8) = 2000.0;
	x0(9) = 0.0;
	x0(10) = 0.0;
	x0(11) = 100.0;
	x0(12) = 0.0;
	x0(13) = 96.32;


% Depending on whether you are using Octave or Matlab,
% you should comment / uncomment one of the following blocks.
% This should also be done for the definition of the function f below.
% Start Matlab code
	tspan=[0:0.01:100];
	opts = odeset('AbsTol',1e-3);
	[t,x]=ode23tb(@f,tspan,x0,opts);
% End Matlab code

% Start Octave code
% 	t=linspace(0,100,100);
% 	x=lsode('f',x0,t);
% End Octave code


	plot(t,x(:,6:7));
end



% Depending on whether you are using Octave or Matlab,
% you should comment / uncomment one of the following blocks.
% This should also be done for the definition of the function f below.
% Start Matlab code
function xdot=f(t,x)
% End Matlab code

% Start Octave code
%function xdot=f(x,t)
% End Octave code

% Compartment: id = Extracellular, name = Extracellular, constant
	compartment_Extracellular=0.277777777777778;
% Compartment: id = PM, name = PM, constant
	compartment_PM=0.5555555555556;
% Compartment: id = Cytosol, name = Cytosol, constant
	compartment_Cytosol=1.0;
% Compartment: id = NM, name = NM, constant
	compartment_NM=0.111111111111111;
% Compartment: id = Nucleus, name = Nucleus, constant
	compartment_Nucleus=0.111111111111111;
% Parameter:   id =  KMOLE, name = KMOLE
	global_par_KMOLE=0.00166112956810631;
% Parameter:   id =  PIP_basal_PIPSyn, name = PIP_basal_PIPSyn
	global_par_PIP_basal_PIPSyn=2857.0;
% Parameter:   id =  kBasalSynPIP_PIPSyn, name = kBasalSynPIP_PIPSyn
	global_par_kBasalSynPIP_PIPSyn=0.0055;
% Parameter:   id =  kStimSynPIP_PIPSyn, name = kStimSynPIP_PIPSyn
	global_par_kStimSynPIP_PIPSyn=0.019;
% Parameter:   id =  tauPIPsyn_PIPSyn, name = tauPIPsyn_PIPSyn
	global_par_tauPIPsyn_PIPSyn=0.05;
% Parameter:   id =  PIPsyndecay_PIPSyn, name = PIPsyndecay_PIPSyn
	global_par_PIPsyndecay_PIPSyn=1.0;
% Parameter:   id =  Ratebasal_PIPsyn_PIPSyn, name = Ratebasal_PIPsyn_PIPSyn
% Parameter:   id =  Ratestim_PIPsyn_PIPSyn, name = Ratestim_PIPsyn_PIPSyn
% Parameter:   id =  tau0_PLCact, name = tau0_PLCact
	global_par_tau0_PLCact=0.05;
% Parameter:   id =  stimdecay_PLCact, name = stimdecay_PLCact
	global_par_stimdecay_PLCact=1.0;
% Parameter:   id =  signal_PLCact, name = signal_PLCact
% Parameter:   id =  kf_PIP2PH_PIP2_PH, name = kf_PIP2PH_PIP2_PH
	global_par_kf_PIP2PH_PIP2_PH=0.12;
% Parameter:   id =  KdPIP2PH_PIP2_PH, name = KdPIP2PH_PIP2_PH
	global_par_KdPIP2PH_PIP2_PH=2.0;
% Parameter:   id =  kr_PIP2PH_PIP2_PH, name = kr_PIP2PH_PIP2_PH
% Parameter:   id =  kStimSynPIP2_PIP2Syn, name = kStimSynPIP2_PIP2Syn
	global_par_kStimSynPIP2_PIP2Syn=0.92;
% Parameter:   id =  tauPIP2syn_PIP2Syn, name = tauPIP2syn_PIP2Syn
	global_par_tauPIP2syn_PIP2Syn=0.05;
% Parameter:   id =  PIP2syndecay_PIP2Syn, name = PIP2syndecay_PIP2Syn
	global_par_PIP2syndecay_PIP2Syn=1.0;
% Parameter:   id =  PIP2_basal_PIP2Syn, name = PIP2_basal_PIP2Syn
	global_par_PIP2_basal_PIP2Syn=4000.0;
% Parameter:   id =  kBasalSynPIP2_PIP2Syn, name = kBasalSynPIP2_PIP2Syn
	global_par_kBasalSynPIP2_PIP2Syn=0.048;
% Parameter:   id =  Rate_PIP2Synbasal_PIP2Syn, name = Rate_PIP2Synbasal_PIP2Syn
% Parameter:   id =  Rate_PIP2SynStim_PIP2Syn, name = Rate_PIP2SynStim_PIP2Syn
% Parameter:   id =  kf_IP3PH_IP3_PHGFP, name = kf_IP3PH_IP3_PHGFP
	global_par_kf_IP3PH_IP3_PHGFP=10.0;
% Parameter:   id =  KdIP3PH_IP3_PHGFP, name = KdIP3PH_IP3_PHGFP
	global_par_KdIP3PH_IP3_PHGFP=2.0;
% Parameter:   id =  kr_IP3PH_IP3_PHGFP, name = kr_IP3PH_IP3_PHGFP
% assignmentRule: variable = Ratebasal_PIPsyn_PIPSyn
	global_par_Ratebasal_PIPsyn_PIPSyn=piecewise(0.581*global_par_kBasalSynPIP_PIPSyn*(-1+exp((global_par_PIP_basal_PIPSyn+(-x(7)))*1/global_par_PIP_basal_PIPSyn)), x(7) < global_par_PIP_basal_PIPSyn, 0);
% assignmentRule: variable = Ratestim_PIPsyn_PIPSyn
	global_par_Ratestim_PIPsyn_PIPSyn=piecewise(global_par_kStimSynPIP_PIPSyn*exp(-(t+(-global_par_tauPIPsyn_PIPSyn))*1/global_par_PIPsyndecay_PIPSyn), t > global_par_tauPIPsyn_PIPSyn, 0);
% assignmentRule: variable = signal_PLCact
	global_par_signal_PLCact=piecewise(exp(-(t+(-global_par_tau0_PLCact))*1/global_par_stimdecay_PLCact), t > global_par_tau0_PLCact, 0);
% assignmentRule: variable = kr_PIP2PH_PIP2_PH
	global_par_kr_PIP2PH_PIP2_PH=global_par_kf_PIP2PH_PIP2_PH*global_par_KdPIP2PH_PIP2_PH;
% assignmentRule: variable = Rate_PIP2Synbasal_PIP2Syn
	global_par_Rate_PIP2Synbasal_PIP2Syn=piecewise(0.581*global_par_kBasalSynPIP2_PIP2Syn*(-1+exp((global_par_PIP2_basal_PIP2Syn+(-x(6)))*1/global_par_PIP2_basal_PIP2Syn)), x(6) < global_par_PIP2_basal_PIP2Syn, 0);
% assignmentRule: variable = Rate_PIP2SynStim_PIP2Syn
	global_par_Rate_PIP2SynStim_PIP2Syn=piecewise(global_par_kStimSynPIP2_PIP2Syn*exp(-(t+(-global_par_tauPIP2syn_PIP2Syn))*1/global_par_PIP2syndecay_PIP2Syn), t > global_par_tauPIP2syn_PIP2Syn, 0);
% assignmentRule: variable = kr_IP3PH_IP3_PHGFP
	global_par_kr_IP3PH_IP3_PHGFP=global_par_kf_IP3PH_IP3_PHGFP*global_par_KdIP3PH_IP3_PHGFP;

% Reaction: id = PIPSyn, name = PIPSyn	% Local Parameter:   id =  I, name = I
	reaction_PIPSyn_I=0.0;

	reaction_PIPSyn=(global_par_Ratebasal_PIPsyn_PIPSyn+global_par_Ratestim_PIPsyn_PIPSyn)*x(3)*compartment_PM;

% Reaction: id = PIP2_hyd, name = PIP2_hyd	% Local Parameter:   id =  I, name = I
	reaction_PIP2_hyd_I=0.0;
	% Local Parameter:   id =  k_PIP2hyd, name = k_PIP2hyd
	reaction_PIP2_hyd_k_PIP2hyd=2.4;

	reaction_PIP2_hyd=reaction_PIP2_hyd_k_PIP2hyd*x(6)*x(12)*compartment_PM;

% Reaction: id = PLCact, name = PLCact	% Local Parameter:   id =  I, name = I
	reaction_PLCact_I=0.0;
	% Local Parameter:   id =  KfPLCact, name = KfPLCact
	reaction_PLCact_KfPLCact=5.0E-4;
	% Local Parameter:   id =  krPLCact, name = krPLCact
	reaction_PLCact_krPLCact=0.1;

	reaction_PLCact=(reaction_PLCact_KfPLCact*x(11)*x(4)*global_par_signal_PLCact+(-reaction_PLCact_krPLCact*x(12)))*compartment_PM;

% Reaction: id = PIP2_PH_hyd, name = PIP2_PH_hyd	% Local Parameter:   id =  I, name = I
	reaction_PIP2_PH_hyd_I=0.0;
	% Local Parameter:   id =  k_PIP2PHhyd, name = k_PIP2PHhyd
	reaction_PIP2_PH_hyd_k_PIP2PHhyd=0.0;

	reaction_PIP2_PH_hyd=reaction_PIP2_PH_hyd_k_PIP2PHhyd*x(12)*x(1)*compartment_PM;

% Reaction: id = PIP2_PH, name = PIP2_PH	% Local Parameter:   id =  I, name = I
	reaction_PIP2_PH_I=0.0;

	reaction_PIP2_PH=(global_par_kf_PIP2PH_PIP2_PH*0.00166112956810631*x(2)*x(6)+(-global_par_kr_PIP2PH_PIP2_PH*x(1)))*compartment_PM;

% Reaction: id = IP3deg, name = IP3deg	% Local Parameter:   id =  kIP3deg, name = kIP3deg
	reaction_IP3deg_kIP3deg=0.08;
	% Local Parameter:   id =  IP3_basal, name = IP3_basal
	reaction_IP3deg_IP3_basal=0.16;

	reaction_IP3deg=reaction_IP3deg_kIP3deg*(0.00166112956810631*x(13)+(-reaction_IP3deg_IP3_basal))*compartment_Cytosol*1*1/global_par_KMOLE;

% Reaction: id = PIP2Syn, name = PIP2Syn	% Local Parameter:   id =  I, name = I
	reaction_PIP2Syn_I=0.0;

	reaction_PIP2Syn=(global_par_Rate_PIP2Synbasal_PIP2Syn+global_par_Rate_PIP2SynStim_PIP2Syn)*x(7)*compartment_PM;

% Reaction: id = IP3_PHGFP, name = IP3-PHGFP
	reaction_IP3_PHGFP=(global_par_kf_IP3PH_IP3_PHGFP*0.00166112956810631*x(13)*0.00166112956810631*x(2)+(-global_par_kr_IP3PH_IP3_PHGFP*0.00166112956810631*x(5)))*compartment_Cytosol*1*1/global_par_KMOLE;

	xdot=zeros(13,1);
	
% Species:   id = PIP2_PHGFP_PM, name = PIP2_PHGFP_PM, affected by kineticLaw
	xdot(1) = (1/(compartment_PM))*((-1.0 * reaction_PIP2_PH_hyd) + ( 1.0 * reaction_PIP2_PH));
	
% Species:   id = PH_GFP_Cyt, name = PH_GFP_Cyt, affected by kineticLaw
	xdot(2) = (1/(compartment_Cytosol))*(( 1.0 * reaction_PIP2_PH_hyd) + (-1.0 * reaction_PIP2_PH) + (-1.0 * reaction_IP3_PHGFP));
	
% Species:   id = PI_PM, name = PI_PM, affected by kineticLaw
	xdot(3) = (1/(compartment_PM))*((-1.0 * reaction_PIPSyn));
	
% Species:   id = stim_PM, name = stim_PM
%WARNING speciesID: stim_PM, constant= false  , boundaryCondition = true but is not involved in assignmentRule, rateRule or events !
	xdot(4) = 0.0;
	
% Species:   id = IP3_PHGFP_Cyt, name = IP3_PHGFP_Cyt, affected by kineticLaw
	xdot(5) = (1/(compartment_Cytosol))*(( 1.0 * reaction_IP3_PHGFP));
	
% Species:   id = PIP2_PM, name = PIP2_PM, affected by kineticLaw
	xdot(6) = (1/(compartment_PM))*((-1.0 * reaction_PIP2_hyd) + (-1.0 * reaction_PIP2_PH) + ( 1.0 * reaction_PIP2Syn));
	
% Species:   id = PIP_PM, name = PIP_PM, affected by kineticLaw
	xdot(7) = (1/(compartment_PM))*(( 1.0 * reaction_PIPSyn) + (-1.0 * reaction_PIP2Syn));
	
% Species:   id = DAG_PM, name = DAG_PM, affected by kineticLaw
	xdot(8) = (1/(compartment_PM))*(( 1.0 * reaction_PIP2_hyd) + ( 1.0 * reaction_PIP2_PH_hyd));
	
% Species:   id = hv_Cytosol, name = hv_Cytosol
%WARNING speciesID: hv_Cytosol, constant= false  , boundaryCondition = true but is not involved in assignmentRule, rateRule or events !
	xdot(9) = 0.0;
	
% Species:   id = IP3X_Cytosol, name = IP3X_Cytosol
%WARNING speciesID: IP3X_Cytosol, constant= false  , boundaryCondition = true but is not involved in assignmentRule, rateRule or events !
	xdot(10) = 0.0;
	
% Species:   id = PLC_PM, name = PLC_PM, affected by kineticLaw
	xdot(11) = (1/(compartment_PM))*((-1.0 * reaction_PLCact));
	
% Species:   id = PLC_act_PM, name = PLC_act_PM, affected by kineticLaw
	xdot(12) = (1/(compartment_PM))*(( 1.0 * reaction_PLCact));
	
% Species:   id = IP3_Cyt, name = IP3_Cyt, affected by kineticLaw
	xdot(13) = (1/(compartment_Cytosol))*(( 1.0 * reaction_PIP2_hyd) + ( 1.0 * reaction_PIP2_PH_hyd) + (-1.0 * reaction_IP3deg) + (-1.0 * reaction_IP3_PHGFP));
end

% adding few functions representing operators used in SBML but not present directly 
% in either matlab or octave. 
function z=pow(x,y),z=x^y;end
function z=root(x,y),z=y^(1/x);end
function z = piecewise(varargin)
	numArgs = nargin;
	result = 0;
	foundResult = 0;
	for k=1:2: numArgs-1
		if varargin{k+1} == 1
			result = varargin{k};
			foundResult = 1;
			break;
		end
	end
	if foundResult == 0
		result = varargin{numArgs};
	end
	z = result;
end


