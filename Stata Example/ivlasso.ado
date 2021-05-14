*! ivlasso 1.0.06 10feb2018
*! authors aa/cbh/ms
* ivlasso:		main program
* s_ivparse:	Mata basic parse of IV-type lists (endog, dexog, exexog, xctrl)
* _pdsparse:	further parsing into partial, notpen, amelioriation set, separately by X and Z vars
* _ivlasso:		estimation program

* Updates (release date):
* 1.0.05  (30jan2018)
*         First public release.
*         Display name for low-dim IVs no longer prefixed with "rho_".
*         Display name for endogenous d no longer prefixed with "hat_".
*         Display dots and title when estimating cset using user-supplied sup-score grid.
*         Display dot at end of a single loop.
*         Added ssomitgrid option. Changed default sup-score method to abound.
*         Promoted to require version 13 or higher.
*         Fixed bug (missing touse) in generation of post-lasso resids with aset option in SelectControls.
* 1.0.06  (xxx)
*         Fixed bug in SelectControls - was reintroducing nocons in call to rlasso when partial(.) also specified.
*         Support for Sergio Correia's FTOOLS FE transform (if installed).

program define ivlasso, eclass sortpreserve
	syntax [anything] [if] [in] ,		///
		[								///
		cmdname(name)					///
		VERsion							///
		* ]
		
	version 13
	local lversion 1.0.05
	if "`version'" != "" {							//  Report program version number, then exit.
		di in gr "`lversion'"
		ereturn clear
		ereturn local version `lversion'
		exit
	}

	if "`cmdname'"=="" {
		local cmdname ivlasso						//  not called by pdslasso so this is an ivlasso estimation
	}
	
	if ~replay() {									//  not replay so estimate

		// check for rlasso (eclass-command) only when estimating anew
		cap rlasso, version
		if _rc != 0 {
			di as err "Error: `cmdname' requires rlasso to run"
			di as err "To install, type " _c
			di in smcl "{stata ssc install rlasso :ssc install rlasso}"
			exit 601
		}
		mata: s_ivparse("`anything'")

// Note that these varlists may have abbreviations, wildcards, etc.
		local depvar	`s(depvar)'
		local dendog	`s(dendog)'
		local dexog		`s(dexog)'
		local xctrl		`s(xctrl)'
		local exexog	`s(exexog)'

		*** Do IV or OLS PDS
		_ivlasso `depvar' `if' `in',		///
			dendog(`dendog')				///
			dexog(`dexog')					///
			xctrl(`xctrl')					///
			exexog(`exexog')				///
			cmdname(`cmdname')				///
			`options'
		
	}
	else if "`e(cmd)'"~="ivlasso" & "`e(cmd)'"~="pdslasso" {	//  replay, so check that ivlasso/pdslasso results exist
		di as err "last estimates not found"
		exit 301
	}

	if "`e(firstvar)'" == "" {									//  display ivlasso/pdslasso main results
		DisplayResults											//  this trap is needed because est replay is used to replay
	}															//  the first stage results
	else {
		ereturn di
	}

end		// end wrapper 

// parses basic IV varlists into PDS-type varlists
program define _pdsparse, sclass
	syntax,						/// varlists passed as strings to stop b/n/o operators being inserted
		depvar(string)			///
	[							///
		dexog(string)			///
		dendog(string)			///
		xctrl(string)			///
		exexog(string)			///
		notpen(string)			///
		partial(string)			///
		aset(string)			///
		debugflag(int 0)		///
	]

	*** Initial syntax checks - overlapping lists
	// check that depvar is not in other lists
	foreach vlist in dendog dexog notpen partial aset {
		local checklist	: list depvar & `vlist'
		local checknum	: word count `checklist'
		if `checknum' {
			di as err "syntax error - cannot also include dependent variable `checklist' in `vlist'"
			exit 198
		}
	}
	// check that dexog vars are not in other lists
	foreach vlist in depvar dendog notpen partial aset {
		local checklist	: list dexog & `vlist'
		local checknum	: word count `checklist'
		if `checknum' {
			di as err "syntax error - cannot also include exogenous variable `checklist' in `vlist'"
			exit 198
		}
	}
	// check that dendog vars are not in other lists
	foreach vlist in depvar dexog notpen partial aset {
		local checklist	: list dendog & `vlist'
		local checknum	: word count `checklist'
		if `checknum' {
			di as err "syntax error - cannot also include endogenous variable `checklist' in `vlist'"
			exit 198
		}
	}
	// check that notpen vars are not in other lists
	foreach vlist in depvar dendog dexog partial aset {
		local checklist	: list notpen & `vlist'
		local checknum	: word count `checklist'
		if `checknum' {
			di as err "syntax error - cannot also include unpenalized `checklist' in `vlist'"
			exit 198
		}
	}
	// check that partial vars are not in other lists
	foreach vlist in depvar dendog dexog notpen aset {
		local checklist	: list partial & `vlist'
		local checknum	: word count `checklist'
		if `checknum' {
			di as err "syntax error - cannot also include partialled-out `checklist' in `vlist'"
			exit 198
		}
	}
	// check that xctrl vars are not in other lists
	foreach vlist in depvar dendog dexog exexog {
		local checklist	: list xctrl & `vlist'
		local checknum	: word count `checklist'
		if `checknum' {
			di as err "syntax error - cannot also include control variable `checklist' in `vlist'"
			exit 198
		}
	}
	// check that IVs are not in other lists
	foreach vlist in depvar dendog dexog xctrl {
		local checklist	: list exexog & `vlist'
		local checknum	: word count `checklist'
		if `checknum' {
			di as err "syntax error - cannot also include IV `checklist' in `vlist'"
			exit 198
		}
	}
	*

	*** separate out the X controls and Z IVs
	// 2 types: (1) high-dim (penalized); (2) unpenalized or partialled-out
	// xhighdim starts as xctrl; zhighdim starts as exexog
	local xnotpen	: list notpen	- exexog	//  = explicit notpen Xs
	local xhighdim	: list xctrl	- notpen	//  = Xs not specified as notpen
	local znotpen	: list notpen	- xctrl		//  = explicit notpen Zs
	local zhighdim	: list exexog	- notpen	//  = Zs not specified as notpen
	// here use xhighdim and zhighdim already started
	local xpartial	: list partial	- exexog	//  = explicit partial Xs
	local xhighdim	: list xhighdim	- partial	//  = Xs not specified as partial
	local zpartial	: list partial	- xctrl		//  = explicit partial Zs
	local zhighdim	: list zhighdim	- partial	//  = Zs not specified as partial
	// checks
	if "`notpen'"~="" {
		local check		: list notpen	- xctrl		//  check that notpen has no extraneous variables
		local check		: list check	- exexog
		local checknum	: word count `check'
		if `checknum' {
			di as err "error: variables `check' appear in notpen(.) but not as controls or instruments"
			exit 198
		}
	}
	if "`partial'"~="" {
		local check		: list partial	- xctrl		//  check that partial has no extraneous variables
		local check		: list check	- exexog
		local checknum	: word count `check'
		if `checknum' {
			di as err "error: variables `check' appear in partial(.) but not as controls or instruments"
			exit 198
		}
	}
	*
	
	*** separate out the amelioration sets
	local xaset		: list xhighdim	- aset		//  inverse xaset
	local xaset		: list xhighdim - xaset		//  reverse back to get xaset
	local zaset		: list zhighdim	- aset		//  inverse zaset
	local zaset		: list zhighdim - zaset		//  reverse back to get zaset
	local checklist	: list aset - xhighdim
	local checklist	: list checklist - zhighdim
	local checknum	: word count `checklist'
	if `checknum' {
		di as err "syntax error - `checklist' in aset(.) but not in list of high-dim controls or instruments"
		exit 198
	}
	*

	if `debugflag' {
		di
		di as text "_pdsparse parsing:"
		local listlist depvar dexog dendog xhighdim xnotpen xpartial xaset zhighdim znotpen zpartial zaset
		foreach list of local listlist {
			di as res "`list': ``list''"
		}
	}

	// depvar, dendog, dexog are automatically low-dim/model variables
	// but return anyway
	sreturn local depvar_o		`depvar'
	sreturn local dendog_o		`dendog'
	sreturn local dexog_o		`dexog'
	sreturn local xhighdim_o	`xhighdim'
	sreturn local xnotpen_o		`xnotpen'
	sreturn local xaset_o		`xaset'
	sreturn local xpartial_o	`xpartial'
	sreturn local zhighdim_o	`zhighdim'
	sreturn local znotpen_o		`znotpen'
	sreturn local zaset_o		`zaset'
	sreturn local zpartial_o	`zpartial'

end

// Main program
program define _ivlasso, eclass sortpreserve
	syntax varlist(numeric fv ts min=1 max=1) [if] [in],	///
		[													///
		dexog(varlist numeric fv ts)						/// low-dim set of exog regressors of interest
		xctrl(varlist numeric fv ts)						/// all controls incl notpen/partial
		exexog(varlist numeric fv ts)						/// all IVs incl notpen/partial
		dendog(varlist numeric fv ts)						/// low-dim set of endog regressors of interest
		PNOTPen(varlist numeric fv ts)						/// unpenalized variables
		partial(string)										/// string so that list can contain "_cons"
		aset(varlist numeric fv ts)							/// amelioration set
															///
		fe													/// requires data to be xtset
															///
		NOIsily												///
		NOCONStant											///
		CLuster(varlist max=1)								/// shared option with lasso and IV estimation
		Robust												/// shared option with lasso and IV estimation
		sqrt												/// shared option with lasso and IV estimation
		first												/// report first-stage regressions
		idstats												/// report weak ID stats
		post(string)										/// which results to ereturn post
															///
		sscset												/// report sup-score CIs
		ssgamma(real 0.05)									/// gamma (significance) for sup-score test (default 5%)
		ssgridmin(numlist missingok)						/// gridmin for sup-score test
		ssgridmax(numlist missingok)						/// gridmax for sup-score test
		ssgridpoints(integer 100)							/// default gridpoints for sup-score test
		ssgridmat(name)										/// optional grid for sup-score test
		ssomitgrid											/// supress display of user-supplied grid
		ssmethod(name)										/// simulate, abound, select; default = abound
															///
		rlasso												/// store rlasso results?
		RLASSO0(name)										/// store rlasso results with prefix `name'
		LOPTions(string)									/// options passed to rlasso
		IVOPTions(string)									/// options passed to IV or OLS estimation
		cmdname(name)										/// ivlasso or pdslasso
		debug												/// triggers debugging output and calcs
		fvstrip												/// alternative undocumented parsing of factor vars
		NOFTOOLS											/// provisional
		]

	*** rlasso-specific
	//  to distinguish between lasso2 and rlasso treatment of notpen,
	//  rlasso option is called pnotpen
	//  rename to notpen here and at end of program save macros as pnotpen
	//  temporary measure until lasso2 and rlasso code is merge
	local notpen	`pnotpen'
	*
	*** debug mode; create flag
	local debugflag	=("`debug'"~="")
	*

	*** post option - which results to post
	// default is PDS; mixed case allowed, convert to lower case.
	local post	=strlower("`post'")
	if "`post'"=="" {
		local post	pds
	}
	else if "`post'"~="pds" & "`post'"~="lasso" & "`post'"~="plasso" {
		di as err "error: post(.) must be one of pds, lasso or plasso"
		exit 198
	}
	*
		
	*** FEs. Change flag to 1/0.
	// Get panel id
	local feflag=("`fe'"~="")
	if `feflag' {
		cap _xt
		if _rc ~= 0 {
			di as err "Error: fe option requires data to be xtset"
			exit 459
		}
		else {
			local ivar	`r(ivar)'
		}
		local noconstant	noconstant					//  if fe, nocons is automatic
	}
	*
	
	*** Shared options
	// robust => ivoptions(robust) + loptions(robust)
	if "`robust'" ~= "" {
		local ivoptions		`ivoptions' robust
		local loptions		`loptions' robust
	}
	// cluster
	if "`cluster'" ~= "" {
		local ivoptions		`ivoptions' cluster(`cluster')
		local loptions		`loptions' cluster(`cluster')
	}
	// sqrt lasso
	if "`sqrt'" ~= "" {
		local loptions		`loptions' sqrt
		local method		sqrt-lasso
	}
	else {
		local method		lasso
	}
	*

	// rename dep var
	local depvar		`varlist'
	*

	*** Record which observations have non-missing values
	marksample touse
	markout `touse' `depvar' `dendog' `dexog' `xctrl' `exexog' `cluster' `ivar'
	sum `touse' if `touse', meanonly		//  will sum weight var when weights are used
	local N		= r(N)
	*

	*** constant, partial, etc.
	local cons=("`noconstant'"=="")						//  if fe, then cons=0
	local partialflag	=("`partial'"~="")				//  =1 even if just cons being partialled out
	// "_cons" allowed as an argument to partial(.) - remove it
	local partial		: subinstr local partial "_cons" "", all word count(local pconscount)
	local notpen		: subinstr local notpen "_cons" "", all word count(local notpenconscount)
	if `cons' {
		local consname "_cons"							//  needed for display purposes
	}
	// if partial list has factor vars, will need to be replaced with tempvars
	cap _fv_check_depvar `partial'
	local partialfvflag	=(_rc==198)
	*
	
	*** Initialize N_clust=0; will be updated to >0 by estimating programs.
	local N_clust	=0		//  Refers to IV/OLS #clusters (not lasso #clusters)
	*
	
	if "`noisily'" == "" {
		local qui "qui"
	}
	*

	*** define various parameters/locals/varlists
	// remove duplicates from varlist
	// do this separately, list-by-list, so that factor var bases/omitteds/etc kept separate
	// _o list is vars with original names, omitteds dropped, b/n/o stripped off
	if "`fvstrip'"~="" {
		foreach vlist in depvar dendog dexog xctrl exexog notpen aset partial {
			if "``vlist''" ~= "" {												//  lists can be empty
				fvstrip ``vlist'' if `touse', expand dropomit
				local `vlist'_o			`r(varlist)'
				// check for duplicates has to follow expand
				local dups			: list dups `vlist'_o
				if "`dups'"~="" {
					di as text "Dropping duplicates: `dups'"
				}
				local `vlist'_o		: list uniq `vlist'_o
			}
		}
	}
	else {
		foreach vlist in depvar dendog dexog xctrl exexog notpen aset partial {
			if "``vlist''" ~= "" {												//  lists can be empty
				fvexpand ``vlist'' if `touse'
				local `vlist'_o			`r(varlist)'
				// check for duplicates has to follow expand
				local dups			: list dups `vlist'_o
				if "`dups'"~="" {
					di as text "Dropping duplicates: `dups'"
				}
				local `vlist'_o		: list uniq `vlist'_o
			}
		}
	}
	// Complete parsing
	_pdsparse,							///
				depvar(`depvar_o')		///
				dendog(`dendog_o')		///
				dexog(`dexog_o')		///
				xctrl(`xctrl_o')		///
				exexog(`exexog_o')		///
				notpen(`notpen_o')		///
				aset(`aset_o')			///
				partial(`partial_o')	/// note that partial gets special treatment - see below
				debugflag(`debugflag')
	// Y and endog automatically all low-dim so these are aliases
	// depvar_o, dendog_o, dexog_o already exist
	local varY_o		`depvar_o'			//  same as already-existing depvar_o
	// returned by PDS parser
	local xhighdim_o	`s(xhighdim_o)'
	local xnotpen_o		`s(xnotpen_o)'
	local xpartial_o	`s(xpartial_o)'
	local xaset_o		`s(xaset_o)'
	local zhighdim_o	`s(zhighdim_o)'
	local znotpen_o		`s(znotpen_o)'
	local zpartial_o	`s(zpartial_o)'
	local zaset_o		`s(zaset_o)'
	// and now set xpartial flag
	local xpartialflag	=("`xpartial_o'"~="") | `pconscount'				//  =1 even if just cons being partialled out
	*

	*** Create _t varlists: Y, X, Z, notpen, partial
	// _o list is vars with original names, omitteds dropped, b/n/o stripped off
	// _t list is temp vars if transform needed, original vars if not
	// _f list is _o list with FVs etc. replaced by temps using fvrevar
	// FEs or partialling-out FVs => everything has to be transformed
	// Partialling out Xs => everything except xpartial has to be transformed
	// Otherwise need temp vars only for TS and FV vars.
	// note not necessary to list aset vars since they already appear in highdim lists
	local allvars_o			`varY_o' `dexog_o' `dendog_o' `zhighdim_o' `znotpen_o' `zpartial_o' `xhighdim_o' `xnotpen_o' `xpartial_o'
	local dict_o			`allvars_o'
	local allbutxpartial_o	`varY_o' `dexog_o' `dendog_o' `zhighdim_o' `znotpen_o' `zpartial_o' `xhighdim_o' `xnotpen_o'
	if `feflag' {													//  everything needs to be transformed including xpartial
																	//  so create allvars_t and allvars_f
		local temp_ct : word count `allvars_o'
		mata: s_maketemps(`temp_ct')								//  uninitialized tempvars
		local allvars_t		`r(varlist)'
		foreach var in `allvars_o' {								//  we also need an _f list for allvars
			fvrevar `var' if `touse'								//  in case of factor vars or time-series operators
			local allvars_f		`allvars_f' `r(varlist)'			//  one-by-one because otherwise omitted base vars are set etc.
		}															//  fvrevar creates temps only when needed
		local dict_t		`allvars_t'								//  dictionary is identical to allvars lists
		local dict_f		`allvars_f'
		if `xpartialflag' {
			matchnames "`allbutxpartial_o'" "`dict_o'" "`dict_t'"
			local allbutxpartial_t	`r(names)'						//  needed in case of FE + partialling-out
		}
	}
	else if `xpartialflag' {										//  everything except xpartial_o needs to be transformed
		local temp_ct : word count `allbutxpartial_o'
		mata: s_maketemps(`temp_ct')								//  uninitialized tempvars excluding xpartial
		local allbutxpartial_t	`r(varlist)'
		foreach var in `allbutxpartial_o' {							//  we also need an _f list for allbutxpartial and xpartial
			fvrevar `var' if `touse'								//  in case of factor vars or time-series operators
			local allbutxpartial_f	`allbutxpartial_f' `r(varlist)'	//  one-by-one because otherwise omitted base vars are set etc.
		}															//  fvrevar creates temps only when needed
		foreach var in `xpartial_o' {								//  same for xpartial
			fvrevar `var' if `touse'								//  
			local xpartial_f	`xpartial_f' `r(varlist)'			//  
		}															//  
		local dict_t	`allbutxpartial_t' `xpartial_o'				//  xpartial_t = xpartial_o (untransformed)
		local dict_f	`allbutxpartial_f' `xpartial_f'
	}
	else {															//  no transformation needed but still need temps
		foreach var in `allvars_o' {								//  in case of factor vars or time-series operators
			fvrevar `var' if `touse'								//  one-by-one because otherwise omitted base vars are set etc.			
			local allvars_t		`allvars_t' `r(varlist)'			//  fvrevar creates temps only when needed
		}
		local allvars_f		`allvars_t'								//  and in this case the _f and _t lists are identical
		local dict_t		`allvars_t'								//  dictionary is identical to allvars lists
		local dict_f		`allvars_f'
	}
	// dictionary is now dict_o / dict_t; dict_f is dict_o but with fv and TS vars replaced with temps
	// create individual lists of _t variables with corresponding dictionary
	foreach vlist in varY dendog dexog xhighdim xnotpen xpartial xaset zhighdim znotpen zpartial zaset {
		matchnames "``vlist'_o'" "`dict_o'" "`dict_t'"
		local `vlist'_t		`r(names)'		//  corresponding tempnames of vars
	}
	*
	
	*** More useful flags
	local npxflag		=("`xnotpen_o'"~="")
	local npzflag		=("`znotpen_o'"~="")
	local hdzflag		=("`zhighdim_o'"~="")
	local hdxflag		=("`xhighdim_o'"~="")
	*** misc IV options. Change flag to 1/0. Ignore if no endog.
	local idstatsflag	=("`idstats'"~="" & "`dendog_o'"~="")
	local firstflag		=("`first'"~="" & "`dendog_o'"~="")
	local sscsetflag	=("`sscset'"~="" & `: word count `dendog_o''==1)
	local ssgridmatflag	=("`ssgridmat'"~="")
	local ssomitgridflag=("`ssomitgrid'"~="")
	local rlassoflag	=("`rlasso'`rlasso0'"~="")
	if `rlassoflag' {
		if "`rlasso0'"=="" {
			local rlassoprefix	_`cmdname'
		}
		else {
			local rlassoprefix	`rlasso0'
		}
	}
	*

	*** Create blank variables for orthogonalized Y and Xs and optimal IVs
	// All below have already had X ctrls partialled out
	// rho_y		depvar, residual after high-dim Xs partialled out
	// rho_d		low-dim exog X, residual after high-dim Xs partialled out
	// rho_e		low-dim endog X, residual after high-dim Xs partialled out
	// rho_z		notpen IV, residual after high-dim Xs partialled out
	// iv_e			optimal IV for low-dim endog X
	// dhat			fitted value of endog d on hi-dim Xs and all Zs
	// Notes:
	// r1			= dendog - dhat
	// dhathat		= fitted value of dhat on hi-dim Xs
	// r2			= dhat - dhathat = iv_e = vhat in CHS paper
	// rho_e		= dendog - dhathat
	//				= dendog - (dhat - r2) = dendog - dhat + r2
	//				= r1 + r2
	
	tempvar rho_y_l
	tempvar rho_y_pl
	qui gen double	`rho_y_l' = .
	qui gen double	`rho_y_pl' = .
//	local rho_y_o	rho_`varY_o'						//  prefix "rho_"
	local rho_y_o	`varY_o'							//  no prefix
	char `rho_y_l'[vname]  `varY_o'
	char `rho_y_pl'[vname] `varY_o'
	// low-dim (model) exogenous regressors
	local numvars	: word count `dexog_o'
	forvalues i=1/`numvars' {							//  create only if #vars > 0
		local var_o		: word `i' of `dexog_o'
		tempvar tvar1 tvar2
		qui gen double `tvar1' = .
		qui gen double `tvar2' = .
		local rho_d_l	`rho_d_l' `tvar1'
		local rho_d_pl	`rho_d_pl' `tvar2'
//		local rho_d_o	`rho_d_o' rho_`var_o'			//  prefix "rho_"
		local rho_d_o	`rho_d_o' `var_o'				//  no prefix
		char `tvar1'[vname] `var_o'
		char `tvar2'[vname] `var_o'
	}
	// low-dim (model) IVs
	local numvars	: word count `znotpen_o'
	forvalues i=1/`numvars' {							//  create only if #vars > 0
		local var_o		: word `i' of `znotpen_o'
		tempvar tvar1 tvar2
		qui gen double `tvar1' = .
		qui gen double `tvar2' = .
		local rho_z_l	`rho_z_l' `tvar1'
		local rho_z_pl	`rho_z_pl' `tvar2'
//		local rho_z_o	`rho_z_o' rho_`var_o'			//  prefix "rho_"
		local rho_z_o	`rho_z_o' `var_o'				//  no prefix
		char `tvar1'[vname] `var_o'
		char `tvar2'[vname] `var_o'
	}
	// low-dim (model) endogenous regressors
	local numvars	: word count `dendog_o'
	forvalues i=1/`numvars' {							//  create only if #vars > 0
		local var_o		: word `i' of `dendog_o'
		tempvar tvar1 tvar2
		qui gen double `tvar1' = .
		qui gen double `tvar2' = .
		local rho_e_l	`rho_e_l' `tvar1'
		local rho_e_pl	`rho_e_pl' `tvar2'
//		local rho_e_o	`rho_e_o' rho_`var_o'			//  prefix "rho_"
		local rho_e_o	`rho_e_o' `var_o'				//  no prefix
		char `tvar1'[vname] `var_o'
		char `tvar2'[vname] `var_o'
	}
	// fitted values for endogenous regressors d
	local numvars	: word count `dendog_o'
	forvalues i=1/`numvars' {							//  create only if #vars > 0
		local var_o		: word `i' of `dendog_o'
		tempvar tvar1 tvar2
		qui gen double `tvar1' = .
		qui gen double `tvar2' = .
		local dhat_l	`dhat_l' `tvar1'
		local dhat_pl	`dhat_pl' `tvar2'
//		local dhat_o	`dhat_o' hat_`var_o'			//  prefix "hat_"
		local dhat_o	`dhat_o' `var_o'				//  no prefix
		char `tvar1'[vname] `var_o'
		char `tvar2'[vname] `var_o'
	}
	// optimal IVs
	local numvars	: word count `dendog_o'
	local optivflag	=`numvars'							//  awkward but useful
	forvalues i=1/`numvars' {							//  create only if #vars > 0
		local var_o		: word `i' of `dendog_o'
		tempvar tvar1 tvar2
		qui gen double `tvar1' = .
		qui gen double `tvar2' = .
		local iv_e_l	`iv_e_l' `tvar1'
		local iv_e_pl	`iv_e_pl' `tvar2'
		local iv_e_o	`iv_e_o' iv_`var_o'				//  prefix "iv_"
//		local iv_e_o	`iv_e_o' `var_o'				//  no prefix
		char `tvar1'[vname] `var_o'
		char `tvar2'[vname] `var_o'
	}
	// create dictionary for new variables
	// (don't add to dict_o/dict_t since lasso and post-lasso entries are different)
	local newvars_o		`rho_y_o'  `rho_d_o'  `rho_z_o'  `rho_e_o'  `dhat_o'  `iv_e_o'
	local newvars_l		`rho_y_l'  `rho_d_l'  `rho_z_l'  `rho_e_l'  `dhat_l'  `iv_e_l'
	local newvars_pl	`rho_y_pl' `rho_d_pl' `rho_z_pl' `rho_e_pl' `dhat_pl' `iv_e_pl'
	*

	******************* Partialling out ***********************************************
	
	//  If FE:    partial-out FEs from temp variables, then preserve,
	//            then partial-out low-dim ctrls from temp variables
	//            restore will restore all temp vars with only FEs partialled-out
	//  If no FE: leave original variables unchanged.
	//            partial-out low-dim ctrls from temp variables.
	//            if no FE/low-dim ctrls, no transform needed

	// Initialize count of FEs to default of 0.
	local N_g			= 0
	if `feflag' {
		di as text "Fixed effects transformation..."
		// transform everything
		_fe `allvars_f',									/// _f list is _o list but with FV and TS vars replaced with temps
						touse(`touse')						///
						tvarlist(`allvars_t')				/// overwrite/initialize these
						`noftools'							/// provisional
						fe(`ivar')							//  
		local N_g	=r(N_g)									//  N_g will be 0 if no FEs
		local noftools `r(noftools)'						//  either not installed or user option
		local dofopt "dofminus(`N_g')"						//  to pass to post-lasso estimation; empty if no FEs
		// And then partial out any additional X vars	
		if `xpartialflag' {
			preserve										//  preserve the original values of tempvars
			_partial `allbutxpartial_t',					/// transform the tempvars
							touse(`touse')					///
							partial(`xpartial_t')			/// partial out X vars only; xpartial_t vars are FE-transformed
							cons(0)							//  FE => nocons
		}
	}
	else if `xpartialflag' {								//  Just partial out
		di as text "Partialling out unpenalized controls..."
		// transform everything except xpartial vars
		_partial `allbutxpartial_f',						/// _f list is _o list but with FV and TS vars replaced with temps
						touse(`touse')						///
						partial(`xpartial_f')				/// partial out X vars only; xpartial_f have FV and TS vars replaced with temps
						tvarlist(`allbutxpartial_t')		/// overwrite/initialize these
						cons(`cons')						//  cons possible
	}

	// Tell lassoshooting if cons has been partialled out or there isn't one in the first place
	if `feflag' | `xpartialflag' | (~`cons') {
		local lscons	0
	}
	else {
		local lscons	1
	}

	************* partialling out END ***********************************************

	************* Selection of vars ************************************************	
	// SelectControls selects HD vars and also creates either immunized LD vars (resids)
	// or creates optimal IVs (xb option => fitted values).
	
	local startcol =11							//  used for display of messages

	*** Step 1: orthogonalize depvar w.r.t. HD Xs
	// All estimations but only if there are HD X vars; if no HD Xs, selected=empty and rho=y.
	// Corresponds to CHS paper steps of obtaining theta^ (coefs on X) and rho^_y (resids).
	// Also used in PDS method.
	local msg	"1.  (PDS/CHS) Selecting HD controls for dep var"
	SelectControls,								///
					resid						/// partial out hi-dim Xs from depvar
					touse(`touse')				///
					modelvars_o(`varY_o')		///
					modelvars_t(`varY_t')		///
					genvars_l(`rho_y_l')		/// updates to resid or to y if nothing selected
					genvars_pl(`rho_y_pl')		/// updates to resid or to y if nothing selected
					highdim_o(`xhighdim_o')		///
					highdim_t(`xhighdim_t')		///
					dict_o(`dict_o')			/// dictionary
					dict_t(`dict_t')			/// dictionary
					aset_o(`xaset_o')			///
					aset_t(`xaset_t')			///
					notpen_o(`xnotpen_o')		///
					notpen_t(`xnotpen_t')		///
					loptions(`loptions')		///
					msg("`msg'")				///
					startcol(`startcol')		///
					debugflag(`debugflag')		///
					lscons(`lscons')			///
					`qui'
	local xselected1_o	`r(allselected0_o)'
	if `rlassoflag' & "`xhighdim_o'`xnotpen_o'"~="" {
		capture est store `rlassoprefix'_step1, title("lasso step 1")
		if _rc > 0 {
			di as text "Warning: unable to store lasso estimation"
		}
		else {
			local rlassolist `rlassolist' `rlassoprefix'_step1
		}
	}
	*

	*** Step 2: orthogonalize exog regressors w.r.t. HD Xs
	// Only if there are exog regressors d; if no HD Xs, selected=empty and rho=d.
	if "`dexog_o'" ~= "" {
		local msg	"2.  (PDS/CHS) Selecting HD controls for exog regressor"
		SelectControls,								///
						resid						/// partial out hi-dim Xs from exog regressors d
						touse(`touse')				///
						modelvars_o(`dexog_o')		///
						modelvars_t(`dexog_t')		///
						genvars_l(`rho_d_l')		/// updates to resid or to d if nothing selected
						genvars_pl(`rho_d_pl')		/// updates to resid or to d if nothing selected
						highdim_o(`xhighdim_o')		///
						highdim_t(`xhighdim_t')		///
						dict_o(`dict_o')			/// dictionary
						dict_t(`dict_t')			/// dictionary
						aset_o(`xaset_o')			///
						aset_t(`xaset_t')			///
						notpen_o(`xnotpen_o')		///
						notpen_t(`xnotpen_t')		///
						loptions(`loptions')		///
						startcol(`startcol')		///
						msg("`msg'")				///
						debugflag(`debugflag')		///
						lscons(`lscons')			///
						`qui'
		local xselected2_o	`r(allselected0_o)'
		if `rlassoflag' & "`xhighdim_o'`xnotpen_o'"~="" {
			capture est store `rlassoprefix'_step2, title("lasso step 2")
			if _rc > 0 {
				di as text "Warning: unable to store lasso estimation"
			}
			else {
				local rlassolist `rlassolist' `rlassoprefix'_step2
			}
		}
	}
	*

	*** Step 3: (PDS only) endog Xs need HD Xs
	// PDS only; selection only, no variables generated
	if "`dendog_o'" ~= "" {
		local msg	"3.  (PDS) Selecting HD controls for endog regressor"
		SelectControls,								/// variable selection only, no vars generated
						/* resid */					/// irrelevant since nothing generated
						touse(`touse')				///
						modelvars_o(`dendog_o')		///
						modelvars_t(`dendog_t')		///
						highdim_o(`xhighdim_o')		///
						highdim_t(`xhighdim_t')		///
						dict_o(`dict_o')			/// dictionary
						dict_t(`dict_t')			/// dictionary
						aset_o(`xaset_o')			///
						aset_t(`xaset_t')			///
						notpen_o(`xnotpen_o')		///
						notpen_t(`xnotpen_t')		///
						loptions(`loptions')		///
						msg("`msg'")				///
						startcol(`startcol')		///
						debugflag(`debugflag')		///
						lscons(`lscons')			///
						`qui'
		local xselected3_o	`r(allselected0_o)'
		if `rlassoflag' & "`xhighdim_o'`xnotpen_o'"~="" {
			capture est store `rlassoprefix'_step3, title("lasso step 3")
			if _rc > 0 {
				di as text "Warning: unable to store lasso estimation"
			}
			else {
				local rlassolist `rlassolist' `rlassoprefix'_step3
			}
		}
	}
	*

	*** Step 4: (PDS only) If no HD Zs, LD Zs need HD Xs
	// PDS only; selection only, no variables generated
	// If no HD Zs, lasso of LD Zs on HD Xs
	// NB: if HD Zs exist, do with opt IV step - lasso of endog X on HD Xs, LD and HD Zs with LD Zs unpenalized.
	if "`dendog_o'" ~= "" & (~`hdzflag') {						//  no HD IVs
		local msg	"4.  (PDS) Selecting HD controls for IV"
		SelectControls,											/// variable selection only, no vars generated
						/* resid */								/// irrelevant since nothing generated
						touse(`touse')							///
						modelvars_o(`znotpen_o' `zpartial_o')	/// notpen Zs are LD Zs and need HD controls
						modelvars_t(`znotpen_t' `zpartial_t')	///
						highdim_o(`xhighdim_o')					///
						highdim_t(`xhighdim_t')					///
						dict_o(`dict_o')						/// dictionary
						dict_t(`dict_t')						/// dictionary
						aset_o(`xaset_o')						///
						aset_t(`xaset_t')						///
						notpen_o(`xnotpen_o')					///
						notpen_t(`xnotpen_t')					///
						loptions(`loptions')					///
						msg("`msg'")							///
						startcol(`startcol')					///
						debugflag(`debugflag')					///
						lscons(`lscons')						///
						`qui'
		local xselected4_o	`r(allselected0_o)'
		if `rlassoflag' & "`xhighdim_o'`xnotpen_o'"~="" {
			capture est store `rlassoprefix'_step4, title("lasso step 4")
			if _rc > 0 {
				di as text "Warning: unable to store lasso estimation"
			}
			else {
				local rlassolist `rlassolist' `rlassoprefix'_step4
			}
		}
	}
	*

	*** Step 5: endog d on hi-dim Xs and all Zs
	// only if there are endogenous regressors d; get fitted values
	// if no IVs selected, fitted values = 0 and model is unidentified
	// if notpen Zs and HD Zs:
	//		selects HD Xs and all IVs
	//		generates fitted values `dhat_l' & `dhat_pl' for CHS opt IV construction
	//		selected vars used for PDS algo
	// if notpen Zs and no HD Zs:
	//		selects HD Xs and notpen Zs
	//		generates fitted values `dhat_l' & `dhat_pl' for opt IV construction
	//		no need for selected vars
	// note that PDS algo uses EITHER step 4 OR step 5 selected vars
	if "`dendog_o'" ~= "" {
		if `hdzflag' {
			local msg	"5.  (PDS/CHS) Selecting HD controls/IVs for endog regressor"
		}
		else {
			local msg	"5.  (CHS) Selecting HD controls and IVs for endog regressor"
		}
		SelectControls,											///
						xb										/// get fitted values
						touse(`touse')							///
						modelvars_o(`dendog_o')					///
						modelvars_t(`dendog_t')					///
						genvars_l(`dhat_l')						/// updates to fitted value or 0 if nothing selected
						genvars_pl(`dhat_pl')					/// updates to fitted value or 0 if nothing selected
						highdim_o(`xhighdim_o' `zhighdim_o')	///
						highdim_t(`xhighdim_t' `zhighdim_t')	///
						dict_o(`dict_o')						/// dictionary
						dict_t(`dict_t')						/// dictionary
						aset_o(`xaset_o' `zaset_o')				///
						aset_t(`xaset_t' `zaset_t')				///
						notpen_o(`xnotpen_o' `znotpen_o')		///
						notpen_t(`xnotpen_t' `znotpen_t')		///
						partial_o(`zpartial_o')					///
						partial_t(`zpartial_t')					///
						loptions(`loptions')					///
						msg("`msg'")							///
						startcol(`startcol')					///
						debugflag(`debugflag')					///
						lscons(`lscons')						///
						`qui'
		local zselected5_o	`r(allselected0_o)'
		if `rlassoflag' & "`xhighdim_o'`zhighdim_o'`xnotpen_o'`znotpen_o'`zpartial_o'"~="" {
			capture est store `rlassoprefix'_step5, title("lasso step 5")
			if _rc > 0 {
				di as text "Warning: unable to store lasso estimation"
			}
			else {
				local rlassolist `rlassolist' `rlassoprefix'_step5
			}
		}
	}
	*

	*** Step 6: (CHS only) create optimal instruments based on fitted values of endog Xs
	// only if there are endogenous regressors; get residuals = optimal IVs;
	// iv_e called vhat in CHS paper; separately for lasso and post-lasso
	if "`dendog_o'" ~= "" {
		local msg	"6a. (CHS) Selecting `method' HD controls and creating optimal IV for endog regressor"
		SelectControls,										///
						resid								/// partial out hi-dim Xs from dhat=fitted values of endog
						touse(`touse')						///
						modelvars_o(`dhat_o')				///
						modelvars_t(`dhat_l')				/// selecting HD controls for lasso-fitted values of endog
						genvars_l(`iv_e_l')					/// updates to resid or to x if nothing selected
						highdim_o(`xhighdim_o')				///
						highdim_t(`xhighdim_t')				///
						dict_o(`dict_o' `newvars_o')		/// need newvars in dictionary
						dict_t(`dict_t' `newvars_l')		/// must be the LASSO newvars
						/* aset_o(`xaset_o') */				/// aset irrelevant for lasso
						/* aset_t(`xaset_t') */				/// aset irrelevant for lasso
						notpen_o(`xnotpen_o')				///
						notpen_t(`xnotpen_t')				///
						loptions(`loptions')				///
						msg("`msg'")						///
						startcol(`startcol')				///
						debugflag(`debugflag')				///
						lscons(`lscons')					///
						`qui'
		local xselected6a_o		`r(allselected0_o)'
		if `rlassoflag' & "`xhighdim_o'`xnotpen_o'"~="" {
			capture est store `rlassoprefix'_step6a, title("lasso step 6a")
			if _rc > 0 {
				di as text "Warning: unable to store lasso estimation"
			}
			else {
				local rlassolist `rlassolist' `rlassoprefix'_step6a
			}
		}
	}
	if "`dendog_o'" ~= "" {
		local msg	"6b. (CHS) Selecting post-`method' HD controls and creating optimal IV for endog regressor"
		SelectControls,										///
						resid								/// partial out hi-dim Xs from fitted values of endog
						touse(`touse')						///
						modelvars_o(`dhat_o')				///
						modelvars_t(`dhat_pl')				/// selecting HD controls for post-lasso-fitted values of endog
						genvars_pl(`iv_e_pl')				/// updates to resid or to x if nothing selected
						highdim_o(`xhighdim_o')				///
						highdim_t(`xhighdim_t')				///
						dict_o(`dict_o' `newvars_o')		/// need newvars in dictionary
						dict_t(`dict_t' `newvars_pl')		/// must be the POST-LASSO newvars
						aset_o(`xaset_o')					///
						aset_t(`xaset_t')					///
						notpen_o(`xnotpen_o')				///
						notpen_t(`xnotpen_t')				///
						loptions(`loptions')				///
						msg("`msg'")						///
						startcol(`startcol')				///
						debugflag(`debugflag')				///
						lscons(`lscons')					///
						`qui'
		local xselected6b_o		`r(allselected0_o)'
		if `rlassoflag' & "`xhighdim_o'`xnotpen_o'"~="" {
			capture est store `rlassoprefix'_step6b, title("lasso step 6b")
			if _rc > 0 {
				di as text "Warning: unable to store lasso estimation"
			}
			else {
				local rlassolist `rlassolist' `rlassoprefix'_step6b
			}
		}
	}
	*

	*** Step 7: Create orthogonalized endogenous regressors rho_e
	// only if there are endogenous regressors dendog (so we created dhat from lasso of dendog on Xs and Zs)
	// = dendog  -  dhathat (predicted value from lasso of dhat on Xs)
	// = dendog  -  (dhat - resid from dhat on Xs)
	// = dendog  -  (dhat - optimal IV)
	if "`dendog_o'" ~= "" {
		local numvars	: word count `dendog_o'
		forvalues i=1/`numvars' {
			local dvar_o	: word `i' of `dendog_o'
			local dvar_t	: word `i' of `dendog_t'
			local xbvar_l	: word `i' of `dhat_l'
			local xbvar_pl	: word `i' of `dhat_pl'
			local ivvar_l	: word `i' of `iv_e_l'
			local ivvar_pl	: word `i' of `iv_e_pl'
			local rhovar_l	: word `i' of `rho_e_l'
			local rhovar_pl	: word `i' of `rho_e_pl'
			
			di as text "7.  (CHS) Creating orthogonalized endogenous regressor `dvar_o'..."
			qui replace `rhovar_l'	= `dvar_t' - ( `xbvar_l' -  `ivvar_l'  )
			qui replace `rhovar_pl'	= `dvar_t' - ( `xbvar_pl' - `ivvar_pl' )
		}
	}
	*
	*********************** sup-score stat **********************************
	// IV estimation only
	// code uses tempvars so partialling or FE has been done
	// execute only if idstats or supscore CI requested
	if (`idstatsflag' | `sscsetflag' | `ssgridmatflag') & "`dendog_o'"~="" & "`dexog_o'"=="" & "`znotpen_o'`zpartial_o'`zaset_o'"=="" {
		di "8.  (PDS/CHS) Performing sup-score test(s)"
		doSupScore,												/// 
						touse(`touse')							///
						varY_o(`varY_o')						///
						varY_t(`varY_t')						///
						rho_y_o(`rho_y_o')						/// orthogonalized dep var
						rho_y_l(`rho_y_l')						/// orthogonalized dep var
						rho_y_pl(`rho_y_pl')					/// orthogonalized dep var
						dendog_o(`dendog_o')					///
						dendog_t(`dendog_t')					///
						xhighdim_o(`xhighdim_o')				///
						xhighdim_t(`xhighdim_t')				///
						zhighdim_o(`zhighdim_o')				///
						zhighdim_t(`zhighdim_t')				///
						xnotpen_o(`xnotpen_o')					///
						xnotpen_t(`xnotpen_t')					///
						xaset_o(`xaset_o')						///
						xaset_t(`xaset_t')						///
						dict_o(`dict_o')						/// dictionary
						dict_t(`dict_t')						/// dictionary
						loptions(`loptions')					///
						ssgamma(`ssgamma')						///
						sscsetflag(`sscsetflag')				///
						ssmethod(`ssmethod')					///
						gridmin(`ssgridmin')					///
						gridmax(`ssgridmax')					///
						gridpoints(`ssgridpoints')				///
						gridmat(`ssgridmat')					///
						startcol(`startcol')					///
						debugflag(`debugflag')					///
						lscons(`lscons')						///
						dofopt(`dofopt')						/// needed for OLS
						ivoptions(`ivoptions')					/// needed for OLS
						cons(`cons')							/// needed for OLS
						`qui'
		local ss_null			`r(ss_null)'
		local ss_null_l			`r(ss_null_l)'
		local ss_null_pl		`r(ss_null_pl)'
		local ss_cset			`r(ss_cset)'
		local ss_cset_l			`r(ss_cset_l)'
		local ss_cset_pl		`r(ss_cset_pl)'
		local ss_gridmin		=r(gridmin)
		local ss_gridmax		=r(gridmax)
		local ss_gridpoints		=`ssgridpoints'					//  rename here
		if `ssgridmatflag' {
			tempname ss_gridmat
			mat `ss_gridmat'	=r(gridresults)
		}
		if `sscsetflag' {
			tempname ss_citable
			mat `ss_citable'	=r(citable)
		}
		local ss_method			`r(ss_method)'
	}
	*

	*************************************************************************
	// newline before main results
	di
	*************************************************************************

	*** Consolidate list of PDS HD controls
	// "selected" locals from above include partial (Zs only), notpen and cons
	// list may include Zs so remove them
	// necessary only if there are controls
	local xselected_o	`xselected1_o' `xselected2_o' `xselected3_o' `xselected4_o' `zselected5_o'
	local xselected_o	`xselected_o' `xnotpen_o'			//  re-add notpen in case FV base vars in notpen list (=> not selected)
	local xselected_o	: list xselected_o - consname		//  remove "_cons" if present
	local xselected_o	: list uniq xselected_o				//  remove duplicates
	local xselected_o	: list xselected_o - zhighdim_o		//  remove any IVs from list of controls
	local xselected_o	: list xselected_o - znotpen_o
	local xselected_o	: list xselected_o - zpartial_o
	// put in original order (as in dictionary macro)
	local unselected_o	: list dict_o - xselected_o
	local xselected_o	: list dict_o - unselected_o
	matchnames "`xselected_o'" "`dict_o'" "`dict_t'"
	local xselected_t	`r(names)'
	// xselected now has all selected X variables including notpen but not cons
	*
	*** Consolidate list of PDS/CHS IVs
	// "selected" locals from above include partial (Zs only), notpen and cons
	// list may include Zs so remove them
	local zselected_o	`zselected5_o'
	local zselected_o	`zselected_o' `znotpen_o'			//  re-add notpen in case FV base vars in notpen list (=> not selected)
	local zselected_o	: list zselected_o - consname		//  remove "_cons" if present
	local zselected_o	: list uniq zselected_o				//  remove duplicates
	local zselected_o	: list zselected_o - xhighdim_o		//  remove HD controls
	local zselected_o	: list zselected_o - xnotpen_o
	// put in original order
	local unselected_o	: list dict_o - zselected_o
	local zselected_o	: list dict_o - unselected_o
	matchnames "`zselected_o'" "`dict_o'" "`dict_t'"
	local zselected_t	`r(names)'
	// zselected now has all selected Z variables including notpen, partial but not cons
	*
	
	*** Consolidate list of CHS HD controls
	local xselected_chs_l_o		`xselected1_o' `xselected2_o' `zselected5_o' `xselected6a_o'
	local xselected_chs_l_o		: list xselected_chs_l_o - zhighdim_o		//  remove HD IVs
	local xselected_chs_l_o		: list xselected_chs_l_o - znotpen_o		//  remove notpen IVs
	local xselected_chs_l_o		: list xselected_chs_l_o - zpartial_o		//  remove partialled-out IVs
	local xselected_chs_l_o		: list uniq xselected_chs_l_o				//  remove duplicates
	local xselected_chs_pl_o	`xselected1_o' `xselected2_o' `zselected5_o' `xselected6b_o' `xaset_o'
	local xselected_chs_pl_o	: list xselected_chs_pl_o - zhighdim_o		//  remove HD IVs
	local xselected_chs_pl_o	: list xselected_chs_pl_o - znotpen_o		//  remove notpen IVs
	local xselected_chs_pl_o	: list xselected_chs_pl_o - zpartial_o		//  remove partialled-out IVs
	local xselected_chs_pl_o	: list uniq xselected_chs_pl_o				//  remove duplicates
	// put in original order (_t versions not needed)
	local unselected_o			: list dict_o - xselected_chs_l_o
	local xselected_chs_l_o		: list dict_o - unselected_o
	local unselected_o			: list dict_o - xselected_chs_pl_o
	local xselected_chs_pl_o	: list dict_o - unselected_o
	*

	*** Variable counts
	local xhighdim_ct			: word count `xhighdim_o'
	local zhighdim_ct			: word count `zhighdim_o'
	local dexog_ct				: word count `dexog_o'
	local dendog_ct				: word count `dendog_o'
	local xselected_ct			: word count `xselected_o'
	local xselected_chs_l_ct	: word count `xselected_chs_l_o'
	local xselected_chs_pl_ct	: word count `xselected_chs_pl_o'
	local xaset_ct				: word count `xaset_o'
	local xnotpen_ct			: word count `xnotpen_o'
	local xpartial_ct			: word count `xpartial_o'
	local zselected_ct			: word count `zselected_o'
	local zaset_ct				: word count `zaset_o'
	local znotpen_ct			: word count `znotpen_o'
	local zpartial_ct			: word count `zpartial_o'
	if `cons' {														//  adjust xnotpen_ct and xpartial_ct; won't enter if FE
		if `xpartialflag' {											//  any partialling-out => _cons partialled out too
			local xpartial_ct	= `xpartial_ct' + 1					//  partialling-out => add 1 to partial
		}
		else {
			local xnotpen_ct	= `xnotpen_ct' + 1					//  no partialling-out  => add 1 to notpen
		}
	}
	*

	****************************************************************************
	************* Estimation results********************************************
	****************************************************************************

	************* First-stage estimations **************************************	

	*** First-stage regression reports only PDS-type first-stage
	// data may still be transformed by FE/partialling
	// regressors are xselected (including notpen, doesn't include X partials)
	// and zselected (includes notpen, partial)
	if `firstflag' {
		foreach var_t of varlist `dendog_t' {
			matchnames "`var_t'" "`dict_t'" "`dict_o'"
			local var_o			`r(names)'
			local regressors	`dexog_t' `xselected_t' `zselected_t'
			doregress `var_t' if `touse',								///
				regressors(`regressors')								///
				instruments(`regressors')								///
				dict_o(`dict_o')										///
				dict_t(`dict_t')										///
				cons(`lscons')											/// nocons if partialled out or FE
				`dofopt'												///
				`ivoptions'												///
				debugflag(`debugflag')
			ereturn local firstvar `var_o'								//  also serves as flag for replay
			capture est store _ivlasso_`var_o', title("First-stage equation for `var_o'")
			if _rc > 0 {
				di
				di as text "Warning: unable to store first-stage eqn for `var_o'"
				di
			}
			else {
				local firstlist `firstlist' _ivlasso_`var_o'
			}
		}
	}
	*
	
	*********************** Identification stats **********************************
	// code uses tempvars so partialling or FE has been done
	// estimation uses xselected (including notpen, doesn't include X partials)
	// and zselected (includes notpen, partial)
	// xpartial used in inexog count
	if `idstatsflag' {
		// collinearities possible - affects dofs and ranktest, so remove
		if `lscons' {
			local rmcons noconstant
		}
		local exexog		`zselected_t'
		local inexog		`dexog_t' `xselected_t'
		// cap since would crash with empty varlist (nothing selected)
		cap _rmcollright `inexog' `exexog' if `touse', `rmcons'
		local dropped		`r(dropped)'
		local exexog		: list exexog - dropped
		local inexog		: list inexog - dropped
		local inexog_ct		: word count `inexog'
		fvstrip `xpartial_o' if `touse', dropomit						//  partial-out list may include omitteds etc.
		local xpartial_ct	: word count `r(varlist)'
		local inexog_ct		= `inexog_ct' + `xpartial_ct' + `N_g' + `cons'
		doIDstats if `touse',											///
			endog(`dendog_t')											///
			endog_l(`rho_y_l')											///
			endog_pl(`rho_y_pl')										///
			exexog(`exexog')											/// after collinearities removed
			exexog_l(`iv_e_l')											///
			exexog_pl(`iv_e_pl')										///
			inexog(`inexog')											///	after collinearities removed
			inexog_ct(`inexog_ct')										///
			cons(`lscons')												/// no cons if FE or partialled out or no cons in the first place
			`robust'													///
			cluster(`cluster')
		local weakid			=r(weakid)
		if `optivflag' {
			local weakid_l		=r(weakid_l)
			local weakid_pl		=r(weakid_pl)
		}
		if "`robust'`cluster'" != "" {		
			local weakidr		=r(weakidr)
			if `optivflag' {
				local weakidr_l	=r(weakidr_l)
				local weakidr_pl	=r(weakidr_pl)
			}
		}		
	}		// End of identification stats

	************* STRUCTURAL EQUATION RESULTS *********************************************
	
	tempname beta_pds beta_lasso beta_plasso V_pds V_lasso V_plasso
	
	************* CHS estimation using orthogonalized vars ************************

	// CHS lasso-based
	doregress `rho_y_l' if `touse',					///
		regressors(`rho_e_l' `rho_d_l')				///
		instruments(`iv_e_l' `rho_d_l')				///
		dict_o(`newvars_o')							///
		dict_t(`newvars_l')							///
		cons(0)										///
		`dofopt'									///
		`ivoptions'									///
		debugflag(`debugflag')
	mat `beta_lasso'=e(b)
	mat `V_lasso'	=e(V)
	// set here since doregress with CHS more graceful if model unidentified
	local N_clust	=e(N_clust)						//  doregress returns 0 if no clustering
	// CHS post-lasso-based
	doregress `rho_y_pl' if `touse',				///
		regressors(`rho_e_pl' `rho_d_pl')			///
		instruments(`iv_e_pl' `rho_d_pl')			///
		dict_o(`newvars_o')							///
		dict_t(`newvars_pl')						///
		cons(0)										///
		`dofopt'									///
		`ivoptions'									///
		debugflag(`debugflag')
	mat `beta_plasso'=e(b)
	mat `V_plasso'	=e(V)
	*

	****************** PDS estimation with full regressor set *****************

	if `feflag' & `xpartialflag' {					//  FE case and there are partialled-out notpen vars
		restore										//  Restores dataset with tempvars after FE transform but before notpen partialled out
	}
	// xselected has all selected X except for xpartial, so need to add that
	// zselected has all selected Z including notpen and partial, so nothing else needed
	if `feflag' {									//  if FE, all vars have been FE-transformed so use _t vars
		local varY_pds		`varY_t'
		local regressors	`dendog_t' `dexog_t' `xselected_t' `xpartial_t'
		local instruments	`dexog_t' `xselected_t' `xpartial_t' `zselected_t'
		local pds_dict_o	`dict_o'
		local pds_dict_t	`dict_t'
	}
	else {											//  if not FE, use _f list = original vars with FVs replaced by temps
		local varY_pds		`varY_o'
		local regressors	`dendog_o' `dexog_o' `xselected_o' `xpartial_o'
		local instruments	`dexog_o' `xselected_o' `xpartial_o' `zselected_o'
		foreach vlist in varY_pds regressors instruments {
			matchnames "``vlist''" "`dict_o'" "`dict_f'"
			local `vlist'		`r(names)'		//  corresponding tempnames of vars
		}
		local pds_dict_o	`dict_o'
		local pds_dict_t	`dict_f'
	}

	doregress `varY_pds' if `touse',				///
		regressors(`regressors')					///
		instruments(`instruments')					///
		dict_o(`pds_dict_o')						///
		dict_t(`pds_dict_t')						///
		cons(`cons')								/// if fe, nocons is automatic
		`dofopt'									///
		`ivoptions'
	mat `beta_pds'	=e(b)
	mat `V_pds'		=e(V)

	************************************************************************************
	****************** Post and save results *******************************************
	************************************************************************************

	// Need this unfortunately because ereturn post means b & V cease to exist
	tempname b V
	if "`post'"=="pds" {
		mat `b'		=`beta_pds'
		mat `V'		=`V_pds'
	}
	else if "`post'"=="lasso" {
		mat `b'		=`beta_lasso'
		mat `V'		=`V_lasso'
	}
	else if "`post'"=="plasso" {
		mat `b'		=`beta_plasso'
		mat `V'		=`V_plasso'
	}
	ereturn post `b' `V', obs(`N') depname(`varY_o') esample(`touse')
		
	ereturn mat beta_pds	=`beta_pds'
	ereturn mat V_pds		=`V_pds'
	ereturn mat beta_lasso	=`beta_lasso'
	ereturn mat V_lasso		=`V_lasso'
	ereturn mat beta_plasso	=`beta_plasso'
	ereturn mat V_plasso	=`V_plasso'
	
	ereturn local noftools	`noftools'
	ereturn local predict	ivlasso_p
	ereturn local post		`post'				//  pds, lasso, plasso
	ereturn local method	`method'			//  lasso, sqrt-lasso
	ereturn local cmd 		ivlasso

	if "`robust'`cluster'"~="" {
		ereturn local vce		robust
		ereturn local vcetype	Robust
	}
	if `N_clust' {
		ereturn local vce		cluster
	}
	ereturn local clustvar			`cluster'
	ereturn scalar N_clust			=`N_clust'

	if `xpartialflag' & `cons' {
		ereturn local xpartial		`xpartial_o' _cons
		ereturn local xnotpen		`xnotpen_o'
	}
	else if `xpartialflag' {
		ereturn local xpartial		`xpartial_o'
		ereturn local xnotpen		`xnotpen_o'
	}
	else if `cons' {
		ereturn local xnotpen		`xnotpen_o' _cons
	}
	else {
		ereturn local xnotpen		`xnotpen_o'
	}
	
	ereturn local znotpen			`znotpen_o'
	ereturn local zpartial			`zpartial_o'
	ereturn local zaset				`zaset_o'
	ereturn local zselected			`zselected_o'
	ereturn local zhighdim			`zhighdim_o'

	ereturn local xaset				`xaset_o'
	ereturn local xselected_chs_l	`xselected_chs_l_o'
	ereturn local xselected_chs_pl	`xselected_chs_pl_o'
	ereturn local xselected			`xselected_o'
	ereturn local xselected_y		`xselected1_o'
	ereturn local xselected_dex		`xselected2_o'
	ereturn local xselected_dend	`xselected3_o'
	ereturn local xselected_z		`xselected4_o'
	ereturn local xzselected_dend	`zselected5_o'
	ereturn local xselected_z_l		`xselected6a_o'
	ereturn local xselected_z_pl	`xselected6b_o'
	ereturn local xhighdim			`xhighdim_o'

	ereturn local dendog			`dendog_o'	
	ereturn local dexog				`dexog_o'

	// FE estimation
	ereturn local ivar				`ivar'
	ereturn scalar fe				=`feflag'
	ereturn scalar N_g				=`N_g'

	ereturn local rlassolist		`rlassolist'
	ereturn local firstlist			`firstlist'
	
	// Currently all 3 estimates produced
	ereturn scalar pds					=1
	ereturn scalar chs_l				=1
	ereturn scalar chs_pl				=1
	
	ereturn scalar cons					=`cons'
	ereturn scalar dexog_ct				=`dexog_ct'
	ereturn scalar dendog_ct			=`dendog_ct'
	ereturn scalar xhighdim_ct			=`xhighdim_ct'
	ereturn scalar xselected_ct			=`xselected_ct'
	ereturn scalar xselected_chs_l_ct	=`xselected_chs_l_ct'
	ereturn scalar xselected_chs_pl_ct	=`xselected_chs_pl_ct'
	ereturn scalar xaset_ct				=`xaset_ct'
	ereturn scalar xnotpen_ct			=`xnotpen_ct'
	ereturn scalar xpartial_ct			=`xpartial_ct'
	ereturn scalar zhighdim_ct			=`zhighdim_ct'
	ereturn scalar zselected_ct			=`zselected_ct'
	ereturn scalar zaset_ct				=`zaset_ct'
	ereturn scalar znotpen_ct			=`znotpen_ct'
	ereturn scalar zpartial_ct			=`zpartial_ct'

	if `idstatsflag' {
		ereturn local ss_null		`ss_null'
		ereturn local ss_null_l		`ss_null_l'
		ereturn local ss_null_pl	`ss_null_pl'
		ereturn scalar ss_level		=(1-`ssgamma')*100
		ereturn scalar ss_gamma		=`ssgamma'
	}
	if "`ss_cset'`ss_cset_l'`ss_cset_pl'" ~= "" {
		ereturn local ss_cset			`ss_cset'
		ereturn local ss_cset_l			`ss_cset_l'
		ereturn local ss_cset_pl		`ss_cset_pl'
		ereturn scalar ss_gridmin		=`ss_gridmin'
		ereturn scalar ss_gridmax		=`ss_gridmax'
		ereturn scalar ss_gridpoints	=`ss_gridpoints'
	}
	if `ssgridmatflag' {
		ereturn matrix ss_gridmat		=`ss_gridmat'
		ereturn scalar ss_omitgrid		=`ssomitgridflag'
	}
	if `sscsetflag' {
		ereturn matrix ss_citable		=`ss_citable'
	}
	ereturn local ss_method				`ss_method'

	if `idstatsflag' {
		ereturn scalar weakid			=`weakid'
		ereturn scalar weakid_chs_l		=`weakid_l'
		ereturn scalar weakid_chs_pl	=`weakid_pl'
		if "`robust'`cluster'"~="" {
			ereturn scalar weakidr			=`weakidr'
			ereturn scalar weakidr_chs_l	=`weakidr_l'
			ereturn scalar weakidr_chs_pl	=`weakidr_pl'
		}
	}
	ereturn scalar idstats				=`idstatsflag'

	// do last so any overwriting doesn't matter
	ereturn local loptions			`loptions'
	// trick to parse lasso options
	local 0 ", `loptions'"
	syntax [, robust cluster(varlist) sqrt * ]
	if "`cluster'" ~= "" {
		ereturn local pload				cluster-lasso
	}
	else if "`robust'" ~= "" {
		ereturn local pload				heteroskedastic
	}
	else {
		ereturn local pload				homoskedastic
	}
	if "`sqrt'" ~= "" {
		ereturn local method			sqrt-lasso
	}
	*
	
end

*********************************** End main programm ********************

program define DisplayResults, eclass
	version 13
	
	if e(dendog_ct) {
		local estimator IV
	}
	else {
		local estimator OLS
	}
	
	di
	di as text "Estimation results:"
	
	// Display model specification.
	local startcol 40
	di
	di as text "Specification:"
	// how many spaces needed for #obs?
	local spaces	=strlen(strtrim("`: di %-12.0fc `e(N)''"))
	local fmt		%`spaces'.0fc
	di as text "Regularization method:" _col(`startcol') as res "`e(method)'"
	di as text "Penalty loadings:" _col(`startcol') as res "`e(pload)'"
	di as text "Number of observations:" _col(`startcol') as res `fmt' `e(N)'
	if e(N_clust) {
		di as text "Number of clusters:" _col(`startcol') as res `fmt' `e(N_clust)'
	}
	if e(N_g) {
		di as text "Number of fixed effects:" _col(`startcol') as res `fmt' `e(N_g)'
	}
	if e(dexog_ct) {
		di as text "Exogenous (`e(dexog_ct)'):" _c
		DispVars `e(dexog)', _col(`startcol')
	}
	if e(dendog_ct) {
		di as text "Endogenous (`e(dendog_ct)'):" _c
		DispVars `e(dendog)', _col(`startcol')
	}
	if e(xhighdim_ct) {
		di as text "High-dim controls (`e(xhighdim_ct)'):"_c
		DispVars `e(xhighdim)', _col(`startcol')
	}
	if e(dendog_ct)==0 {									//  no endogenous so HD ctrls same for all
		if e(xhighdim_ct) {
			di as text "Selected controls (`e(xselected_ct)'):"_c
			DispVars `e(xselected)', _col(`startcol')
		}
	}
	else {													//  endog => HD ctrls can differ
		if e(pds) & e(xhighdim_ct) {
			di as text "Selected controls, PDS (`e(xselected_ct)'):"_c
			DispVars `e(xselected)', _col(`startcol')
		}
		if e(chs_l) & e(xhighdim_ct) {
			di as text "Selected controls, CHS-L (`e(xselected_chs_l_ct)'):" _c
			DispVars `e(xselected_chs_l)', _col(`startcol')
		}
		if e(chs_pl) & e(xhighdim_ct) {
			di as text "Selected controls, CHS-PL (`e(xselected_chs_pl_ct)'):" _c
			DispVars `e(xselected_chs_pl)', _col(`startcol')
		}
	}
	if e(xaset_ct) {
		di as text "  of which, amelioration set (`e(xaset_ct)'):" _c
		DispVars `e(xaset)', _col(`startcol')
	}
	if e(xnotpen_ct) {
		di as text "Unpenalized controls (`e(xnotpen_ct)'):" _c
		DispVars `e(xnotpen)', _col(`startcol')
	}
	if e(xpartial_ct) {
		di as text "Partialled-out controls (`e(xpartial_ct)'):" _c
		DispVars `e(xpartial)', _col(`startcol')
	}
	if e(zhighdim_ct) {
		di as text "High-dim instruments (`e(zhighdim_ct)'):"_c
		DispVars `e(zhighdim)', _col(`startcol')
	}
	if e(pds) & e(zhighdim_ct) {
		di as text "Selected instruments (`e(zselected_ct)'):" _c
		DispVars `e(zselected)', _col(`startcol')
	}
	if e(zaset_ct) {
		di as text " of which, amelioration set (`e(zaset_ct)'):" _c
		DispVars `e(zaset)', _col(`startcol')
	}
	if e(znotpen_ct) {
		di as text "Unpenalized instruments (`e(znotpen_ct)'):" _c
		DispVars `e(znotpen)', _col(`startcol')
	}
	if e(zpartial_ct) {
		di as text "Partialled-out instruments (`e(zpartial_ct)'):" _c
		DispVars `e(zpartial)', _col(`startcol')
	}
	di

	// Display rlasso estimations.
	if "`e(rlassolist)'" ~= "" {
		di as text "lasso estimation(s):"
		foreach eqn in `e(rlassolist)' {
			di
			di as smcl "{stata estimates replay `eqn':`eqn'}" _c			//  avoid blank line
			est replay `eqn', noheader
		}
		di
	}
	// Display first-stage estimations.
	if "`e(firstlist)'" ~= "" {
		di as text "First-stage estimation(s):"
		foreach eqn in `e(firstlist)' {
			est replay `eqn', noheader
		}
		di
	}
	
	// Display weak-ID stats.
	if e(idstats) {
		di as text "Weak identification F stats (i.i.d.):"
		if e(weakid_chs_l)<. | e(weakid_chs_pl)<. {
			di as text "  Optimal Lasso IV(s):"        _col(28) as res %8.2f `e(weakid_chs_l)'
			di as text "  Optimal Post-Lasso IV(s):"   _col(27) as res %8.2f `e(weakid_chs_pl)'
		}		
		di as text "  Full IV set:"                _col(28) as res %8.2f `e(weakid)'
		if "`e(vce)'" == "robust" | "`e(vce)'" == "cluster" {
			if "`e(vce)'" == "robust" {
				di as text "Weak identification F stats (robust):"
			}
			else {
				di as text "Weak identification F stats (cluster-robust):"
			}
			if e(weakidr_chs_l)<. | e(weakidr_chs_pl)<. {
				di as text "  Optimal Lasso IV(s):"        _col(28) as res %8.2f `e(weakidr_chs_l)'
				di as text "  Optimal Post-Lasso IV(s):"   _col(27) as res %8.2f `e(weakidr_chs_pl)'
			}		
			di as text "  Full IV set:"                _col(28) as res %8.2f `e(weakidr)'
		}
		di
	}
	
	// Display sup-score stats.
	if "`e(ss_null)'`e(ss_null_l)'`e(ss_null_pl)'" ~= "" {
		di as text "Sup-score weak-identification-robust tests (method=" as res "`e(ss_method)'" as text ")"
		di as text "H0: b(`e(dendog)')=0   (" 100*e(ss_gamma) "% significance level)"
		if "`e(xhighdim)'"~="" {
			di as text "  Lasso-orthogonalized:      " as res "`e(ss_null_l)'"
			di as text "  Post-lasso-orthogonalized: " as res "`e(ss_null_pl)'"
			di as text "  No orthogonalization:      " _c
			if "`e(ss_null)'" ~= "" {
				di as res "`e(ss_null)'"
			}
			else {
				di as res "n.a."
			}
		}
		else {
			di as res "  `e(ss_null)'"
		}
		if "`e(ss_cset)'`e(ss_cset_l)'`e(ss_cset_pl)'" ~= "" {
			di as text "`e(ss_level)'% confidence set (grid min=" %-3.2f e(ss_gridmin) ", max=" %-3.2f e(ss_gridmax) ", points=" e(ss_gridpoints) "):"
			if "`e(xhighdim)'"~="" {
				di as text "  Lasso-orthogonalized:      " as res "`e(ss_cset_l)'"
				di as text "  Post-lasso-orthogonalized: " as res "`e(ss_cset_pl)'"
				di as text "  No orthogonalization:      " _c
				if "`e(ss_cset)'" ~= "" {
					di as res "`e(ss_cset)'"
				}
				else {
					di as res "n.a."
				}
			}
			else {
				di as res "  `e(ss_cset)'"
			}
		}
		di
	}
	if ~e(ss_omitgrid) {
		tempname sstable
		mat `sstable'=e(ss_gridmat)
		if `sstable'[1,1]~=. {			// table exists
			di as text "Sup-score test results from user-supplied grid of hypothesized betas:"
			_matrix_table `sstable'
			di
		}
	}

	// Display structural equation results.
	local method		`e(method)'
	local depname		`e(depvar)'
	local N				=e(N)
	local vce			`e(vce)'
	local vcetype		`e(vcetype)'
	local clustvar		`e(clustvar)'
	local N_clust		=e(N_clust)
	local N_g			=e(N_g)

	tempname beta_pds beta_lasso beta_plasso V_pds V_lasso V_plasso
	mat `beta_pds'		=e(beta_pds)
	mat `V_pds'			=e(V_pds)
	mat `beta_lasso'	=e(beta_lasso)
	mat `V_lasso'		=e(V_lasso)
	mat `beta_plasso'	=e(beta_plasso)
	mat `V_plasso'		=e(V_plasso)
	// Check here for id failure for L and PL, before posting.
	// Can arise if opt IVs are collinear.
	// ID check for PDS is separate - id failure means ALL coefs unidentified
	local d0ct_l		=diag0cnt(`V_lasso')
	local d0ct_pl		=diag0cnt(`V_plasso')
	local d0ct_pds		=diag0cnt(`V_pds')
	local dendog_ct		=e(dendog_ct)
	local zselected_ct	=e(zselected_ct)
	tempname esthold
	_estimates hold `esthold'				//  temporarily store existing results
	if `N_g' {
		di as text "Structural equation (fixed effects, #groups=" as res `N_g' as text "):"
	}
	else {
		di as text "Structural equation:"
	}
	if el(`beta_lasso',1,1)<. {
		di
		di as text "`estimator' using CHS `method'-orthogonalized vars"
		ereturn post `beta_lasso' `V_lasso', obs(`N') depname(`depname')
		ereturn local vce		`vce'
		ereturn local vcetype	`vcetype'
		ereturn local clustvar	`clustvar'
		ereturn scalar N_clust	=`N_clust'
		ereturn di
		if `d0ct_l' {
			di as text "Warning: some coefficients unidentified."
		}
	}
	if el(`beta_plasso',1,1)<. {
		di
		di as text "`estimator' using CHS post-`method'-orthogonalized vars"
		ereturn post `beta_plasso' `V_plasso', obs(`N') depname(`depname')
		ereturn local vce		`vce'
		ereturn local vcetype	`vcetype'
		ereturn local clustvar	`clustvar'
		ereturn scalar N_clust	=`N_clust'
		ereturn di
		if `d0ct_pl' {
			di as text "Warning: some coefficients unidentified."
		}
	}
	if el(`beta_pds',1,1)<. {
		di
		di as text "`estimator' with PDS-selected variables and full regressor set"
		ereturn post `beta_pds' `V_pds', obs(`N') depname(`depname')
		ereturn local vce		`vce'
		ereturn local vcetype	`vcetype'
		ereturn local clustvar	`clustvar'
		ereturn scalar N_clust	=`N_clust'
		ereturn di
		// #endog > #excluded => PDS equation unidentified (CHS may still be identified)
		if `dendog_ct' > `zselected_ct' {
			di as text "Warning: not enough instruments selected; model unidentified"
		}
	}
	_estimates unhold `esthold'

end


*************** Select controls/IVs subroutine *********************

// returns r(allselected0_o) - all selected INCLUDING PARTIAL, NOTPEN, CONS
//         r(allselected_o)  - all selected penalized only, EXCLUDING PARTIAL, NOTPEN, CONS
// also generates fitted values or residuals
program define SelectControls, rclass
	version 13
	syntax ,						///
		touse(string)				/// pass all varlists as strings to stop fvs processed etc.
		modelvars_o(string)			/// can be y, x or z
		modelvars_t(string)			///
		dict_o(string)				/// dictionary
		dict_t(string)				/// dictionary
		[							///
		genvars_l(string)			/// can be a varlist if multiple explanatory vars; returns either
		genvars_pl(string)			/// predicted value or residual; empty => no xb/resid needed
		highdim_o(string)			/// can be hi-dim Xs or hi-dim Zs
		highdim_t(string)			///
		resid						/// generate residuals
		xb							/// generate fitted values
		aset_o(string)				///
		aset_t(string)				///
		notpen_o(string)			///
		notpen_t(string)			///
		partial_o(string)			///
		partial_t(string)			///
		loptions(string)			///
		msg(string)					///
		startcol(int 15)			///
		qui							///
		debugflag(int 0)			///
		lscons(int 0)				///
		]
	
	// create flags etc., then syntax check
	local residflag	=("`resid'"~="" | "`resid'`xb'"=="")		//  default is resid
	local xbflag	="`xb'"~=""
	// convenience macro
	if `lscons' {
		local consname _cons
	}
	else if "`partial_o'"=="" {									//  nocons and partial not compatible
		local noconstant noconstant
	}
	if `residflag' & `xbflag' {									//  can't both be ==1
		di as err "internal ivlasso/SelectControls error"
		exit 198
	}
	*

	tempname beta betaOLS

	// loop through model variables, selecting controls and constructing gen vars.
	// index i used because several lists: modelvars_o, modelvars_t, gen_l, gen_pl
	local numvars	: word count `modelvars_o'
	forvalues i=1/`numvars' {

		local var_o		: word `i' of `modelvars_o'
		local var_t		: word `i' of `modelvars_t'
		local gen_l		: word `i' of `genvars_l'
		local gen_pl	: word `i' of `genvars_pl'
		// initialize count of selected
		local sL	=0
		local sPL	=0

		if "`highdim_o'`notpen_o'`partial_o'"~="" & "`msg'"~="" {		//  report rlasso selections if highdim list not empty
			di as text "`msg' `var_o'..."
		}
		if "`qui'"=="" {
				di as text "rlasso estimation:"							//  Note in noi output that results are lasso
		}
		
		if "`highdim_o'`notpen_o'`partial_o'"~="" {						//  Enter if any regressors provided
			// rlasso reports OLS if highdim empty and only notpen or partial provided
			`qui' rlasso `var_t' `highdim_t' `notpen_t' `partial_t'		///
									if `touse',							///
									pnotpen(`notpen_t')					/// may be empty
									partial(`partial_t')				/// may be empty
									displaynames_o(`dict_t')			/// original names for rlasso are those in varlist
									displaynames_d(`dict_o')			/// display names are ivlasso original names
									`noconstant'						///
									`loptions'
			mat `beta'				=e(beta)							//  rlasso saves row vectors
			mat `betaOLS'			=e(betaOLS)							//  rlasso saves row vectors
			local betavars_o		`e(selected0)'						//  e(selected0) macro INCLUDES PARTIAL, NOTPEN AND CONS
			local sL				=e(s0)
			local sPL				=`sL'								//  will augment post-lasso count below
			local allselected0_o	`allselected0_o' `betavars_o'		//  add to full set of selected variables INCL PARTIAL, NOTPEN, CONS
				
			if "`highdim_o'`notpen_o'`partial_o"~="" {					//  report rlasso selections if anything to report
				// Report selected (penalized)
				local selected_o	`e(selected)'						//  e(selected) macro EXCLUDES notpen, partial, cons
				local allselected_o	`allselected_o' `selected_o'		//  allselected_o EXCLUDES notpen, partial, cnos
				di as text "Selected: " _c
				DispVars `selected_o', _col(`startcol')
				// Report additional amelioration set vars.
				local also_aset		: list aset_o - selected_o			//  also_aset = additional aset vars not in selected
				local allselected0_o `allselected0_o' `also_aset'		//  and add to allselected0 list
				local alsocount		: word count `also_aset'
				if `alsocount' {
					local sPL	=`sPL'+`alsocount'						//  augment count of selected vars for post-lasso
					di as text "Also inc: " _c							//  don't add aset to selected_o since not in beta
					DispVars `also_aset', _col(`startcol')
				}
				// Report additional unpenalized vars
				local alsocount		: word count `notpen_o' `partial_o'
				if `alsocount' {
					di as text "Also inc: " _c							//  notpen is not included in e(selected)
					DispVars `notpen_o' `partial_o', _col(`startcol')
				}
			}
		}
		else {															//  no regressors provided
			local sL	=0
			local sPL	=0
		}

		// Summary: beta and betaOLS available unless no hi-dim selected and no notpen vars.
		//          sL  = count of hi-dim selected + notpen
		//          sPL = count of hi-dim selected + notpen + aset
		//          selected_o = lasso-selected vars only
		//          betavars_o = vars corresponding to vars in beta and betaOLS, including notpen but NOT aset
		
		// first generate lasso fitted value, then replace with resid if requested
		if "`gen_l'" ~= "" {												//  empty means nothing to generate
			if `sL'	{														//  lasso selected something; could be hi-dim/notpen
				matchnames "`betavars_o'" "`dict_o'" "`dict_t'"				//  set up lasso coeffs for mat score
				mat colnames `beta'		= `r(names)'						//  for mat score to work, need to have temp names
				// first generate fitted value, then replace with resid if requested
				qui matrix score `gen_l'	= `beta' if `touse', replace
				if `residflag' {
					qui replace `gen_l'		= `var_t' - `gen_l'
				}
			}
			else {															//  lasso selected nothing or no hi-dims/notpens
				if `residflag' {
					qui replace `gen_l'		= `var_t'						//  "residual" is orig variable
				}
				else {														//  "fitted value" is zero vector
					qui replace `gen_l'		= 0
				}
			}
		}
		// first generate post-lasso fitted value, then replace with resid if requested
		if "`gen_pl'" ~= "" {												//  empty means nothing to generate
			if `sPL'	{													//  post-lasso selected something; could be hi-dim/aset/notpen
				if `sPL'==`sL' {											//  no aset vars so no need to re-run post-lasso
					matchnames "`betavars_o'" "`dict_o'" "`dict_t'"			// Set up post-lasso coeffs for mat score
					mat colnames `betaOLS'		= `r(names)'				//  for mat score to work, need to have temp names
					qui matrix score `gen_pl'	= `betaOLS' if `touse', replace
				}
				else {														//  need to re-run post-lasso including aset vars
					matchnames "`betavars_o' `also_aset'" "`dict_o'" "`dict_t'"
					local regressors `r(names)'
					local regressors : list regressors - consname
					qui regress `var_t' `regressors' `aset_t' if `touse', `noconstant'
					tempvar xb
					qui predict double `xb' if `touse', xb
					qui replace `gen_pl' = `xb'
				}
				if `residflag' {
					qui replace `gen_pl'	= `var_t' - `gen_pl'
				}
			}
			else {															//  post-lasso selected nothing or no hi-dims
				if `residflag' {
					qui replace `gen_pl'	= `var_t'						//  "residual" is orig variable
				}
				else {														//  "fitted value" is zero vector
					qui replace `gen_pl'		= 0
				}
			}
		}

	}	//  end of loop through modelvars

	// remove duplicates
	local allselected_o				: list uniq allselected_o
	local allselected0_o			: list uniq allselected0_o
	// return selected lists and counts
	return local allselected_o		`allselected_o'							//  EXCLUDING NOTPEN, PARTIAL, CONS
	return local allselected0_o		`allselected0_o'						//  INCLUDING NOTPEN, PARTIAL, CONS
	return scalar s					=`: word count `allselected_o''
	return scalar s0				=`: word count `allselected0_o''

end		// End program SelectControls

// idstats - returns weak ID stats for 3 IV estimations (CHS L and PL, PDS)
// return classical and robust/cluster ("r") versions
// expects tempvars etc. - no FV or TS vars allowed (ranktest won't take them)
// inexog and exexog arrive after collinearities/omitted/etc. removed
prog define doIDstats, rclass
	syntax [if] [in],						///
		endog(string)						/// string (ensures no processing)
		[									///
		endog_l(string)						///
		endog_pl(string)					///
		exexog(string)						///
		exexog_l(string)					///
		exexog_pl(string)					///
		inexog(string)						///
		inexog_ct(integer 0)				///
		cons(int 0)							///
		robust								///
		cluster(string)						///
		]

	// ranktest is a necessary component
	capture ranktest, version
	if _rc != 0 {
		di as err "Error: must have ranktest installed for idstats option"
		di as err "To install, from within Stata type " _c
		di in smcl "{stata ssc install ranktest :ssc install ranktest}"
		exit 601
	}
		
	marksample touse
	
	if ~`cons' {
		local noconstant noconstant
	}
	
	// PDS is full IV set so constant is possible
	local exexog_ct	: word count `exexog'
	local endog_ct	: word count `endog'
	if `endog_ct' > `exexog_ct' {
		di as text "Warning - idstats error; model unidentified"
	}
	else {
		cap ranktest 	(`endog') (`exexog') if `touse'			///
						, partial(`inexog') 					///
						full									///
						wald									///
						`noconstant'
		if _rc {
			di as text "Warning - idstats/ranktest error, singular matrices possible"
		}
		else {
			local weakid		=r(chi2)/r(N)*(r(N)-`inexog_ct'-`exexog_ct')/`exexog_ct'
		}
		if "`robust'`cluster'"~="" {
			cap ranktest 	(`endog') (`exexog') if `touse'			///
							, partial(`inexog') 					///
							full									///
							wald									///
							`robust' cluster(`cluster')				///
							`noconstant'
			if _rc {
				di as text "Warning - idstats/ranktest error, singular matrices possible"
			}
			else if "`robust'"~="" {
				local weakidr	=r(chi2)/r(N)*(r(N)-`inexog_ct'-`exexog_ct')/`exexog_ct'
			}
			else if "`cluster'"~="" {
				local weakidr	=r(chi2)/r(N)*(r(N)-`inexog_ct'-`exexog_ct')/`exexog_ct' * (r(N_clust)-1)/r(N_clust)
			}
		}
	}
	// endog_l is optimal lasso IV so always nocons
	// optimal IVs always provided so order condition always satisfied
	local exexog_l_ct	: word count `exexog_l'
	cap ranktest 	(`endog_l') (`exexog_l') if `touse'		///
					, partial(`inexog') 					///
					full									///
					wald									///
					nocons
	if _rc {
		di as text "Warning - idstats/ranktest error, singular matrices possible"
	}
	else {
		local weakid_l		=r(chi2)/r(N)*(r(N)-`inexog_ct'-`exexog_l_ct')/`exexog_l_ct'
	}
	if "`robust'`cluster'"~="" {
		cap ranktest 	(`endog_l') (`exexog_l') if `touse'		///
						, partial(`inexog') 					///
						full									///
						wald									///
						`robust' cluster(`cluster')				///
						nocons
		if _rc {
			di as text "Warning - idstats/ranktest error, singular matrices possible"
		}
		else if "`robust'"~="" {
			local weakidr_l	=r(chi2)/r(N)*(r(N)-`inexog_ct'-`exexog_l_ct')/`exexog_l_ct'
		}
		else if "`cluster'"~="" {
			local weakidr_l	=r(chi2)/r(N)*(r(N)-`inexog_ct'-`exexog_l_ct')/`exexog_l_ct' * (r(N_clust)-1)/r(N_clust)
		}
	}
	
	// endog_pl is optimal post-lasso IV so always nocons
	// optimal IVs always provided so order condition always satisfied
	local exexog_pl_ct	: word count `exexog_pl'
	cap ranktest 	(`endog_pl') (`exexog_pl') if `touse'	///
					, partial(`inexog') 					///
					full									///
					wald									///
					nocons
	if _rc {
		di as text "Warning - idstats/ranktest error, singular matrices possible"
	}
	else {
		local weakid_pl		=r(chi2)/r(N)*(r(N)-`inexog_ct'-`exexog_pl_ct')/`exexog_pl_ct'
	}
	if "`robust'`cluster'"~="" {
		cap ranktest 	(`endog_pl') (`exexog_pl') if `touse'	///
						, partial(`inexog') 					///
						full									///
						wald									///
						`robust' cluster(`cluster')				///
						nocons
		if _rc {
			di as text "Warning - idstats/ranktest error, singular matrices possible"
		}
		else if "`robust'"~="" {
			local weakidr_pl		=r(chi2)/r(N)*(r(N)-`inexog_ct'-`exexog_pl_ct')/`exexog_pl_ct'
		}
		else else if "`cluster'"~="" {
			local weakidr_pl		=r(chi2)/r(N)*(r(N)-`inexog_ct'-`exexog_pl_ct')/`exexog_pl_ct' * (r(N_clust)-1)/r(N_clust)
		}
	}

	return local weakid		`weakid'	
	return local weakid_l	`weakid_l'
	return local weakid_pl	`weakid_pl'
	return local weakidr	`weakidr'	
	return local weakidr_l	`weakidr_l'
	return local weakidr_pl	`weakidr_pl'
	
end

// regress utility
// uses old-style reg y x1 x2 (x1 z1) to do IV estimation
program define doregress, eclass
	version 13
	syntax varlist(numeric fv ts min=1 max=1) [if] [in],	/// dep var
		regressors(varlist numeric fv ts min=1)				/// all regressors
		instruments(varlist numeric fv ts min=1)			/// all exogenous
		[													///
		dict_o(string)										/// original varnames
		dict_t(string)										/// temporary varnames
		cons(int 1)											/// flag for constant
		dofminus(integer 0)									/// large-sample dof adjustment (e.g. FEs)
		debugflag(int 0)									///
		*													///  all options
		]

	marksample touse
	local varY_t	`varlist'
	
	if ~`debugflag' {
		local qui qui
	}

	if ~`cons' {
		local noconstant noconstant
	}
	
	tempname b V
	// capture in case not identified
	cap regress `varY_t'								///
		`regressors' (`instruments')					///
	 	if `touse',										///
	 	`noconstant'									///
		`options'
	
	// eqn identified
	if _rc==0 {
		if `debugflag' {
			ereturn list
		}
	
		// Extract options to be passed back or used below
		local N			=`e(N)'
		local df_r		=`e(df_r)'
		local df_m		=`e(df_m)'
		local rank		=`e(rank)'
		local vcetype	`e(vcetype)'
		local vce		`e(vce)'
		local clustvar	`e(clustvar)'
		if "`clustvar'"=="" {
			local N_clust	=0
		}
		else {
			local N_clust	=`e(N_clust)'
		}
	
		mat `b'	=e(b)
		mat `V'	=e(V)
	
		// Fix (remove) regress dof
		if "`clustvar'"=="" {							//  remove basic N-K dof adjustment
			mat `V' = `V' * (`N'-`rank') / (`N'-`dofminus')
		}
		else {											//  remove cluster dof adjustment
			mat `V' = `V' * (`N'-`rank') / (`N'-1)		///
						* (`N_clust'-1) / `N_clust'
		}
		
		// Fix tempnames if dictionary supplied.
		if "`dict_o'" ~= "" {
			matchnames "`varY_t'" "`dict_t'" "`dict_o'"
			local varY_o	`r(names)'
			local cnames_t	: colnames `b'					//  tempnames so no FV notation (otherwise b/n/o would get inserted!)
			matchnames "`cnames_t'" "`dict_t'" "`dict_o'"
			local cnames_o	`r(names)'
			mat colnames `b' = `cnames_o'
			mat colnames `V' = `cnames_o'
			mat rownames `V' = `cnames_o'
		}
		else {
			local varY_o	`e(depvar)'
		}
		// build FV info
		_ms_build_info	`b' if `touse'
	
		// Post b and V
		ereturn post `b' `V', esample(`touse') depname(`varY_o') obs(`N')
	
		// Needed to get est store to work
		ereturn local cmd 		ivlasso
		// Other results to return
		ereturn local vcetype	`vcetype'
		ereturn local vce		`vce'
		ereturn local clustvar	`clustvar'
		ereturn scalar N		=`N'
		ereturn scalar N_clust	=`N_clust'
	}
	else {
		// eqn unidentified

		qui count if `touse'
		local N		=r(N)
		
		matchnames "`varY_t'" "`dict_t'" "`dict_o'"
		local varY_o	`r(names)'
		matchnames "`regressors'" "`dict_t'" "`dict_o'"
		local colnames	`r(names)'
		if `cons' {
			local colnames	`colnames' _cons
		}
		local k	: word count `colnames'
				
		mat `b'	=J(1,`k',0)
		mat `V'	=J(`k',`k',0)
		local cnames_o	`r(names)'
		mat colnames `b' = `cnames_o'
		mat colnames `V' = `cnames_o'
		mat rownames `V' = `cnames_o'
		
		// build FV info
		_ms_build_info	`b' if `touse'
	
		// Post b and V
		ereturn post `b' `V', esample(`touse') depname(`varY_o') obs(`N')
	
		// Needed to get est store to work
		ereturn local cmd 		ivlasso
	
	}

end

// sup-score tests
program define doSupScore, rclass
	version 13
	syntax [,								///
			touse(string)					///
			varY_o(string)					///
			varY_t(string)					///
			rho_y_o(string)					///
			rho_y_l(string)					///
			rho_y_pl(string)				///
			dendog_o(string)				///
			dendog_t(string)				///
			xhighdim_o(string)				///
			xhighdim_t(string)				///
			zhighdim_o(string)				///
			zhighdim_t(string)				///
			xnotpen_o(string)				///
			xnotpen_t(string)				///
			xaset_o(string)					///
			xaset_t(string)					///
			dict_o(string)					///
			dict_t(string)					///
			loptions(string)				///
			sscsetflag(integer 0)			///
			ssmethod(name)					///
			ssgamma(real 0.05)				///
			gridmin(numlist missingok)		///
			gridmax(numlist missingok)		///
			gridpoints(integer 100)			///
			gridmat(name)					///
			startcol(integer 1)				///
			debugflag(integer 0)			///
			lscons(integer 0)				///
			dofopt(string)					///
			ivoptions(string)				///
			cons(integer 1)					///
			qui								///
			]

		if "`ssmethod'"=="" {
			local ssmethod			abound						//  default
		}
		if "`ssmethod'"=="abound" {								//  use conservative asymptotic bound, collect rejections
			local sssimulateflag	=0
			local ssaboundflag		=1
			local ssselectflag		=0
		}
		else if "`ssmethod'"=="simulate" {						//  simulate distrib of sup-score stat
			local sssimulateflag	=1
			local ssaboundflag		=0
			local ssselectflag		=0
			local sslevel			level(`ssgamma')			//  collect pvalues
		}
		else if "`ssmethod'"=="select" {						//  reject if any vars selected, collect rejections
			local sssimulateflag	=0
			local ssaboundflag		=0
			local ssselectflag		=1
		}
		else {
			di as err "error - invalid ssmethod(.) option"
			exit 198
		}
		
		// trick to parse the lasso options
		local 0 ", `loptions'"
		// this strips out gamma and gammad options in loptions macro
		// and leaves remainder in macro `options'
		syntax [, gamma(real 0.05) gammad(real 1) * ]
		// set options for rlasso
		if `sssimulateflag' {
			local ssoptions		testonly gamma(`ssgamma') `options'
		}
		else if `ssaboundflag' {
			local ssoptions		testonly ssnumsim(0) gamma(`ssgamma') `options'
		}
		else if `ssselectflag' {
			local ssoptions		gamma(`ssgamma') gammad(1) `options'
		}
		else {
			di as err "error - internal ivlasso sup-score error"
			exit 198
		}

		tempname supscoremat
		// initialize
		local dendog_ct		: word count `dendog_o'
		tempvar ytilde
		// ytilde is a single variable
		qui gen double `ytilde'=.
		// rho_e_l and rho_e_pl can be varlists
		forvalues i=1/`dendog_ct' {
			tempvar v_l v_pl
			qui gen double `v_l'=.
			qui gen double `v_pl'=.
			local rho_e_l	`rho_e_l'  `v_l'
			local rho_e_pl	`rho_e_pl' `v_pl'
		}
		
		************** user-supplied grid for search ******************
		if "`gridmat'"~="" {
			local gridrows		=rowsof(`gridmat')
			local gridcols		=colsof(`gridmat')
			// check dimensions
			if `gridcols'~=`dendog_ct' {
				di as err "Error - number of gridmat columns does not match number of endogenous regressors"
				exit 503
			}
			// vector for rejections (1/0)
			tempname gridreject gridreject_l gridreject_pl

			// grid search
			_dots 0 0, title(Estimating sup-score confidence set over `gridrows' points in grid `gridmat')
			local counter 1
			forvalues rownum=1/`gridrows' {
				// standard, using unorthogonalized y and d
				// high-dim controls may be present
				qui replace `ytilde'	= `varY_t' if `touse'
				tokenize `dendog_t'
				forvalues colnum=1/`gridcols' {
					qui replace `ytilde' = `ytilde'-el(`gridmat',`rownum',`colnum')*``colnum'' if `touse'
				}
				qui SelectControls,										/// variable selection only, no vars generated, xaset unnecessary
								touse(`touse')							///
								modelvars_o(`varY_o')					/// just to give it a name
								modelvars_t(`ytilde')					///
								highdim_o(`xhighdim_o' `zhighdim_o')	///
								highdim_t(`xhighdim_t' `zhighdim_t')	///
								notpen_o(`xnotpen_o')					///
								notpen_t(`xnotpen_t')					///
								dict_o(`dict_o')						/// dictionary
								dict_t(`dict_t')						/// dictionary
								loptions(`ssoptions')					/// use sup-score options
								startcol(`startcol')					///
								debugflag(`debugflag')					///
								lscons(`lscons')						///
								`qui'
				if `sssimulateflag' {
					mat `gridreject'	= nullmat(`gridreject') \ e(supscore_p)
				}
				else if `ssaboundflag' {
					local r				= e(supscore)>e(supscore_cv)
					mat `gridreject'	= nullmat(`gridreject') \ `r'
				}
				else {
					local ss_selected	`r(allselected_o)'					//  use only penalized IVs
					local ss_selected	: list ss_selected - xhighdim_o
					local ss_selected	: list ss_selected - xnotpen_o
					local ss_s			: word count `ss_selected'
					local r				= `ss_s'>0
					mat `gridreject'	= nullmat(`gridreject') \ `r'
				}
				
				if "`xhighdim_o'" ~= "" {
					// using lasso- and post-lasso-orthogonalized y and d
					// first orthogonalize d (using same settings as for y)
					qui SelectControls,										///
									touse(`touse')							///
									modelvars_o(`dendog_o')					///
									modelvars_t(`dendog_t')					///
									genvars_l(`rho_e_l')					/// updates to resid or to y if nothing selected
									genvars_pl(`rho_e_pl')					/// updates to resid or to y if nothing selected
									highdim_o(`xhighdim_o')					///
									highdim_t(`xhighdim_t')					///
									notpen_o(`xnotpen_o')					/// xnotpen needed for equivalent with partialling-out
									notpen_t(`xnotpen_t')					///
									dict_o(`dict_o')						/// dictionary
									dict_t(`dict_t')						/// dictionary
									aset_o(`xaset_o')						///
									aset_t(`xaset_t')						///
									loptions(`loptions')					/// use same settings as for main estimation
									startcol(`startcol')					///
									debugflag(`debugflag')					///
									lscons(`lscons')

					// sup-score using lasso-orthogonalized variables
					qui replace `ytilde'	= `rho_y_l' if `touse'
					tokenize `rho_e_l'
					forvalues colnum=1/`gridcols' {
						qui replace `ytilde' = `ytilde'-el(`gridmat',`rownum',`colnum')*``colnum'' if `touse'
					}
					qui SelectControls,										/// variable selection only, no vars generated, xaset unnecessary
									touse(`touse')							///
									modelvars_o(`varY_o')					///
									modelvars_t(`ytilde')					///
									highdim_o(`zhighdim_o')					///
									highdim_t(`zhighdim_t')					///
									notpen_o(`xnotpen_o')					/// xnotpen needed for equivalent with partialling-out
									notpen_t(`xnotpen_t')					///
									dict_o(`dict_o')						/// dictionary
									dict_t(`dict_t')						/// dictionary
									loptions(`ssoptions')					/// use sup-score options
									startcol(`startcol')					///
									debugflag(`debugflag')					///
									lscons(`lscons')						///
									`qui'
					if `sssimulateflag' {
						mat `gridreject_l'	= nullmat(`gridreject_l') \ e(supscore_p)
					}
					else if `ssaboundflag' {
						local r				= e(supscore)>e(supscore_cv)
						mat `gridreject_l'	= nullmat(`gridreject_l') \ `r'
					}
					else {
						local r				= e(s)>0
						mat `gridreject_l'	= nullmat(`gridreject_l') \ `r'
					}

					// sup-score using post-lasso-orthogonalized variables
					qui replace `ytilde'	= `rho_y_pl' if `touse'
					tokenize `rho_e_pl'
					forvalues colnum=1/`gridcols' {
						qui replace `ytilde' = `ytilde'-el(`gridmat',`rownum',`colnum')*``colnum'' if `touse'
					}
					qui SelectControls,										/// variable selection only, no vars generated, xaset unnecessary
									touse(`touse')							///
									modelvars_o(`rho_y_o')					///
									modelvars_t(`ytilde')					///
									highdim_o(`zhighdim_o')					///
									highdim_t(`zhighdim_t')					///
									notpen_o(`xnotpen_o')					/// xnotpen needed for equivalent with partialling-out
									notpen_t(`xnotpen_t')					///
									dict_o(`dict_o')						/// dictionary
									dict_t(`dict_t')						/// dictionary
									loptions(`ssoptions')					/// use sup-score options
									startcol(`startcol')					///
									debugflag(`debugflag')					///
									lscons(`lscons')						///
									`qui'
					if `sssimulateflag' {
						mat `gridreject_pl'	= nullmat(`gridreject_pl') \ e(supscore_p)
					}
					else if `ssaboundflag' {
						local r				= e(supscore)>e(supscore_cv)
						mat `gridreject_pl'	= nullmat(`gridreject_pl') \ `r'
					}
					else {
						local r				= e(s)>0
						mat `gridreject_pl'	= nullmat(`gridreject_pl') \ `r'
					}
				}
				// display dot upon completing loop; increment counter
				_dots `counter' 0
				local ++counter
			}		// end loop over user grid
		
			// gridmat augmented by rejection column(s)
			tempname gridresults
			if "`xhighdim_o'" == "" {
				mat `gridresults'			= `gridmat', `gridreject'
				if `sssimulateflag' {
					mat colnames `gridresults'	= `dendog_o' pvalue_L
				}
				else {
					mat colnames `gridresults'	= `dendog_o' reject_L
				}
			}
			else {
				mat `gridresults'			= `gridmat', `gridreject_l', `gridreject_pl', `gridreject'
				if `sssimulateflag' {
					mat colnames `gridresults'	= `dendog_o' pvalue_L pvalue_PL pvalue
				}
				else {
					mat colnames `gridresults'	= `dendog_o' reject_L reject_PL reject
				}
			}
			forvalues i=1/`gridrows' {
				local rnames `rnames' Null_`i'
			}
			mat rownames `gridresults'		= `rnames'

		}	// end user-supplied grid code
		*
		
		****************** test default null that beta=0 **************************
		// Use un-orthognalized vars. Sup-score statistic unavailable if any
		// high-dim Xs present, unless select method is used.
		if `ssselectflag' | ("`xhighdim_o'"=="") {
				qui SelectControls,										/// no vars generated, xaset unnecessary
								touse(`touse')							///
								modelvars_o(`varY_o')					///
								modelvars_t(`varY_t')					///
								highdim_o(`xhighdim_o' `zhighdim_o')	///
								highdim_t(`xhighdim_t' `zhighdim_t')	///
								notpen_o(`xnotpen_o')					///
								notpen_t(`xnotpen_t')					///
								dict_o(`dict_o')						/// dictionary
								dict_t(`dict_t')						/// dictionary
								loptions(`ssoptions')					/// include gammad=1 option
								startcol(`startcol')					///
								debugflag(`debugflag')					///
								lscons(`lscons')						///
								`qui'
			if `sssimulateflag' {
				if e(supscore_p)<`ssgamma' {
					local ss_null reject
				}
				else {
					local ss_null fail to reject
				}
			}
			else if `ssaboundflag' {
				if e(supscore)>e(supscore_cv) {
					local ss_null reject
				}
				else {
					local ss_null fail to reject
				}
			}
			else {
				local ss_selected	`r(allselected_o)'					//  use only penalized
				local ss_selected	: list ss_selected - xhighdim_o		//  remove high-dim Xs from selected list
				local ss_selected	: list ss_selected - xnotpen_o		//  unnecessary?
				local ss_s			: word count `ss_selected'
				local r				= `ss_s'>0
				if `r' {
					local ss_null reject
				}
				else {
					local ss_null fail to reject
				}
			}
		}	// end un-orthognalized vars block. locals will be empty if sup-score + xhighdim present.

		// Orthogonalization approach if X highdim present;
		// test is a0=0 so no orthgonalization of d_endog is needed.
		if "`xhighdim_o'" ~= "" {
			// lasso-orthogonalized
			qui SelectControls,										/// variable selection only, no vars generated
							touse(`touse')							///
							modelvars_o(`rho_y_o')					///
							modelvars_t(`rho_y_l')					///
							highdim_o(`zhighdim_o')					///
							highdim_t(`zhighdim_t')					///
							notpen_o(`xnotpen_o')					/// xnotpen needed for equivalent with partialling-out
							notpen_t(`xnotpen_t')					///
							dict_o(`dict_o')						/// dictionary
							dict_t(`dict_t')						/// dictionary
							loptions(`ssoptions')					/// include gammad=1 option
							startcol(`startcol')					///
							debugflag(`debugflag')					///
							lscons(`lscons')						///
							`qui'
			if `sssimulateflag' {
				if e(supscore_p)<`ssgamma' {
					local ss_null_l reject
				}
				else {
					local ss_null_l fail to reject
				}
			}
			else if `ssaboundflag' {
				if e(supscore)>e(supscore_cv) {
					local ss_null_l reject
				}
				else {
					local ss_null_l fail to reject
				}
			}
			else {
				if e(s) {											//  count only penalized IVs
					local ss_null_l reject
				}
				else {
					local ss_null_l fail to reject
				}
			}

			// post-lasso-orthogonalized
			qui SelectControls,											/// variable selection only, no vars generated
							touse(`touse')							///
							modelvars_o(`rho_y_o')					///
							modelvars_t(`rho_y_pl')					///
							highdim_o(`zhighdim_o')					///
							highdim_t(`zhighdim_t')					///
							notpen_o(`xnotpen_o')					/// xnotpen needed for equivalent with partialling-out
							notpen_t(`xnotpen_t')					///
							dict_o(`dict_o')						/// dictionary
							dict_t(`dict_t')						/// dictionary
							loptions(`ssoptions')					/// include gammad=1 option
							startcol(`startcol')					///
							debugflag(`debugflag')					///
							lscons(`lscons')						///
							`qui'
			if `sssimulateflag' {
				if e(supscore_p)<`ssgamma' {
					local ss_null_pl reject
				}
				else {
					local ss_null_pl fail to reject
				}
			}
			else if `ssaboundflag' {
				if e(supscore)>e(supscore_cv) {
					local ss_null_pl reject
				}
				else {
					local ss_null_pl fail to reject
				}
			}
			else {
				if e(s) {											//  count only penalized IVs
					local ss_null_pl reject
				}
				else {
					local ss_null_pl fail to reject
				}
			}

		}  // end test of default null

		************** construct CI if #endog=1 and requested ***********************
		if `sscsetflag' {
			// gridlimits are numlists - use first in list if it exists
			local gridmin	: word 1 of `gridmin'
			local gridmax	: word 1 of `gridmax'
			// ols coef to get default grid limits
			if "`gridmin'"=="" | "`gridmax'"=="" {
				doregress `varY_t' if `touse',							///
					regressors(`dendog_t' `xnotpen_t')					///
					instruments(`dendog_t' `xnotpen_t')					///
					dict_o(`dict_o')									///
					dict_t(`dict_t')									///
					cons(`cons')										///
					`dofopt'											///
					`ivoptions'											///
					debugflag(`debugflag')
				if "`gridmin'"=="" {
					local gridmin		= _b[`dendog_o']-20*_se[`dendog_o']
				}
				if "`gridmax'"=="" {
					local gridmax		= _b[`dendog_o']+20*_se[`dendog_o']
				}
			}
			local gridinterval = .999999999*(`gridmax'-`gridmin')/(`gridpoints'-1)
			local grid "`gridmin'(`gridinterval')`gridmax'"				//  grid is in numlist form
			tempname citable citable_l citable_pl
			_dots 0 0, title(Estimating sup-score confidence set over `gridpoints' grid points)
			local counter 1
			forvalues g=`grid' {
				// Use un-orthognalized vars. Sup-score statistic unavailable if any
				// high-dim Xs present, unless select method is used.
				if `ssselectflag' | ("`xhighdim_o'"=="") {
					qui replace `ytilde' = `varY_t'-`g'*`dendog_t' if `touse'
					qui SelectControls,										/// variable selection only, no vars generated, xaset unnecessary
									touse(`touse')							///
									modelvars_o(`varY_o')					/// just to give it a name
									modelvars_t(`ytilde')					///
									highdim_o(`xhighdim_o' `zhighdim_o')	///
									highdim_t(`xhighdim_t' `zhighdim_t')	///
									notpen_o(`xnotpen_o')					///
									notpen_t(`xnotpen_t')					///
									dict_o(`dict_o')						/// dictionary
									dict_t(`dict_t')						/// dictionary
									loptions(`ssoptions')					/// include gammad=1 option
									startcol(`startcol')					///
									debugflag(`debugflag')					///
									lscons(`lscons')						///
									`qui'
					if `sssimulateflag' {
						mat `citable'		= nullmat(`citable') \ (`g', e(supscore_p))
					}
					else if `ssaboundflag' {
						local r				= e(supscore)>e(supscore_cv)
						mat `citable'		= nullmat(`citable') \ (`g', `r')
					}
					else {
						local ss_selected	`r(allselected_o)'					//  use only penalized
						local ss_selected	: list ss_selected - xhighdim_o		//  remove high-dim Xs from selected list
						local ss_selected	: list ss_selected - xnotpen_o		//  unnecessary?
						local ss_s			: word count `ss_selected'
						local r				= `ss_s'>0
						mat `citable'	= nullmat(`citable') \ (`g', `r')
					}
				}	// end un-orthognalized vars block. locals will be empty if sup-score + xhighdim present.
				
				if "`xhighdim_o'" ~= "" {
					// using lasso- and post-lasso-orthogonalized y and d
					// first orthogonalize d (using same settings as for y)
					qui SelectControls,										///
									touse(`touse')							///
									modelvars_o(`dendog_o')					///
									modelvars_t(`dendog_t')					///
									genvars_l(`rho_e_l')					/// updates to resid or to y if nothing selected
									genvars_pl(`rho_e_pl')					/// updates to resid or to y if nothing selected
									highdim_o(`xhighdim_o')					///
									highdim_t(`xhighdim_t')					///
									notpen_o(`xnotpen_o')					/// xnotpen needed for equivalent with partialling-out
									notpen_t(`xnotpen_t')					///
									dict_o(`dict_o')						/// dictionary
									dict_t(`dict_t')						/// dictionary
									aset_o(`xaset_o')						///
									aset_t(`xaset_t')						///
									loptions(`loptions')					/// use same settings as for main estimation
									startcol(`startcol')					///
									debugflag(`debugflag')					///
									lscons(`lscons')
					// using lasso-orthogonalized variables
					qui replace `ytilde' = `rho_y_l'-`g'*`rho_e_l' if `touse'
					qui SelectControls,										/// variable selection only, no vars generated, xaset unnecessary
									touse(`touse')							///
									modelvars_o(`rho_y_o')					///
									modelvars_t(`ytilde')					///
									highdim_o(`zhighdim_o')					///
									highdim_t(`zhighdim_t')					///
									notpen_o(`xnotpen_o')					/// xnotpen needed for equivalent with partialling-out
									notpen_t(`xnotpen_t')					///
									dict_o(`dict_o')						/// dictionary
									dict_t(`dict_t')						/// dictionary
									loptions(`ssoptions')					/// include gammad=1 option
									startcol(`startcol')					///
									debugflag(`debugflag')					///
									lscons(`lscons')						///
									`qui'
						if `sssimulateflag' {
							mat `citable_l'		= nullmat(`citable_l') \ (`g', e(supscore_p))
						}
						else if `ssaboundflag' {
							local r				= e(supscore)>e(supscore_cv)
							mat `citable_l'		= nullmat(`citable_l') \ (`g', `r')
						}
						else {
							local r				= e(s)>0
							mat `citable_l'		= nullmat(`citable_l') \ (`g', `r')
						}

					// using post-lasso-orthogonalized variables
					qui replace `ytilde' = `rho_y_pl'-`g'*`rho_e_pl' if `touse'
					qui SelectControls,										/// variable selection only, no vars generated, xaset unnecessary
									touse(`touse')							///
									modelvars_o(`rho_y_o')					///
									modelvars_t(`ytilde')					///
									highdim_o(`zhighdim_o')					///
									highdim_t(`zhighdim_t')					///
									notpen_o(`xnotpen_o')					/// xnotpen needed for equivalent with partialling-out
									notpen_t(`xnotpen_t')					///
									dict_o(`dict_o')						/// dictionary
									dict_t(`dict_t')						/// dictionary
									loptions(`ssoptions')					/// include gammad=1 option
									startcol(`startcol')					///
									debugflag(`debugflag')					///
									lscons(`lscons')						///
									`qui'
						if `sssimulateflag' {
							mat `citable_pl'	= nullmat(`citable_pl') \ (`g', e(supscore_p))
						}
						else if `ssaboundflag' {
							local r				= e(supscore)>e(supscore_cv)
							mat `citable_pl'	= nullmat(`citable_pl') \ (`g', `r')
						}
						else {
							local r				= e(s)>0
							mat `citable_pl'	= nullmat(`citable_pl') \ (`g', `r')
						}
				}
				// display dot upon completing loop; increment counter
				_dots `counter' 0
				local ++counter
			}

			// get confidence set from table
			// if `sslevel' is empty, it's a table of rejections,
			// otherwise it's a table of p-values and rejection depends on test level
			// un-orthogonalized version unavailable unless select method used or no high-dim Xs
			if `ssselectflag' | ("`xhighdim_o'"=="") {
				get_ci_from_table, citable(`citable') `sslevel'
				local ss_cset		`r(ss_cset)'
			}
			// orthogonalized version available if there are high-dim Xs
			if "`xhighdim_o'" ~= "" {
				get_ci_from_table, citable(`citable_l') `sslevel'
				local ss_cset_l	`r(ss_cset)'
				get_ci_from_table, citable(`citable_pl') `sslevel'
				local ss_cset_pl	`r(ss_cset)'
			}
		}  // end construction of CI

		// finish up
		return local ss_null	`ss_null'
		return local ss_null_l	`ss_null_l'
		return local ss_null_pl	`ss_null_pl'
		return local ss_cset	`ss_cset'
		return local ss_cset_l	`ss_cset_l'
		return local ss_cset_pl	`ss_cset_pl'
		return local gridmin	`gridmin'		// can be missing/empty
		return local gridmax	`gridmax'		// can be missing/empty
		if `sscsetflag' {
			if "`ss_cset'"~="" & "`xhighdim_o'" ~= "" {
				// all 3 available
				mat `citable'				= `citable_l', `citable_pl'[1...,2], `citable'[1...,2]
				if "`ssmethod'"=="simulate" {
					mat colnames `citable'		= `dendog_o' pvalue_L pvalue_PL pvalue
				}
				else {
					mat colnames `citable'		= `dendog_o' reject_L reject_PL reject
				}
			}
			else if "`ss_cset'"=="" {
				// 2 orthogonalized versions available
				mat `citable'				= `citable_l', `citable_pl'[1...,2]
				if "`ssmethod'"=="simulate" {
					mat colnames `citable'		= `dendog_o' pvalue_L pvalue_PL
				}
				else {
					mat colnames `citable'		= `dendog_o' reject_L reject_PL
				}
			}
			else {
				// just the un-orthognalized version available, already in `citable'
				if "`ssmethod'"=="simulate" {
					mat colnames `citable'		= `dendog_o' pvalue
				}
				else {
					mat colnames `citable'		= `dendog_o' reject
				}
			}
			
			local rnames
			local rownum				= rowsof(`citable')
			forvalues i=1/`rownum' {
				local rnames `rnames' Null_`i'
			}
			mat rownames `citable'		= `rnames'
			return matrix citable		`citable'
		}
		if "`gridmat'"~="" {
			return matrix gridresults	`gridresults'
		}
		return local ss_method	`ssmethod'		// will be non-empty

end


// Adapted from weakiv
program define get_ci_from_table, rclass
	version 13
	syntax [,								///
				citable(name local)			///
				level(real 0)				///
			]

* Macro `citable' is name of Mata matrix with table of either p-values or rejections.
* Rejections table is then collapsed in Mata.
* "collapsed" = rows which are not on either side of CI border for any test are dropped.
* CIs are then constructed by looping through Stata rejections table `rtable'.

	tempname rtable												//  table of collapsed rejections
																//  name use for both Mata and Stata objects
	
	mata: collapse_citable("`citable'", `level')				//  create table of rejections and collapse (delete unneeded rows)
	mat `rtable' = r(rtable)

* create macros for storing confidence sets
	local ss_cset ""
	local ss_rbegin=0
	local ss_rend=0
	local ss_rbegin_null=0
	local ss_rend_null=0
	local ss_flag=0

	local rows		=rowsof(`rtable')

	forvalue i=1/`rows' {
		local gridnull	= el(`rtable',`i',1)
* write out confidence set from rejection indicators
		local ss_r = el(`rtable',`i',2)					//  rejection indicator for test				
		
		if `ss_r'< . {
			local ss_flag=1								//  at least one value not missing
		}
		if `ss_r'==0 {
			if `ss_rbegin'==0 {
				local ss_rbegin=`i'
				if `i'==1 {
					local ss_rbegin_null "   ...  "
				}
				else {
					local ss_rbegin_null : di %8.0g `gridnull'
				}
			}
			local ss_rend=`i'
			if `i'==`rows' {
				local ss_rend_null "   ...  "
				}
				else {
					local ss_rend_null : di %8.0g `gridnull'
				}
		}
		if `ss_r'==1 | (`ss_r'==0 & `i'==`rows') {
			if `ss_rbegin'>0 & `ss_rend'>0 & (`ss_rbegin'==`ss_rend' & `i'<`rows') {
				local rnull : di %8.0g "`ss_rbegin_null'"
				if length("`ss_cset'")==0	local ss_cset "`rnull'"
				else							local ss_cset "`ss_cset' U `rnull'"
				local ss_rbegin=0
				local ss_rend=0
			}
			else if `ss_rbegin'>0 & `ss_rend'>0 & (`ss_rbegin'<`ss_rend' | `i'==`rows') {
				local rnull1 "`ss_rbegin_null'"
				local rnull2 "`ss_rend_null'"
				if length("`ss_cset'")==0	local ss_cset "[`rnull1',`rnull2']"
				else							local ss_cset "`ss_cset' U [`rnull1',`rnull2']"
				local ss_rbegin=0
				local ss_rend=0
			}
		}

	}	// end loop over grid points
	if `ss_flag'==0 {
		local ss_cset "."					//  all grid values missing, cset==.
	}
	else if length("`ss_cset'")==0 {
		local ss_cset "null set"			//  never failed to reject
	}
	tokenize "`ss_cset'", parse(",[] ")
	local wcount : word count `*'
* If cset is "[   ...  ,   ...  ]" then it has 5 tokenized elements and #2=#4="..."
	if `wcount'==5 & "`2'"=="..." & "`4'"=="..." {
		local ss_cset "entire grid"			//  never rejected
	}
	
	return local ss_cset "`ss_cset'"


end		//  end get_ci_from_table


*************************** Stata utilities ******************************

// copied from lassoutils 1.0.03 17/12/2017
// subroutine for partialling out
program define _partial, rclass sortpreserve
	version 13
	syntax anything ,							/// anything is actually a varlist but this is more flexible - see below
		touse(varlist numeric max=1)			/// required `touse' variable (so marksample unneeded)
		[										/// 
		PARtial(string)							/// string so that fv operators aren't inserted
		tvarlist(string)						/// optional list of temp vars - may be "unfilled" (NaNs)
		cons(int 1)								///
		solver(string)							/// svd, qr or empty
		]

	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename.
	local varlist `anything'
	*

	*** Error check - if tvarlist provided, should have same number of elements as varlist
	local v_ct	: word count `varlist'
	local tv_ct	: word count `tvarlist'
	if `tv_ct' & (`v_ct'~=`tv_ct') {
		di as err "internal partialling error - mismatched lists"
		exit 198
	}
	*

	*** recode solver macro into numeric 1/2/3
	if "`solver'"=="" {
		local solver 0			//  default method (SVD, QR if collinearities - see below)
	}
	else if "`solver'"=="svd" {
		local solver 1			//  force use of SVD solver
	}
	else if "`solver'"=="qr" {
		local solver 2			//  force use of QR solver
	}
	else {
		di as err "Syntax error: solver `solver' not allowed"
		exit 198
	}
	*
	
	*** partial out
	mata: s_partial("`varlist'", 			/// y and X
					"`partial'",          	/// to be partialled out
					"`tvarlist'",			///
					"`touse'",				/// touse
					`cons', 				/// treatment of constant
					`solver')				//  choice of solver (optional)
	return scalar rank	=`r(rank)'			//  rank of matrix P of partialled-out variables 
	return local dlist	`r(dlist)'			//  list of dropped collinear variables in P
	*

end


// copied from lassoutils 1.0.03 17/12/2017
// subroutine for fe transformation
program define _fe, rclass sortpreserve
	version 13
	syntax anything ,							/// anything is actually a varlist but this is more flexible - see below
		touse(varlist numeric max=1)			/// required `touse' variable (so marksample unneeded)
		tvarlist(string)						/// optional list of temp vars - may be "unfilled" (NaNs)
		FE(varlist numeric min=1 max=1) 		/// fe argument is ivar
		[										///
		NOFTOOLS								///
		]

	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename.
	local varlist `anything'
	*

	*** ftools
	// do not use ftools if (a) specifcally not requested; (b) not installed
	// use "which ftools" to check - faster, plus compile was already triggered
	// in conditional load section at end of this ado
	if "`noftools'"=="" {
		cap which ftools
		if _rc>0 {
			// fails check, likely not installed, so use (slower) Stata code
			local noftools "noftools"
		}
	}
	*
	
	*** Error check - if tvarlist provided, should have same number of elements as varlist
	local v_ct	: word count `varlist'
	local tv_ct	: word count `tvarlist'
	if `tv_ct' & (`v_ct'~=`tv_ct') {
		di as err "internal FE error - mismatched lists"
		exit 198
	}
	*

	if "`noftools'"~="" {
		// timer on 1
		*** Within-transformation / demeaning
		// varlist should be all doubles so no recast needed
		// first generate means
		// note that if data are xtset by panel and time var, this can change the order
		sort `fe' `touse'
		foreach var of varlist `varlist' {
			tempvar `var'_m
			qui by `fe' `touse' : gen double ``var'_m'=sum(`var')/_N if `touse'
			qui by `fe' `touse' : replace    ``var'_m'=``var'_m'[_N] if `touse' & _n<_N
		}
		// FE count; uses same sort
		tempvar T_i
		qui by `fe' `touse': gen long `T_i' = _N if `touse'
		qui by `fe' `touse': replace  `T_i' = .  if _n~=_N
		qui count if `T_i' < .
		local N_g	=r(N)
		// if data were xtset by panel and tim var, restore sort so that TS operators will work etc.
		cap xtset
		// now demean
		local i 1
		foreach var of varlist `varlist' {
			local tvar	: word `i' of `tvarlist'
			qui replace `tvar'=`var'-``var'_m' if `touse'
			local ++i
		}
		return scalar N_g = `N_g'
		*
		// timer off 1
	}
	else {
		// Mata routine; uses Sergio Correia's FTOOLS package.
		// timer on 2
		mata: s_fe("`varlist'", 				/// 
					"`tvarlist'",				///
					"`fe'",         		 	///
					"`touse'")
		local N_g	=r(N_g)
		return scalar N_g = `N_g'
		// timer off 2
	}
	return local noftools	`noftools'

end


// internal version of fvstrip 1.01 ms 24march2015
// takes varlist with possible FVs and strips out b/n/o notation
// returns results in r(varnames)
// optionally also omits omittable FVs
// expand calls fvexpand either on full varlist
// or (with onebyone option) on elements of varlist

program define fvstrip, rclass
	version 11.2
	syntax [anything] [if] , [ dropomit expand onebyone NOIsily ]
	if "`expand'"~="" {												//  force call to fvexpand
		if "`onebyone'"=="" {
			fvexpand `anything' `if'								//  single call to fvexpand
			local anything `r(varlist)'
		}
		else {
			foreach vn of local anything {
				fvexpand `vn' `if'									//  call fvexpand on items one-by-one
				local newlist	`newlist' `r(varlist)'
			}
			local anything	: list clean newlist
		}
	}
	foreach vn of local anything {									//  loop through varnames
		if "`dropomit'"~="" {										//  check & include only if
			_ms_parse_parts `vn'									//  not omitted (b. or o.)
			if ~`r(omit)' {
				local unstripped	`unstripped' `vn'				//  add to list only if not omitted
			}
		}
		else {														//  add varname to list even if
			local unstripped		`unstripped' `vn'				//  could be omitted (b. or o.)
		}
	}
// Now create list with b/n/o stripped out
	foreach vn of local unstripped {
		local svn ""											//  initialize
		_ms_parse_parts `vn'
		if "`r(type)'"=="variable" & "`r(op)'"=="" {			//  simplest case - no change
			local svn	`vn'
		}
		else if "`r(type)'"=="variable" & "`r(op)'"=="o" {		//  next simplest case - o.varname => varname
			local svn	`r(name)'
		}
		else if "`r(type)'"=="variable" {						//  has other operators so strip o but leave .
			local op	`r(op)'
			local op	: subinstr local op "o" "", all
			local svn	`op'.`r(name)'
		}
		else if "`r(type)'"=="factor" {							//  simple factor variable
			local op	`r(op)'
			local op	: subinstr local op "b" "", all
			local op	: subinstr local op "n" "", all
			local op	: subinstr local op "o" "", all
			local svn	`op'.`r(name)'							//  operator + . + varname
		}
		else if"`r(type)'"=="interaction" {						//  multiple variables
			forvalues i=1/`r(k_names)' {
				local op	`r(op`i')'
				local op	: subinstr local op "b" "", all
				local op	: subinstr local op "n" "", all
				local op	: subinstr local op "o" "", all
				local opv	`op'.`r(name`i')'					//  operator + . + varname
				if `i'==1 {
					local svn	`opv'
				}
				else {
					local svn	`svn'#`opv'
				}
			}
		}
		else if "`r(type)'"=="product" {
			di as err "fvstrip error - type=product for `vn'"
			exit 198
		}
		else if "`r(type)'"=="error" {
			di as err "fvstrip error - type=error for `vn'"
			exit 198
		}
		else {
			di as err "fvstrip error - unknown type for `vn'"
			exit 198
		}
		local stripped `stripped' `svn'
	}
	local stripped	: list retokenize stripped						//  clean any extra spaces
	
	if "`noisily'"~="" {											//  for debugging etc.
di as result "`stripped'"
	}

	return local varlist	`stripped'								//  return results in r(varlist)
end

// Internal version of matchnames
// Sample syntax:
// matchnames "`varlist'" "`list1'" "`list2'"
// takes list in `varlist', looks up in `list1', returns entries in `list2', called r(names)
program define matchnames, rclass
	version 11.2
	args	varnames namelist1 namelist2

	local k1 : word count `namelist1'
	local k2 : word count `namelist2'

	if `k1' ~= `k2' {
		di as err "namelist error"
		exit 198
	}
	foreach vn in `varnames' {
		local i : list posof `"`vn'"' in namelist1
		if `i' > 0 {
			local newname : word `i' of `namelist2'
		}
		else {
* Keep old name if not found in list
			local newname "`vn'"
		}
		local names "`names' `newname'"
	}
	local names	: list clean names
	return local names "`names'"
end


// Display varlist with specified indentation
program define DispVars
	version 11.2
	syntax [anything] [, _col(integer 15) _lc1(integer 0) _lc2(integer 0) ]
//	local maxlen = 80-`_col'
	local maxlen = c(linesize)-`_col'
	local len = 0
	local first = 1
	foreach vn in `anything' {
		local vnlen		: length local vn
		if `len'+`vnlen' > `maxlen' {
			di
			local first = 1
			local len = `vnlen'
			if `_lc1' {
				di as text _col(`_lc1') "{c |}" _c
			}
			if `_lc2' {
				di as text _col(`_lc2') "{c |}" _c
			}
		}
		else {
			local len = `len'+`vnlen'+1
		}
		if `first' {
			local first = 0
			di as res _col(`_col') "`vn'" _c
			}
		else {
			di as res " `vn'" _c
		}
	}
* Finish with a newline
	di
end

version 13
mata:

// makes "unfilled" tempvars in Mata (faster than initializing as missings)
void s_maketemps(real scalar p)
{
	(void) st_addvar("double", names=st_tempname(p), 1)
	st_global("r(varlist)",invtokens(names))
}

// basic IV parser
void s_ivparse(string scalar cmdline)

{

	// dep var is first word on line
	depvar	= ""
	i 		= 1
	c		= substr(cmdline,1,1)
	while (!((c==" ") | (c=="") | (c=="(") )) {
		depvar	= depvar + c
		i		= i + 1
		c		= substr(cmdline,i,1)
	}

	// loop through remainder, separating into dexog and parenthetical lists
	dexog	= ""
	dendog	= ""
	xctrl	= ""
	exexog	= ""
	pclauselist = J(0, 1, "")

	while (!(c=="")) {
		
		if (c!="(") {
			while (!((c=="(") | (c==""))) {
				dexog = dexog + c
				i	= i + 1
				c	= substr(cmdline,i,1)
			}
		}
		else {
			// start counter and initialize pclause
			pbal=1
			pclause = ""
			i	= i + 1
			c	= substr(cmdline,i,1)
			while (!((pbal==0) | (c==""))) {
				if (c=="(") {
					pbal=pbal+1
				}
				if (c==")") {
					pbal=pbal-1
				}
				if (pbal!=0) {
					pclause = pclause + c
					i	= i + 1
					c	= substr(cmdline,i,1)
				}
			}
			pclauselist = (pclauselist \ pclause)
			i	= i + 1
			c	= substr(cmdline,i,1)
			
			if (strpos(pclause,"=")) {
				epos	= strpos(pclause,"=")
				epos2	= strrpos(pclause,"=")
				if (epos==epos2) {
					dendog = dendog + substr(pclause, 1, epos-1)
					exexog = exexog + substr(pclause, epos+1, .)
				}
				else {
					errprintf("\nsyntax error - too many =s\n")
					exit(198)
				}
			}
			else {
				xctrl = xctrl + pclause + " "
			}
			
		}
	}

	depvar	= strtrim(stritrim(depvar))
	dexog	= strtrim(stritrim(dexog))
	dendog	= strtrim(stritrim(dendog))
	xctrl	= strtrim(stritrim(xctrl))
	exexog	= strtrim(stritrim(exexog))

	st_global("s(exexog)",exexog)
	st_global("s(xctrl)",xctrl)
	st_global("s(dendog)",dendog)
	st_global("s(dexog)",dexog)
	st_global("s(depvar)",depvar)

}

// copied from lassoutils 1.0.03 17/12/2017
// partial out program
void s_partial(	string scalar Ynames,
				string scalar Pnames,
				string scalar tYnames,
				string scalar touse,
				scalar cons,
				scalar solver)

{

// All varnames should be basic form, no FV or TS operators etc.
// Y = structural variables
// P = variables to be partialled out
// touse = sample
// cons = 0 or 1
// solver = 0, 1 or 2
// Strategy is to demean (numerically more stable in case of scaling problems)
// and then use svqrsolve(.) Mata program:
//   svsolve if no collinearites (more accurate),
//   qrsolve if there are collinearities (drops columns/set coeffs to zero),
//   and svsolve if qrsolve can't find the collinearities that svsolve can.

	Ytokens=tokens(Ynames)
	Ptokens=tokens(Pnames)
	st_view(Y, ., Ytokens, touse)
	st_view(P, ., Ptokens, touse)

	if (tYnames ~= "") {
		tYflag=1
		tYtokens=tokens(tYnames)
		st_view(tY, ., tYtokens, touse)
	}
	else {
		tYflag=0
	}

	L = cols(P)

	if (cons & L>0) {					//  Vars to partial out including constant
		Ymeans = mean(Y) 
		Pmeans = mean(P) 
	}
	else if (cons) {					//  Only constant to partial out = demean
		Ymeans = mean(Y) 
	}
	
//	Partial-out coeffs.
//	Not necessary if no vars other than constant. r=rank.
	if (cons & L>0) {
		b = svqrsolve(P :- Pmeans, Y :- Ymeans, r=., solver)	//  partial out P + cons
	}
	else if (!cons & L>0) {
		b = svqrsolve(P, Y, r=., solver)						//  partial out P
	}
	else {
		r=1														//  partial out cons
	}
//	Replace with residuals
	if (cons & L>0) {					//  Vars to partial out including constant
		if (tYflag) {
			tY[.,.] = (Y :- Ymeans) - (P :- Pmeans)*b
		}
		else {
			Y[.,.] = (Y :- Ymeans) - (P :- Pmeans)*b
		}
	}
	else if (!cons & L>0) {				//  Vars to partial out NOT including constant
		if (tYflag) {
			tY[.,.] = Y - P*b
		}
		else {
			Y[.,.] = Y - P*b
		}

	}
	else {								//  Only constant to partial out = demean
		if (tYflag) {
			tY[.,.] = (Y :- Ymeans)
		}
		else {
			Y[.,.] = (Y :- Ymeans)
		}
	}

	
// Return list of tokens of dropped collinear vars, if any
	if (r<cols(P)) {								//  something was dropped
		dsel = (b :== 0)							//  matrix, cols=#Y, rows=L, =1 if dropped
		dsel = (rowsum(dsel) :== cols(Y))			//  col vector, =1 if always dropped
		dlist = invtokens(select(Ptokens, dsel'))	//  string of names (tokens)
	}
	else {
		dlist = ""									//  empty list
	}
	st_global("r(dlist)",dlist)						//	Return list of dropped collinears.
	st_numscalar("r(rank)",r+cons)					//  Include constant in rank
}  
//end program s_partial


// copied from lassoutils 1.0.03 17/12/2017
// Mata utility for sequential use of SVD & QR solvers
// Default is SVD;
// if rank-deficient, use QR (drops columns);
// if QR then doesn't catch reduced rank, use original SVD.
// override: solver=0 => above;
//           solver=1 => always use SVD;
//           solver=2 => always use QR.
function svqrsolve (	numeric matrix A,			///
						numeric matrix B,			///
						|							///
						rank,						///
						real scalar solver			///
						)
{
	if (args()<=3) solver = 0

	real matrix Csv
	real matrix Cqr
	real scalar rsv
	real scalar rqr

	if (solver==0) {
		Csv = svsolve(A, B, rsv)
		if (rsv<cols(A)) {				//  not full rank, try QR
			Cqr = qrsolve(A, B, rqr)
			if (rsv<rqr) {				//  QR failed to detect reduced rank, use SVD after all
				C = Csv
			}
			else {						//  QR and SVD agree on reduced rank
				C = Cqr					//  QR dropped cols, use it
			}
		}
		else {
			C = Csv						//  full rank, use SVD (more accurate than QR)
		}
	}
	else if (solver==1) {
		C = svsolve(A, B, rsv)			//  override default, use SVD no matter what
	}
	else {
		C = qrsolve(A, B, rqr)			//  override default, use QR no mater whatt
	}
	if (solver<=1) {
		rank = rsv						//  rsv has rank
	}
	else {
		rank = rqr						//  forced use of QR so return QR rank
	}
	return(C)
}	// end svqrsolve


// based on code from weakiv
real matrix collapse_citable(									///
							string scalar citable,				///
							|									///
							real scalar level					//
							)
{

	if (args()<2) {
		level		= 0
	}

	m				= st_matrix(citable)
	// if p-values provided, convert to 1/0 rejection indicators
	if (level) {
		m[.,2]		= m[.,2] :<= level
	}
	smat1			= rowsum((m)[.,2..2])
	smat2			= smat1[(2::rows(smat1)),1]
	smat1			= smat1[(1::rows(smat1)-1),1]
	smat			= (smat1-smat2) :~= 0
	smat			= (1 \ smat) :| (smat \ 1)
	st_matrix("r(rtable)",select((m),smat))

}


// END MAIN MATA SECTION
end

*********** CONDITIONAL COMPILATION SECTION *****************
// Section is in Stata environment.

// FTOOLS
// Tell Stata to exit before trying to complile Mata function
// s_fe if required FTOOLS package is not installed.
// Note this runs only once, on loading, and if ftools has not
// been compiled, the check option will trigger compilation.

cap ftools, check
if _rc {
	// fails check, likely not installed, so do not complile
	// _fe Stata program will use (slower) Stata code
	exit
}

// Compile Mata function s_fe.
version 13
mata:


// FE transformation.
// Uses Sergio Correia's FTOOLS package - faster and does not require the data to be sorted.
void s_fe(		string scalar Xnames,
				string scalar tXnames,
				string scalar fe,
				string scalar touse)
{

	class Factor scalar F
	F = factor(fe, touse)
	F.panelsetup()

	Xtokens=tokens(Xnames)
	X = st_data( ., Xtokens, touse)
	tXtokens=tokens(tXnames)
	st_view(tX, ., tXtokens, touse)

	w = J(rows(X),1,1)
	counts = panelsum(w, F.info)
	means = panelsum(F.sort(X), w, F.info) :/ counts
	tX[.,.] = X - means[F.levels, .]

	N_g = F.num_levels
	st_numscalar("r(N_g)",N_g)

}

// End Mata section for s_fe
end

// END CONDITIONAL COMPILATION SECTION

******************** END ALL PROGRAM CODE *******************
exit
