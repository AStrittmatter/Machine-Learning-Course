*! lassoutils 1.0.08 5apr2018
*! authors aa/cbh/ms
* Adapted/expanded from lassoShooting version 12, cbh 25/09/2012
* mm_quantile from moremata package by Ben Jann:
*   version 1.0.8  20dec2007  Ben Jann

* Notes:
* Partialling out, temp vars, factor vars, FE transform all handled
* by calling programs.
* names_o and names_t (via options) are original and (possibly) temp names of Xs
* lassoutils branches to internal _rlasso, _lassopath, _fe, _unpartial, _partial, _std
* current cluster code is memory-intensive since it works with a full n x p matrix instead of cluster-by-cluster

* Updates (release date):
* 1.0.05  (30jan2018)
*         First public release.
*         Added seed(.) option to rlasso to control rnd # seed for xdep & sup-score.
*         Fixed up return code for lassoutils (method, alpha).
*         Promoted to required version 13 or higher.
*         Introduced centerpartial(.) Mata program for use with rlasso; returns X centered and with Xnp partialled-out.
*         Separate fields in datastruct: sdvec, sdvecpnp, sdvecp, sdvecnp. Latter two conformable with Xp and Xnp.
*         Added dots option for simulations (supscore, xdep).
*         Recoding relating to different treatment of cross-validation.
*         Changes to _std, _fe, _partial relating to holdout vs full sample.
*         Recoding of cons flag; now also dmflag to indcate zero-mean data.
* 1.0.06  (10feb2018)
*         Misc Mata speed tweaks following advice at http://scorreia.com/blog/2016/10/06/mata-tips.html,
*         e.g., evaluating for limits before loop; referring to vector elements with a single subscript.
*         Rewrite of FE transform to use Sergio Correia's FTOOLS package (if installed).
* 1.0.07  (17feb2018)
*         Bug fix related to ftools - leaves matastrict set to on after compilation, causing rest of lassoutils
*         to fail to load. Fix is to reset matastrict to what it was before calling ftools,compile.
* 1.0.08  (5apr2018)
*	      Changed the default maximum lambda for the elastic net: it was 2*max(abs(X'y)),
*		  and is now 2*max(abs(X'y))/max(0.001,alpha); see Friedman et al (2010, J of Stats Software). 
*		  Added Mata programs getInfoCriteria() and getMinIC(). getInfoCriteria() calculates information criteria
*		  along with RSS and R-squared. getMinIC() obtains minimum IC values and corresponding lambdas.
*		  Degrees of freedom calculation was added to DoLassoPath() and DoSqrtLassoPath().
* 		  Misc changes to the initial beta (= Ridge estimate), which in some cases didn't account for 
*		  penalty loadings.
* 		  Added lglmnet option to facilitate comparison with glmnet (was dysfunctional 'ladjust' option).
*         Restructured calc of lambda with cluster to match JBES paper; prev used sqrt(nclust) rather than sqrt(N),
*         now always 2c*sqrt(N)*invnormal(1-(gamma/log(N))/(2*p))). Fixed bug in calc of lambda for sqrt-lasso with cluster.
*         Definition of xdep lambda + cluster also changed slightly (by factor of sqrt(nclust/(nclust-1)).
*         Undocumented nclust1 option mimics behavior of CBH lassoCluster.ado; uses (nclust-1)*T rather than nclust*T.

program lassoutils, rclass sortpreserve
	version 13
	syntax [anything] ,							/// anything is actually a varlist but this is more flexible - see below
		[										/// 
		rlasso									/// branches to _rlasso
		path									/// branches to _lassopath
		unpartial								/// branches to _unpartial
		fe(string)								/// branches to _fe
		std										/// branches to _std
		partialflag(int 0)						/// branches to _partial if =1
		partial(string)							/// used by _partial and _unpartial
		tvarlist(string)						/// used by _partial, _fe, _std; optional list of temp vars to fill
												///
		ALPha(real 1) 							/// elastic net parameter
		SQRT									/// square-root-lasso
		ols 									///
		adaptive								///
												///
												/// verbose options will be transformed to 0/1/2 before passing to subs
		VERbose									/// shows additional detail
		VVERbose								/// even more detail
		*										///
		]

	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename.
	local varlist `anything'
	*

	** verbose
	if ("`vverbose'"!="") {
		local verbose=2							//  much detail
	} 
	else if ("`verbose'"!="") {
		local verbose=1							//  some additional detail
	}
	else {
		local verbose=0							//  no additional output
	}
	*

	************* SUBROUTINES: shoot, path, CV **************************************

	if "`rlasso'" ~= "" {
		_rlasso `varlist',							 	/// branch to _rlasso
			verbose(`verbose')							///
			`sqrt'										///
			`options'									//  no penloads etc.
	}
	else if "`path'" != "" {
		_lassopath `varlist',							///  branch to _ lassopath
			verbose(`verbose')							///
			`sqrt'										///
			alpha(`alpha')								///
 			`adaptive'									///
			`ols'										///
			`options'
	}
	else if "`unpartial'" ~= "" {
		_unpartial `varlist',							/// branch to _unpartial
			partial(`partial')							///
			`options'
	}
	else if `partialflag' {
		_partial `varlist',							 	/// branch to _partial
			partial(`partial')							///
			tvarlist(`tvarlist')						///
			`options'
	}
	else if "`fe'" ~= "" {
		_fe `varlist',								 	/// branch to _fe
			tvarlist(`tvarlist')						///
			fe(`fe')									///
			`options'
	}
	else if "`std'" ~= "" {
		_std `varlist',								 	/// branch to _std
			tvarlist(`tvarlist')						///
			`options'
	}
	else {
		di as err "internal lassoutils error"
		exit 499
	}

	*** return
	// so subroutine r(.) detail passed on
	// can add additional items to r(.) here
	return add

	// if estimation, return method etc.
	if "`rlasso'`path'`cv'" ~= "" {
		return scalar alpha			= `alpha'
		return scalar sqrt			= ("`sqrt'"!="")
		
		if ("`sqrt'" ~= "") {
			return local method		"sqrt-lasso"
		}
		else if `alpha'==1 {
			return local method		"lasso"
		}
		else if `alpha'==0 {
			return local method		"ridge"
		}
		else if (`alpha'<1) & (`alpha'>0)  {
			return local method		"elastic net"
		}
		else {
			di as err "internal lassoutils error - alpha outside allowed range"
			exit 499
		}
	}

end

// Subroutine for BCH _rlasso
program define _rlasso, rclass sortpreserve
	version 13
	syntax anything ,							/// anything is actually a varlist but this is more flexible - see below
		touse(varlist numeric max=1)			/// required `touse' variable (so marksample unneeded)
		[										///
		verbose(int 0)							///
												///
		ROBust									/// penalty level/loadings allow for heterosk. [homoskedasticity is the default]
		CLuster(varlist max=1)					/// penalty level/loadings allow for within-panel dependence & heterosk.
		XDEPendent								/// penalty level is estimated depending on X (default: No)
		numsim(integer 5000)					/// number of simulations with xdep (default=5,000)
												/// 
		TOLOpt(real 1e-10)						/// was originally 1e-5 but 1e-10 gives good equiv between xstd and stdloadings
		TOLUps(real 1e-4)						///
		TOLZero(real 1e-4)						///
		MAXIter(integer 10000)					/// max number of iterations in estimating lasso
		MAXUPSIter(int 2)						/// max number of lasso-based iterations in est penalty loadings
		CORRNumber(int 5) 						/// number of regressors used in InitialResiduals(); if =0, then y is used [as in CBH program]
												///
		LASSOUPS								/// use lasso residuals to estimate penalty loadings (post-lasso is default)
												///
		c(real 1.1)								/// "c" in lambda function
		gamma(real 0.1)							/// "gamma" numerator of fraction in lambda function gamma/log(N)
		gammad(real 0)							/// "gamma" denominator of fraction in lambda function gamma/log(N) (default = converted to log N)
		lambda0(real 0)							/// optional user-supplied lambda0
		LALTernative 							/// use alternative, less sharp lambda0 formula
		PMINus(int 0)							/// dim(X) adjustment in lambda0 formula
		nclust0									/// no longer in use as of 1.0.08; replaced by nclust1
		nclust1									/// use (nclust-1)*T instead of nclust*T in cluster-lasso
		center									/// center x_i*e_i or x_ij*e_ij in cluster-lasso
												///
		supscore								///
		ssnumsim(integer 500)					///
		testonly								///
												///
		seed(real -1)							/// set random # seed; relevant for xdep and supscore
		dots									/// display dots for xdep and supscore repetitions
												///
		xnames_o(string)						/// original names in varXmodel if tempvars used (string so varlist isn't processed)
		xnames_t(string)						/// temp names in varXmodel
		notpen_o(string)						/// used in checking
		notpen_t(string)						///
		consflag(int 0)							/// =0 if no cons or already partialled out
		dmflag(int 0)							/// data have been demeaned
												///
		stdy(string)							///
		stdx(string)							///
												///
		SQRT									/// square-root lasso
		]

	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename.
	local varlist `anything'
	*

	*** count number of obs
	_nobs `touse'
	local nobs = r(N)
	*

	*** sort needed if cluster option provided; sortpreserve means sort restored afterwards
	if "`cluster'"~="" {
		sort `cluster'
	}
	*

	*** define various parameters/locals/varlists
	local varY			`varlist'							// now is only element in varlist
	local varXmodel		`xnames_t'							// including notpen if any; may include tempnames
	local pen_t			: list xnames_t - notpen_t			// list of penalized regressors
	local p0			: word count `varXmodel' 			// #regressors incl notpen but excl constant
	local p				= `p0' - `pminus'					// p defined by BCH
	local p0			= `p0' + `consflag'					// #regressors but now incl constant
	*
	
	*** define various flags
	local clustflag		=("`cluster'"!="")					// =1 if cluster
	local hetero		=("`robust'"!="") & ~`clustflag'	// =1 if het-robust but NOT cluster
	local sqrtflag		=("`sqrt'"!="")
	local xdep			=("`xdependent'"!="")
	local lassoups		=("`lassoups'"!="")
	local lalternative	=("`lalternative'"!="")				// =1 if alternative, less sharp lambda0 should be used
	local nclust1		=("`nclust1'"!="")
	local center		=("`center'"!="")
	local supscoreflag	=("`supscore'`testonly'"!="")
	local testonlyflag	=("`testonly'"!="")
	local dotsflag		=("`dots'"!="")
	*

	*** syntax checks (incomplete)
	*

	tempname lambda slambda rmse rmseOLS		//  note lambda0 already defined as a macro
	tempname b bOLS sb sbOLS Ups sUps eUps stdvec
	tempname bAll bAllOLS
	tempname supscoremat BCH_ss BCH_p BCH_cv BCH_gamma hdm_ss hdm_p
	// initialize so that returns don't crash in case values not set
	local k				=0						//  used as flag to indicate betas are non-missing
	local N				=.
	local N_clust		=.
	local s				=.
	local s0			=.
	local niter			=.
	local nupsiter		=.
	scalar `lambda'		=.
	scalar `slambda'	=.
	scalar `rmse'		=.
	scalar `rmseOLS'	=.
	mat `b'				=.
	mat `bAll'			=.
	mat `bOLS'			=.
	mat `bAllOLS'		=.
	mat `sb'			=.
	mat `sbOLS'			=.
	mat `Ups'			=.
	mat `sUps'			=.
	mat `eUps'			=.
	mat `stdvec'		=.
	mat `supscoremat'	=.
	scalar `BCH_ss'		=.
	scalar `BCH_p'		=.
	scalar `BCH_cv'		=.
	scalar `BCH_gamma'	=.
	scalar `hdm_ss'		=.
	scalar `hdm_p'		=.

	if `p' & ~`testonlyflag' {					//  there are penalized model variables; estimate
	
		*** BCH rlasso
		mata:	EstimateRLasso(					///
							"`varY'",			///
							"`varXmodel'",		///
							"`xnames_o'",		///
							"`pen_t'",			///
							"`notpen_t'",		///
							"`touse'",			///
							"`cluster'",		///
							"`stdy'",			///
							"`stdx'",			///
							`sqrtflag',			///
							`hetero',			///
							`xdep',				///
							`numsim',			///
							`lassoups',			///
							`tolopt',			///
							`maxiter',			///
							`tolzero',			///
							`maxupsiter',		///
							`tolups',			///
							`verbose',			///
							`c',				///
							`gamma',			///
							`gammad',			///
							`lambda0',			///
							`lalternative',		///
							`corrnumber',		///
							`pminus',			///
							`nclust1',			///
							`center',			///
							`supscoreflag',		///
							`ssnumsim',			///
							`seed',				///
							`dotsflag',			///
							`consflag',			///
							`dmflag')

		*** Via ReturnResults(.)
		// coefs are returned as column vectors
		// convert to row vectors (Stata standard)
		mat `b'				=r(b)'					//  in original units
		mat `bOLS'			=r(bOLS)'
		mat `sb'			=r(sb)'					//  in standardized units
		mat `sbOLS'			=r(sbOLS)'
		mat `bAll'			=r(bAll)'
		mat `bAllOLS'		=r(bAllOLS)'
		mat `Ups'			=r(Ups)
		mat `sUps'			=r(sUps)
		mat `eUps'			=r(eUps)
		mat `stdvec'		=r(stdvecp)				//  stdvec for penalized vars only
		local selected0		`r(sel)'				//  selected variables INCL NOTPEN, EXCL CONSTANT
		local s0			=r(s)					//	number of selected vars INCL NOTPEN, EXCL CONSTANT; may be =0
		local k				=r(k)					//  number of all vars in estimated parameter vector INC CONSTANT; may be =0
		local niter			=r(niter)
		local nupsiter		=r(nupsiter)
		local N				=r(N)
		local N_clust		=r(N_clust)
		scalar `lambda'		=r(lambda)				//  relates to depvar in original units
		scalar `slambda'	=r(slambda)				//  relates to standardized depvar
		scalar `rmse'		=r(rmse)				//  lasso rmse
		scalar `rmseOLS'	=r(rmsePL)				//  post-lasso rmse
		if `supscoreflag' {
			mat `supscoremat'	=r(supscore)			//  sup-score row vector of results
			scalar `BCH_ss'		=`supscoremat'[1,colnumb(`supscoremat',"BCH_ss")]
			scalar `BCH_p'		=`supscoremat'[1,colnumb(`supscoremat',"BCH_p")]
			scalar `BCH_cv'		=`supscoremat'[1,colnumb(`supscoremat',"BCH_cv")]
			scalar `BCH_gamma'	=`supscoremat'[1,colnumb(`supscoremat',"BCH_gamma")]
			scalar `hdm_ss'		=`supscoremat'[1,colnumb(`supscoremat',"hdm_ss")]
			scalar `hdm_p'		=`supscoremat'[1,colnumb(`supscoremat',"hdm_p")]
		}
		// these may overwrite existing locals
		tempname lambda0 c gamma gammad
		scalar `lambda0'	=r(lambda0)				//  BCH definition of lambda; overwrites existing/default macro
		scalar `c'			=r(c)
		scalar `gamma'		=r(gamma)
		scalar `gammad'		=r(gammad)
		*
	}
	else if ~`testonlyflag' {						//  there are no penalized model vars; just do OLS

		if `p0' {									//  there are unpenalized vars and/or constant
			di as text "warning: no penalized variables; estimates are OLS"
		}
		if ~`consflag' {
			local noconstant noconstant
		}
		else {
			local consname _cons
		}
		qui regress `varY' `varXmodel' if `touse', `noconstant' cluster(`clustvar')
		foreach m in `b' `bOLS' `bAll' `bAllOLS' {
			mat `m'				=e(b)
			mat colnames `m'	=`xnames_o' `consname'
		}
		scalar `rmse'		=e(rmse)
		scalar `rmseOLS'	=e(rmse)
		local selected0		`xnames_o'
		local k				=e(df_m)+`consflag'		//  will be =0 if no variables in model at all
		local N				=e(N)
		if `clustflag' {
			local N_clust	=e(N_clust)
		}
		else {
			local N_clust	=0
		}
	}
	else if `testonlyflag' {					//  supscore test only

		mata:	EstimateSupScore(				///
							"`varY'",			///
							"`varXmodel'",		///
							"`xnames_o'",		///
							"`pen_t'",			///
							"`notpen_t'",		///
							"`touse'",			///
							"`cluster'",		///
							"`stdy'",			///
							"`stdx'",			///
							`hetero',			///
							`verbose',			///
							`ssnumsim',			///
							`c',				///
							`gamma',			///
							`pminus',			///
							`seed',				///
							`dotsflag',			///
							`consflag',			///
							`dmflag')

		mat  `supscoremat'	=r(supscore)
		scalar `BCH_ss'		=`supscoremat'[1,colnumb(`supscoremat',"BCH_ss")]
		scalar `BCH_p'		=`supscoremat'[1,colnumb(`supscoremat',"BCH_p")]
		scalar `BCH_cv'		=`supscoremat'[1,colnumb(`supscoremat',"BCH_cv")]
		scalar `BCH_gamma'	=`supscoremat'[1,colnumb(`supscoremat',"BCH_gamma")]
		scalar `hdm_ss'		=`supscoremat'[1,colnumb(`supscoremat',"hdm_ss")]
		scalar `hdm_p'		=`supscoremat'[1,colnumb(`supscoremat',"hdm_p")]
		local N				=r(N)
		local N_clust		=r(N_clust)

	}
	else {										//  shouldn't reach here
		di as err "internal _rlasso error"
		exit 499
	}

	*** check notpen is in selected0
	fvstrip `selected0'							//  fvstrip removes b/n/o prefixes (shouldn't need dropomit)
	local list1			`r(varlist)'
	fvstrip `notpen_o', dropomit				//  use dropomit option here
	local list2			`r(varlist)'
	local missing		: list list2 - list1
	local missing_ct	: word count `missing'
	if `missing_ct' {
		di as err "internal _rlasso error - unpenalized `missing' missing from selected vars"
		di as err "set tolzero(.) or other tolerances smaller or use partial(.) option"
		exit 499
	}
	*

	*** conventions
	// k			= number of selected/notpen INCLUDING constant; k>0 means estimated coef vector is non-empty
	// s0			= number of selected/notpen EXCLUDING constant
	// s			= number of selected (penalized)
	// p0			= number of all variables in model INCLUDING constant
	// selected0	= varlist of selected and unpenalized EXCLUDING constant; has s0 elements
	// selected		= varlist of selected; has s elements
	// cons			= 1 if constant, 0 if no constant
	// notpen_ct	= number of unpenalized variables in the model EXCLUDING CONSTANT
	// coef vectors b, bOLS, sb, sbOLS have k elements
	// coef vectors bAll, bAllOLS have p0 elements
	// Ups, sUps have p0-cons elements
	// stdvec has p0 - notpen_ct - cons elements (=#penalized)
	// Note that full set of regressors can including omitteds, factor base categories, etc.
	*

	*** fix colnames of beta vectors to include omitted "o." notation
	// includes post-lasso vector in case of collinearity => a selected variable was omitted
	// trick is to use _ms_findomitted utility but give it
	// diag(bAll) as vcv matrix where it looks for zeros
	// also build in fv info
	tempname tempvmat
	if `k' & ~`testonlyflag' {							//  if any vars in est coef vector (k can be zero)
		mat `tempvmat'	= diag(`bOLS')
		_ms_findomitted	`bOLS' `tempvmat'
		_ms_build_info	`bOLS' if `touse'
	}
	if `p0' & ~`testonlyflag' {							//  if any variables in model (p0 can be zero)
		mat `tempvmat'	= diag(`bAll')
		_ms_findomitted	`bAll' `tempvmat'
		_ms_build_info	`bAll' if `touse'
		mat `tempvmat'	= diag(`bAllOLS')
		_ms_findomitted	`bAllOLS' `tempvmat'
		_ms_build_info	`bAllOLS' if `touse'
	}
	*

	*** manipulate variable counts and lists
	// selected0 and s0 are selected+notpen (excluding constant)
	// selected and s are selected only
	local notpen_ct		: word count `notpen_o'			//  number of notpen EXCL CONSTANT
	local selected		: list selected0 - notpen_o		//  selected now has only selected/penalized variables
	local s				: word count `selected'			//  number of selected/penalized vars EXCL NOTPEN/CONSTANT
	*

	*** error checks
	// check col vector of b (selected) vs. k
	local col_ct		=colsof(`b')
	if colsof(`b')==1 & el(`b',1,1)==. & ~`testonlyflag' {
		// coef vector is empty so k should be zero
		if `k'~=0 {
			di as err "internal _rlasso error - r(k)=`k' does not match number of selected vars/coefs=`col_ct'"
			exit 499
		}
	}
	else if `k'~=`col_ct' & ~`testonlyflag' {
		// coef vector not empty so k should match col count
		di as err "internal _rlasso error - r(k)=`k' does not match number of selected vars/coefs=`col_ct'"
		exit 499
	}
	// check col vector of bAll vs. p0
	local col_ct		=colsof(`bAll')
	if `p0'~=`col_ct' & ~`testonlyflag' {
		// full coef vector not empty so p0 should match col count
		di as err "internal _rlasso error - p0=`p0' does not match number of model vars/coefs=`col_ct'"
		exit 499
	}
	*

	*** return results
	// coef vectors are row vectors (Stata standard)
	return matrix beta		=`b'
	return matrix betaOLS	=`bOLS'
	return matrix sbeta		=`sb'
	return matrix sbetaOLS	=`sbOLS'
	return matrix betaAll	=`bAll'
	return matrix betaAllOLS=`bAllOLS'
	return matrix Ups		=`Ups'
	return matrix sUps		=`sUps'
	return matrix eUps		=`eUps'
	return matrix stdvec	=`stdvec'					//  penalized vars only
	return scalar lambda	=`lambda'					//  relates to depvar in original units
	return scalar slambda	=`slambda'					//  relates to standardized depvar
	return scalar lambda0	=`lambda0'					//  BCH definition of lambda
	return scalar rmse		=`rmse'						//  lasso rmse
	return scalar c			=`c'
	return scalar gamma		=`gamma'
	return scalar gammad	=`gammad'
	return scalar rmseOLS	=`rmseOLS'					//  post-lasso rmse
	return local  selected0	`selected0'					//  all selected/notpen vars INCLUDING NOTPEN (but excl constant)
	return local  selected	`selected' 					//  all selected (penalized) vars EXCLUDING NOTPEN & CONS
	return scalar k			=`k'						//  number of all vars in sel/notpen parameter vector INCLUDING CONSTANT
	return scalar s0		=`s0'						//  number of vars selected INCLUDING NOTPEN (but excl constant)
	return scalar s			=`s'						//  number of vars selected EXCLUDING NOTPEN & CONS
	return scalar p0		=`p0'						//  number of all vars in original model including constant
	return scalar p			=`p'						//  p defined by BCH; excludes notpen & cons
	return scalar N_clust	=`N_clust'
	return scalar N			=`N'
	return scalar center	=`center'

	return local  clustvar	`cluster'
	return local  robust	`robust'
	return scalar niter		=`niter'
	return scalar maxiter	=`maxiter'
	return scalar nupsiter	=`nupsiter'
	return scalar maxupsiter=`maxupsiter'
	
	return scalar ssnumsim		=`ssnumsim'
	return scalar supscore		=`BCH_ss'
	return scalar supscore_p	=`BCH_p'
	return scalar supscore_cv	=`BCH_cv'
	return scalar supscore_gamma=`BCH_gamma'
	return scalar hdm_ss		=`hdm_ss'
	return scalar hdm_p			=`hdm_p'
	*

end		// end _rlasso subroutine

// Subroutine for lassopath
program define _lassopath, rclass sortpreserve
	version 13
	syntax [anything] ,							/// anything is actually a varlist but this is more flexible - see below
		toest(varlist numeric max=1)			/// required `touse' variable (so marksample unneeded)
		[										/// 
		verbose(int 0)							///
		consflag(int 1)							/// default is constant specified
		dmflag(int 0)							/// data have been demeaned
		notpen_t(string) 						///
		notpen_o(string) 						///
												///
												/// lambda settings
		ALPha(real 1) 							/// elastic net parameter
		Lambda(string)							/// overwrite default lambda
		LCount(integer 100)						///
		LMAX(real 0) 							///
		LMINRatio(real 1e-4)					/// ratio of maximum to minimum lambda
												///
		TOLOpt(real 1e-10)						/// was originally 1e-5 but 1e-10 gives good equiv between xstd and stdloadings
		TOLZero(real 1e-4)						///
		MAXIter(integer 10000)					/// max number of iterations in estimating lasso
												///
		PLoadings(string)						/// set penalty loadings as Stata row vector
												///
		LGLMnet									/// lambda is divided by 2*n (to make results comparable with glmnet)
												///
		xnames_o(string)						/// original names in varXmodel if tempvars used (string so varlist isn't processed)
		xnames_t(string)						/// temp names in varXmodel
												///
		sqrt 									/// square-root lasso
		ols										///
												///
  		STDLflag(int 0) 						/// use standardisation loadings?
		STDCOEF(int 0) 							/// don't unstandardize coefs
		stdy(string)							///
		stdx(string)							///
												///
		ADAPTive  								/// adaptive lasso
		ADATheta(real 1) 						/// gamma paramater for adapLASSO
		ADALoadings(string)						///
		holdout(varlist) 						///
												///
		NOIC									///
		]
		
	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename to `varlist' (usual Stata convention).
	local varlist `anything'

	*** count number of obs
	_nobs `toest'
	local nobs = r(N)
	*

	*** define various parameters/locals/varlists
	local varY			`varlist'						// now is only element in varlist
	local varXmodel		`xnames_t'						// including notpen if any; may include tempnames
	local pen_t			: list xnames_t - notpen_t		// list of penalized regressors
	local pmodel		: word count `varXmodel' 		//  # regressors (excl. constant) minus partial
	local p0			: word count `varXmodel' 		// # all regressors incl notpen but excl constant
	local p0			= `p0' + `consflag'					// ditto but now incl constant
	*

	*** syntax checks
	if (`lminratio'<0) | (`lminratio'>=1) {
		di as err "lminratio should be greater than 0 and smaller than 1."
	}	
	if ("`lambda'"!="") & ((`lcount'!=100) | (`lmax'>0) | (`lminratio'!=1e-4)) {
		di as err "lcount | lmax | lminratio option ignored since lambda() is specified."
	}
	*
	
	*** useful flags and counts
	local sqrtflag		=("`sqrt'"!="")
	local olsflag		=("`ols'"!="")
	local adaflag 		=("`adaptive'"!="")
	local lglmnetflag	=("`lglmnet'"!="")
	local noicflag		=("`noic'"!="")
	if (~`consflag') {
		local noconstant noconstant
	}
	*
	
	*********************************************************************************
	*** penalty loadings and adaptive lasso 									  ***
	*********************************************************************************
	
	*** syntax checks for penalty loadings
	if (`adaflag') & ("`ploadings'"!="") {
		di as error "ploadings() option not allowed with adaptive."
		exit 198
	}
	if (~`adaflag') & ("`adaloadings'"!="") {
		di as error "adaloadings(`adaloadings') ignored. specify adaptive option."
	}
	if (~`adaflag') & (`adatheta'!=1) {
		di as err "adatheta(`adatheta') ignored. specify adaptive option."
	}
	local pltest = ("`ploadings'"!="")+`adaflag'+`stdlflag'
	if (`pltest'>1) {
		di as error "only one of the following options allowed: ploadings, adaptive, stdloadings."
		exit 198
	}
	** check dimension of penalty loading vector
	if ("`ploadings'"!="") {
		// Check that ploadings is indeed a matrix
		tempname Ups
		cap mat `Ups' = `ploadings'
		if _rc != 0 {
			di as err "invalid matrix `ploadings' in ploadings option"
			exit _rc
		}
		// check dimension of col vector
		if (colsof(`Ups')!=(`pmodel')) | (rowsof(`Ups')!=1) {
			di as err "`ploadings' should be a row vector of length `pmodel'=dim(X) (excl constant)."
			exit 503
		}
	}
	else if (`adaflag') {
		if ("`adaloadings'"=="") {
			// obtain ols coefficients
			sum `toest', meanonly
			local nobs = r(N)
			tempname adaptiveUps
			if (`pmodel'>`nobs') { // is dim(X)> n ? if yes, do univariate OLS
					di as text "Adaptive weights calculated using univariate OLS regressions since dim(X)>#Observation."
					mat `adaptiveUps' = J(1,`pmodel',0)
					local ip=1
					foreach var of varlist `varXmodel' {
						qui _regress `varY' `var' if `toest', `noconstant'
						mat `adaptiveUps'[1,`ip'] = _b[`var']
						local ip=`ip'+1
				}
			}
			else { // is dim(X)>n ? if not, do full OLS
				di as text "Adaptive weights calculated using OLS."
				_regress `varY' `varXmodel' if `toest', `noconstant'
				mat `adaptiveUps' = e(b)
				if (`consflag') {
					mat `adaptiveUps' = `adaptiveUps'[1,1..`pmodel']
				}
			}
		} 
		else {
			tempname adaptiveUps0 adaptiveUps
			cap mat `adaptiveUps0' = `adaloadings'
			if _rc != 0 {
				di as err "invalid matrix `adaloadings' in adaloadings option"
				exit _rc
			}
			local ebnames: colnames `adaptiveUps0'
			local ebnamescheck: list varXmodel - ebnames
			if ("`ebnamescheck'"!="") {
				di as err "`ebnamescheck' do not appear in matrix `adaloadings' as colnames."
				exit 198
			}
			mat `adaptiveUps' = J(1,`pmodel',.)
			matrix colnames `adaptiveUps0' = varXmodel
			local j = 1
			foreach var of local varXmodel {
				local vi : list posof "`var'" in ebnames
				mat `adaptiveUps'[1,`j']=`adaptiveUps0'[1,`vi']
				local j=`j'+1
			}
		}
		// do inversion
		forvalue j = 1(1)`pmodel' {
			mat `adaptiveUps'[1,`j'] = abs((1/`adaptiveUps'[1,`j'])^`adatheta')
		}
		tempname Ups
		mat `Ups' = `adaptiveUps'
	}
	*

	*********************************************************************************
	*** penalty loadings and adaptive lasso END									  ***
	*********************************************************************************
	
	//egen test = rowmiss(`varXmodel') if `toest'
	//sum test

	mata:	EstimateLassoPath(			///
					"`varY'",			///
					"`varXmodel'",		///
					"`xnames_o'",		///
					"`notpen_o'",		///
					"`notpen_t'",		///
					"`toest'",			///
					"`holdout'", 		/// holdout for MSPE
					`consflag',			///
					`dmflag',			///
					"`lambda'",			/// lambda matrix or missing (=> construct default list)
					`lmax',				/// maximum lambda (optional; only used if lambda not specified)
					`lcount',			/// number of lambdas (optional; only used if lambda not specified)
					`lminratio',		/// lmin/lmax ratio (optional; only used if lambda not specified)
					`lglmnetflag',		/// 
					"`Ups'",			///
					"`stdy'",			/// Stata matrix with SD of dep var
					"`stdx'",			/// Stata matrix with SDs of Xs
					`stdlflag',			/// use standardisation loadings  
					`sqrtflag',			/// sqrt lasso  
					`alpha',			/// elastic net parameter  
					`olsflag',			/// post-OLS estimation
					`tolopt',			///
					`maxiter',			///
					`tolzero',			///
					`verbose',			///
					`stdcoef',			///
					`noicflag'			///
					)

	
	if (`r(lcount)'>1) { //------- #lambda > 1 -----------------------------------------------//
	
		tempname Ups betas dof lambdamat l1norm wl1norm stdvec shat
		mat `Ups' 			= r(Ups)
		mat `betas' 		= r(betas)
		mat `dof'			= r(dof)
		mat `shat'			= r(shat)
		mat `lambdamat'		= r(lambdalist)
		mat `l1norm'		= r(l1norm)
		mat `wl1norm'		= r(wl1norm)
		mat `stdvec'		= r(stdvec)
		if ("`holdout'"!="") {
			tempname mspe0
			mat `mspe0' = r(mspe)	
		}
		else {
			tempname rss ess tss rsq aic aicc bic ebic
			mat `rss' = r(rss)
			mat `ess' = r(ess)
			mat `tss' = r(tss)
			mat `rsq' = r(rsq)
			mat `aic' = r(aic)
			mat `bic' = r(bic)
			mat `aicc' = r(aicc)
			mat `ebic' = r(ebic)
			return scalar aicmin = r(aicmin)
			return scalar laicid = r(laicid)
			return scalar aiccmin = r(aiccmin)
			return scalar laiccid = r(laiccid)
			return scalar bicmin = r(bicmin)
			return scalar lbicid = r(lbicid)
			return scalar ebicmin = r(ebicmin)
			return scalar lebicid = r(lebicid)
		}
		
		return scalar N				= r(N)
		return scalar lcount 		= r(lcount)
		return scalar olsflag		= `olsflag'
		return scalar lmax 			= r(lmax)
		return scalar lmin 			= r(lmin)
		return matrix betas 		= `betas'
		return matrix dof 			= `dof'
		return matrix shat 			= `shat'
		return matrix lambdalist 	= `lambdamat'
		return matrix Ups			= `Ups'
		return matrix l1norm		= `l1norm'
		return matrix wl1norm		= `wl1norm'
		return matrix stdvec 		= `stdvec'
		if ("`holdout'"!="") {
			return matrix mspe 		= `mspe0'
		}
		else {
			return matrix rss 		= `rss' 
			return matrix ess 		= `ess' 
			return matrix tss		= `tss' 
			return matrix rsq 		= `rsq' 
			return matrix aic 		= `aic' 
			return matrix aicc 		= `aicc' 
			return matrix bic 		= `bic' 
			return matrix ebic 		= `ebic' 
		}
	}
	else if (`r(lcount)'==1) { 
	
		*** Via ReturnResults(.)
		*** the following code is based on _rlasso
		tempname b bOLS sb sbOLS Ups sUps stdvec
		tempname bAll bAllOLS
		tempname lambda slambda lambda0 rmse rmseOLS
		// coefs are returned as column vectors
		// convert to row vectors (Stata standard)
		mat `b'				=r(b)'					//  in original units
		mat `bOLS'			=r(bOLS)'
		//*//mat `sb'			=r(sb)'					//  in standardized units
		//*//mat `sbOLS'			=r(sbOLS)'
		mat `bAll'			=r(bAll)'
		mat `bAllOLS'		=r(bAllOLS)'
		mat `Ups'			=r(Ups)
		//*//mat `sUps'			=r(sUps)
		mat `stdvec'		=r(stdvec)
		local selected0		`r(sel)'				//  selected variables INCL NOTPEN, EXCL CONSTANT
		local s0			=r(s)					//	number of selected vars INCL NOTPEN, EXCL CONSTANT; may be =0
		local k				=r(k)					//  number of all vars in estimated parameter vector INC CONSTANT; may be =0
		local niter			=r(niter)
		//*//local nupsiter		=r(nupsiter)
		local N				=r(N)
		local N_clust		=r(N_clust)
		scalar `lambda'		=r(lambda)				//  relates to depvar in original units
		//*//scalar `slambda'	=r(slambda)				//  relates to standardized depvar
		//*//scalar `lambda0'	=r(lambda0)				//  BCH definition of lambda
		scalar `rmse'		=r(rmse)				//  lasso rmse
		scalar `rmseOLS'	=r(rmsePL)				//  post-lasso rmse
		*

		*** check notpen is in selected0
		fvstrip `selected0'							// fvstrip removes b/n/o prefixes.
		local list1			`r(varlist)'
		fvstrip `notpen_o', dropomit				// use dropomit option here
		local list2		`r(varlist)'
		local missing		: list list2 - list1
		local missing_ct	: word count `missing'
		if `missing_ct' {
			di as err "internal _rlasso error - unpenalized `missing' missing from selected vars"
			di as err "set tolzero(.) or other tolerances smaller or use partial(.) option"
			exit 499
		}
		*
		
		*** conventions
		// k			= number of selected/notpen INCLUDING constant; k>0 means estimated coef vector is non-empty
		// s0			= number of selected/notpen EXCLUDING constant
		// s			= number of selected (penalized)
		// p0			= number of all variables in model INCLUDING constant
		// selected0	= varlist of selected and unpenalized EXCLUDING constant; has s0 elements
		// selected		= varlist of selected; has s elements
		// cons			= 1 if constant, 0 if no constant
		// notpen_ct	= number of unpenalized variables in the model EXCLUDING CONSTANT
		// coef vectors b, bOLS, sb, sbOLS have k elements
		// coef vectors bAll, bAllOLS have p0 elements
		// Ups, sUps have p0-cons elements
		// stdvec has p0 - notpen_ct - cons elements (=#penalized)
		// Note that full set of regressors can including omitteds, factor base categories, etc.
		*

		*** fix colnames of beta vectors to include omitted "o." notation
		// includes post-lasso vector in case of collinearity => a selected variable was omitted
		// trick is to use _ms_findomitted utility but give it
		// diag(bAll) as vcv matrix where it looks for zeros
		// also build in fv info
		tempname tempvmat
		mat `tempvmat'	= diag(`bOLS')
		_ms_findomitted	`bOLS' `tempvmat'
		_ms_build_info	`bOLS' if `toest'
		mat `tempvmat'	= diag(`bAll')
		_ms_findomitted	`bAll' `tempvmat'
		_ms_build_info	`bAll' if `toest'
		mat `tempvmat'	= diag(`bAllOLS')
		_ms_findomitted	`bAllOLS' `tempvmat'
		_ms_build_info	`bAllOLS' if `toest'
		*

		*** manipulate variable counts and lists
		// selected0 and s0 are selected+notpen (excluding constant)
		// selected and s are selected only
		local notpen_ct		: word count `notpen_o'			//  number of notpen EXCL CONSTANT
		local selected		: list selected0 - notpen_o		//  selected now has only selected/penalized variables
		local s				: word count `selected'			//  number of selected/penalized vars EXCL NOTPEN/CONSTANT
		*

		*** error checks
		// check col vector of b (selected) vs. k
		local col_ct		=colsof(`b')
		if colsof(`b')==1 & el(`b',1,1)==. {	//  coef vector is empty so k should be zero
			if `k'~=0 {
				di as err "internal _rlasso error - r(k)=`k' does not match number of selected vars/coefs=`col_ct'"
				exit 499
			}
		}
		else if `k'~=`col_ct' {					//  coef vector not empty so k should match col count
			di as err "internal _rlasso error - r(k)=`k' does not match number of selected vars/coefs=`col_ct'"
			exit 499
		}
		// check col vector of bAll vs. p0
		local col_ct		=colsof(`bAll')
		if `p0'~=`col_ct' {						// full coef vector not empty so k should match col count
			di as err "internal _rlasso error - p0=`p0' does not match number of model vars/coefs=`col_ct'"
			exit 499
		}
		*

		*** return results
		// coef vectors are row vectors (Stata standard)
		return matrix beta		=`b'
		return matrix betaOLS	=`bOLS'
		//*//return matrix sbeta		=`sb'
		//*//return matrix sbetaOLS	=`sbOLS'
		return matrix betaAll	=`bAll'
		return matrix betaAllOLS=`bAllOLS'
		return matrix Ups		=`Ups'
		//*//return matrix sUps		=`sUps'
		return matrix stdvec	=`stdvec'					//  penalized vars only
		return scalar lambda	=`lambda'					//  relates to depvar in original units
		//*//return scalar slambda	=`slambda'					//  relates to standardized depvar
		//*//return scalar lambda0	=`lambda0'					//  BCH definition of lambda
		return scalar rmse		=`rmse'						//  lasso rmse
		return scalar rmseOLS	=`rmseOLS'					//  post-lasso rmse
		return local  selected0	`selected0'					//  all selected/notpen vars INCLUDING NOTPEN (but excl constant)
		return local  selected	`selected' 					//  all selected (penalized) vars EXCLUDING NOTPEN & CONS
		return scalar k			=`k'						//  number of all vars in sel/notpen parameter vector INCLUDING CONSTANT
		return scalar s0		=`s0'						//  number of vars selected INCLUDING NOTPEN (but excl constant)
		return scalar s			=`s'						//  number of vars selected EXCLUDING NOTPEN & CONS
		return scalar p0		=`p0'						//  number of all vars in original model including constant
		return scalar N_clust	=`N_clust'
		return scalar N			=`N'
		//*//return scalar center	=`center'
		return local  clustvar	`cluster'
		return local  robust	`robust'
		return scalar niter		=`niter'
		return scalar maxiter	=`maxiter'
		//*//return scalar nupsiter	=`nupsiter'
		//*//return scalar maxupsiter=`maxupsiter'
		return scalar lcount 	= 1
		return scalar olsflag	= `olsflag'
	}
	else {
		di as err "internal _lassopath error - lcount=`lcount'"
		exit 499
	}
end

// subroutine for partialling out
program define _partial, rclass sortpreserve
	version 13
	syntax anything ,							/// anything is actually a varlist but this is more flexible - see below
		touse(varlist numeric max=1)			/// required `touse' variable (so marksample unneeded)
		[										/// 
		toest(varlist numeric max=1)			/// optional `toest' variable (subset on which standardization is based)
		PARtial(string)							/// string so that fv operators aren't inserted
		tvarlist(string)						/// optional list of temp vars - may be "unfilled" (NaNs)
		dmflag(int 0)							/// =0 if cons and not demeaned; =1 if nocons and already demeaned
		solver(string)							/// svd, qr or empty
		]

	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename.
	local varlist `anything'
	*

	*** toest vs touse
	// touse = all data
	// toest = estimation sample
	// if toest is missing, set equal to touse
	if "`toest'"=="" {
		tempvar toest
		qui gen byte `toest' = `touse'
	}
	*

	*** Error check - if tvarlist provided, should have same number of elements as varlist
	local v_ct	: word count `varlist'
	local tv_ct	: word count `tvarlist'
	if `tv_ct' & (`v_ct'~=`tv_ct') {
		di as err "internal lassoutils partialling error - mismatched lists"
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
					"`toest'",				/// toest
					`dmflag', 				/// treatment of constant/demeaned or not
					`solver')				//  choice of solver (optional)
	return scalar rank	=`r(rank)'			//  rank of matrix P of partialled-out variables 
	return local dlist	`r(dlist)'			//  list of dropped collinear variables in P
	*

end


// subroutine for fe transformation
program define _fe, rclass sortpreserve
	version 13
	syntax anything ,							/// anything is actually a varlist but this is more flexible - see below
		touse(varlist numeric max=1)			/// required `touse' variable (so marksample unneeded)
		tvarlist(string)						/// optional list of temp vars - may be "unfilled" (NaNs)
		FE(varlist numeric min=1 max=1) 		/// fe argument is ivar
		[										///
		toest(varlist numeric max=1)			/// optional `toest' variable (subset on which standardization is based)
		NOFTOOLS								///
		]

	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename.
	local varlist `anything'
	*

	*** toest vs touse
	// toest = estimation sample
	// touse = all data
	// if toest is missing, set equal to touse
	if ("`toest'"=="") {
		tempvar toest
		qui gen `toest'=`touse'
	}
	*

	*** ftools
	// do not use ftools if (a) specifcally not requested; (b) not installed
	// use "which ftools" to check - faster, plus compile was already triggered
	// in conditional load section at end of this ado
	if "`noftools'"=="" {
		cap which ftools
		if _rc {
			// fails check, likely not installed, so use (slower) Stata code
			local noftools "noftools"
		}
	}
	*
	
	*** Error check - if tvarlist provided, should have same number of elements as varlist
	local v_ct	: word count `varlist'
	local tv_ct	: word count `tvarlist'
	if (`v_ct'~=`tv_ct') {
		di as err "internal lassoutils FE error - mismatched lists"
		exit 198
	}
	*

	if "`noftools'"~="" {
		// timer on 1
		*** Within-transformation / demeaning
		// varlist should be all doubles so no recast needed
		// xtset is required for FEs so this check should never fail
		cap xtset
		if _rc {
			di as err "internal lassoutils xtset error"
			exit 499
		}
		// panelvar always exists; timevar may be empty
		// if data xtset by panelvar only, may not be sorted
		// if data xtset by panelvar and time var, will be sorted by both
		// resorting can be triggers by simple call to xtset
		local panelvar `r(panelvar)'
		local timevar `r(timevar)'
		// sort by panel and estimation sample; put latest observation last
		sort `panelvar' `toest' `timevar'
		// toest_m is indicator that this ob will have the mean
		// N_est is number of obs in panel and in estimation sample
		tempvar toest_m N_est
		// last ob in panel and estimation sample tagged to have mean
		qui by `panelvar' `toest' : gen `toest_m' = (_n==_N) & `toest'
		qui by `panelvar' `toest' : gen `N_est' = sum(`toest')
		// count is of panels used in estimation
		qui count if `toest_m'
		local N_g	=r(N)
		// if xtset by panel and time vars, restore sort
		cap xtset
		// create means for each variable
		foreach var of varlist `varlist' {
			tempvar `var'_m
			local mlist `mlist' ``var'_m'
			// use only training/estimation data to calculate mean (toest)
			qui by `panelvar' : gen double ``var'_m'=sum(`var')/`N_est' if `toest'
			qui by `panelvar' : replace ``var'_m' = . if ~`toest_m'
		}
		// sort so that last ob in each panel has the mean in it
		sort `panelvar' `toest_m'
		// and propagate to the rest of the panel
		foreach var of varlist `mlist' {
			qui by `panelvar' : replace `var' = `var'[_N] if _n<_N
		}
		// if xtset by panel and time vars, restore sort
		// need to do this if e.g. any time-series operators in use
		cap xtset
		// finally, demean data
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
					"`touse'",					///
					"`toest'")
		local N_g	=r(N_g)
		return scalar N_g = `N_g'
		// timer off 2
	}
	// indicate whether FTOOLS not used
	return local noftools `noftools'
end


// subroutine for standardizing in Stata
program define _std, rclass
	version 13
	syntax anything ,							/// anything is actually a varlist but this is more flexible - see below
		touse(varlist numeric max=1)			/// required `touse' variable (so marksample unneeded)
		[										///
		toest(varlist numeric max=1)			/// optional `toest' variable (subset on which standardization is based)
		tvarlist(string)						/// optional list of temp vars - may be "unfilled" (NaNs)
		consmodel(int 1)						/// =1 if constant in model (=> data are to be demeaned)
		dmflag(int 0)							/// =1 if data already demeaned
		NOChange ///
		]

	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename.
	local varlist `anything'
	*
	
	local transform = ("`nochange'"=="")
	
	*** toest vs touse
	// touse = all data
	// toest = estimation sample
	// if toest is missing, set equal to touse
	if "`toest'"=="" {
		tempvar toest
		qui gen byte `toest' = `touse'
	}
	*

	tempname stdvec mvec
	mata: s_standard("`varlist'","`tvarlist'","`touse'","`toest'",`consmodel',`dmflag',`transform') 
	mat `stdvec' = r(stdvec)
	mat `mvec' = r(mvec)
	// these will be tempnames if std program called with tempvars
	mat colnames `stdvec' = `varlist'
	mat colnames `mvec' = `varlist'
	return matrix stdvec = `stdvec'
	return matrix mvec = `mvec'	
end


// subroutine for recovering coefs of partialled-out vars
program define _unpartial, rclass sortpreserve
	version 13
	syntax anything ,							/// anything is actually a varlist but this is more flexible - see below
		touse(varlist numeric max=1)			/// required `touse' variable (so marksample unneeded)
		[										///
		beta(string)							///
		depvar(string)							///
		scorevars(string)						///
		names_t(string)							/// string so that fv operators aren't inserted
		names_o(string)							/// ditto
		PARtial(string)							/// ditto
		consmodel(int 1)						/// include constant when recovering coefs
		]

	*** `anything' is actually a varlist but this is more flexible.
	// avoids problems with uninitialized vars, unwanted parsing of factor vars, etc.
	// ... so immediately rename.
	local varlist `anything'
	*

	local depvar `varlist'
	
	tempname b
	tempvar xb yminus
	if "`scorevars'" ~= "" {
		mat `b' = `beta'
		mat colnames `b' = `scorevars'
		qui mat score double `xb' = `b' if `touse'
		qui gen double `yminus' = `depvar' - `xb' if `touse'
	}
	else {
		qui gen double `yminus' = `depvar'
	}
	//  partial uses _tnames; use _t names and then convert to _o names
	//  _t names will be FE-transformed or just the values of the original vars (without i/b/n etc.)
	if ~`consmodel' {
		local noconstant noconstant
	}
	qui reg `yminus' `partial' if `touse', `noconstant'
	tempname bpartial
	mat `bpartial' = e(b)
	// replace temp names with original names
	local cnames_t	: colnames `bpartial'				//  may have _cons at the end already
	fvstrip `cnames_t'									//  regress output may have omitted vars
	local cnames_t	`r(varlist)'
	matchnames "`cnames_t'" "`names_t'" "`names_o'"
	local cnames_o	`r(names)'
	mat colnames `bpartial' = `cnames_o'
	// may be omitteds so need to reinsert o. notation
	// trick is to use _ms_findomitted utility but give it
	// diag(bAll) as vcv matrix where it looks for zeros
	tempname tempvmat
	mat `tempvmat'	= diag(`bpartial')
	_ms_findomitted `bpartial' `tempvmat'
	// build in fv info
	_ms_build_info `bpartial' if `touse'
	// attach to b matrix if b not empty
	if `beta'[1,1] ~= . {
		mat `b' = `beta' , `bpartial'
	}
	else {
		mat `b' = `bpartial'
	}
	
	return matrix b = `b'
	return matrix bpartial = `bpartial'

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

********************************************************************************
*** Mata section															 ***
********************************************************************************

version 13
mata:

// data structure
struct dataStruct {
	pointer colvector y
	pointer matrix X
	pointer colvector clustid	// cluster id
	string colvector nameX		// names of actual variables (can be tempvars)
	string colvector nameX_o	// original names of variables
	real scalar cons			// =1 if model also has a constant, =0 otherwise
	real scalar dmflag			// =1 if data have mean zero, =0 otherwise
	real scalar n				// number of observations
	real scalar nclust			// number of clusters; 0 if no clustering
	real scalar p				// number of columns in X (may included constant)
	real rowvector sdvec		// vector of SDs of X
	real rowvector varvec		// vector of variances of X
	real scalar ysd				// standard deviation of y
	real scalar prestdflag		// flag indicating data have been pre-standardized
	real rowvector mvec			// vector of means
	real scalar ymvec			// mean of y (not actually a vector)
	real matrix pihat			// used for partialling out with rlasso
	real matrix ypihat			// used for partialling out with rlasso
	real matrix XX				// cross prod of all Xs
	real matrix Xy				// cross prod of all Xs and y
	pointer matrix Xp			// used by rlasso - penalized Xs
	pointer matrix Xnp			// used by rlasso - unpenalized Xs
	string colvector nameXp		// used by rlasso - names of actual variables (can be tempvars)
	string colvector nameXnp	// used by rlasso - names of actual variables (can be tempvars)
	real rowvector selXp		// used by rlasso - selection row vector for penalized vars
	real rowvector selXnp		// used by rlasso - selection row vector for unpenalized vars
	real rowvector selindXp		// used by rlasso - selection index row vector for penalized vars
	real rowvector selindXnp	// used by rlasso - selection index row vector for unpenalized vars
	real rowvector mvecp		// used by rlasso - vector of means of penalized Xs
	real rowvector mvecnp		// used by rlasso - vector of means of unpenalized Xs
	real rowvector sdvecp		// used by rlasso - sdvec of penalized vars after partialling-out unpenalized vars
	real rowvector sdvecpnp		// used by rlasso - same as sdvecp except conformable with full set of Xs (p and np)
	real scalar ysdp			// used by rlasso - SD of y after partialling out unpenalized vars
	}

struct outputStruct {
	real colvector beta			// beta in original units
	//real colvector sbeta		// beta in standardized units
	real colvector betaPL		// OLS (post-lasso) beta in std units 
	//real colvector sbetaPL
	real colvector betaAll		// full beta vector
	real colvector betaAllPL	// full OLS (post-lasso) beta vector
	real colvector beta_init
	real scalar cons			// flag =1 if cons in regression, =0 if no cons
	real scalar intercept		// estimated intercept
	real scalar interceptPL		// estimated OLS (post-lasso) intercept
	string colvector nameXSel	// vector of names of selected Xs
	real colvector index		// index of selected vars
	real rowvector Ups			// penalty loadings
	real rowvector sUps			// standardized penalty loadings
	real rowvector eUps			// estimated penalty loadings (rlasso only)
	real scalar prestdflag		// flag indicating data have been pre-standardized
	real colvector v			// residuals
	real colvector vPL			// OLS (post-lasso) residuals
	real scalar lambda			// penalty scalar
	real scalar slambda			// standardized
	real scalar lambda0			// lambda without estimate of sigma (rlasso only)
	real scalar c				// part of BCH lambda
	real scalar gamma			// part of BCH lambda
	real scalar gammad			// part of BCH lambda
	real scalar rmse			// rmse using estimated beta
	real scalar rmsePL			// rmse using OLS (post-lasso) beta
	real scalar n				// number of obs
	real scalar s				// number of selected vars
	real scalar nclust			// number of clusters (rlasso only)
	real scalar niter			// number of iterations
	real scalar nupsiter		// number of rlasso penalty loadings iterations
	real rowvector supscore		// sup-score stats: BCH stat, p-value, crit-value, signif; rlasso stat, p-value
	}
// end outputStruct
	
struct outputStructPath {
	real colvector lambdalist
	real matrix betas
	real matrix sbetas
	real rowvector Ups
	real rowvector sdvec
	real rowvector sUps
	real colvector dof // "effective" degrees of freedom
	real colvector shat // number of non-zero parameters
	real scalar cons // yes/no
	real rowvector intercept
	real rowvector interceptPL
	real scalar n
	real scalar nclust
	}
// end outputStructPath


struct dataStruct scalar MakeData(	string scalar nameY,
									string scalar nameX,
									string scalar nameX_o,
									string scalar touse,
									real scalar cons,
									real scalar dmflag,			///
									|							///  optional arguments
									string scalar stdymat,		///
									string scalar stdxmat,		///
									string scalar nameclustid,	///  optional arguments - rlasso-specific
									string scalar nameP,		///
									string scalar nameNP)
{

	if (args()<=6) {
		stdymat		= ""
		stdxmat		= ""
		nameclustid	= ""
	}
	if (args()<=9) {
		nameP = ""				//  default list is empty
		nameNP = ""				//  default list is empty
	}

	struct dataStruct scalar d
	

	// dep var
	st_view(y,.,nameY,touse)
	d.y	=&y

	// X vars
	st_view(X,.,nameX,touse)
	d.X			=&X
	d.nameX		=tokens(nameX)
	d.nameX_o	=tokens(nameX_o)
	d.n			=rows(X)
	d.p			=cols(X)
	d.cons		=cons
	d.dmflag	=dmflag

	// cluster var
	if (nameclustid!="") {
		st_view(cid,.,nameclustid,touse)
		d.clustid	= &cid
		info		= panelsetup(cid, 1)
		d.nclust	= rows(info)
	} 
	else {
		d.nclust	= 0
	}

	// mean vectors
	if (dmflag) {
		// already demeaned
		d.mvec		= J(1,cols(*d.X),0)
		d.ymvec		= 0
	}
	else {
		// calculate means
		d.mvec		= mean(*d.X)
		d.ymvec		= mean(*d.y)
	}

	// standardization vectors
	if ((stdymat=="") & (dmflag)) {
		// already demeaned (mean zero)
		d.ysd		= sqrt(mean((*d.y):^2))
		d.varvec	= mean((*d.X):^2)
		d.sdvec	= sqrt(d.varvec)
		d.prestdflag= 0
	}
	else if (stdymat=="") {
		// not mean zero so need to demean
		d.ysd		= sqrt(mean(((*d.y):-d.ymvec):^2))
		d.varvec	= mean(((*d.X):-d.mvec):^2)
		d.sdvec	= sqrt(d.varvec)
		d.prestdflag= 0
	}
	else {
		// data are prestandardized so just store these
		d.ysd		=st_matrix(stdymat)
		d.sdvec	=st_matrix(stdxmat)
		d.varvec	=(d.sdvec):^2
		d.prestdflag= 1
	}

	// X'X and X'y
	if (cons) {
		d.XX = quadcrossdev(*d.X,d.mvec,*d.X,d.mvec)
		d.Xy = quadcrossdev(*d.X,d.mvec,*d.y,d.ymvec)
	}
	else {
		d.XX = quadcross(*d.X,*d.X)
		d.Xy = quadcross(*d.X,*d.y)
	}	

	if (nameNP!="") {
		// unpenalized regressors in rlasso
		st_view(Xp,.,nameP,touse)
		st_view(Xnp,.,nameNP,touse)
		d.Xp		=&Xp
		d.Xnp		=&Xnp
		d.nameXp	=tokens(nameP)
		d.nameXnp	=tokens(nameNP)
		selXp		= J(1,cols(*d.X),1)
		forbound	= cols(d.nameXnp)		//  faster
		for (i=1;i<=forbound;i++) {
			selXp = selXp - (d.nameX :== d.nameXnp[1,i])
		}
		d.selXp		=selXp
		d.selXnp	=1:-selXp
		d.selindXp	=selectindex(d.selXp)
		d.selindXnp	=selectindex(d.selXnp)

		if (cons) {
			// model has a constant (unpenalized)
			d.mvecp		= mean(*d.Xp)
			d.mvecnp	= mean(*d.Xnp)
			d.ypihat	= qrsolve((*d.Xnp):-d.mvecnp,(*d.y):-d.ymvec)
			d.ysdp		= sqrt(mean((((*d.y):-d.ymvec)-((*d.Xnp):-mean(*d.Xnp))*d.ypihat):^2))
			d.pihat		= qrsolve((*d.Xnp):-d.mvecnp,(*d.Xp):-d.mvecp)
			// std vector for just the Xp variables.
			d.sdvecp	= sqrt(mean((((*d.Xp):-d.mvecp)-((*d.Xnp):-mean(*d.Xnp))*d.pihat):^2))
		}
		else {
			// model has no constant
			d.mvecp		= J(1,cols(*d.Xp),0)
			d.mvecnp	= J(1,cols(*d.Xnp),0)
			d.ypihat	= qrsolve(*d.Xnp,*d.y)
			d.ysdp		= sqrt(mean(((*d.y)-(*d.Xnp)*d.ypihat):^2))
			d.pihat		= qrsolve(*d.Xnp,*d.Xp)
			// std vector for just the Xp variables.
			d.sdvecp	= sqrt(mean(((*d.Xp)-(*d.Xnp)*d.pihat):^2))
		}

		// Now create blank full-size std vector of zeros for all X variables.
		d.sdvecpnp					= J(1,cols(*d.X),0)
		// Then insert values in appropriate columns.
		d.sdvecpnp[1,(d.selindXp)]	= d.sdvecp
	}
	else {													//  no unpenalized regressors (in rlasso)
		d.Xp		= d.X									//  penalized = all
		d.nameXp	= d.nameX
		d.selXp		= J(1,cols(*d.X),1)
		d.selindXp	= selectindex(d.selXp)
		d.mvecp		= d.mvec
		d.pihat		= 0
		d.ypihat	= 0
		if (d.prestdflag) {
			d.ysdp		= sqrt(mean(((*d.y):-d.ymvec):^2))
			d.sdvecp	= sqrt(mean(((*d.X):-d.mvec):^2))
			d.sdvecpnp	= d.sdvecp
		}
		else {
			d.ysdp		= d.ysd
			d.sdvecp	= d.sdvec
			d.sdvecpnp	= d.sdvec
		}
	}

	return(d)
}
// end MakeData



// this function calculates lasso path for range of lambda values
struct outputStructPath DoLassoPath(struct dataStruct scalar d,
									real rowvector Ups,					//  vector of penalty loadings (L1 norm)
									real rowvector Ups2,				//  vector of penalty loadings (L2 norm)
									real rowvector lvec,				//  vector of lambdas (L1 norm)
									real rowvector lvec2,				//  vector of lambdas (L2 norm)
									real scalar post,
									real scalar verbose,
									real scalar optTol,
									real scalar maxIter,
									real scalar zeroTol,
									real scalar alpha,
									real scalar lglmnet,
									real scalar noic)
{

		struct outputStructPath scalar t

		p = cols(*d.X)
		n = rows(*d.y)
		
		XX = d.XX
		Xy = d.Xy
		
		if (verbose>=2) {
			printf("Lambda list: %s\n",invtokens(strofreal(lvec)))
		}
		
		lmax=max(lvec)
		lcount=cols(lvec)
		beta=lusolve(XX+lmax/2*diag(Ups2),Xy) // beta start	
		if (verbose==3) {
			printf("Initial beta:\n")
			beta'
		}
		
		XX2=XX*2
		Xy2=Xy*2
		
		lpath = J(lcount,p,.) // create empty matrix which stores coef path
		for (k = 1;k<=lcount;k++) { // loop over lambda
			
			// Separate blocks for lasso, ridge, elastic net.
			// Separate lambdas and penalty loadings for L1 and L2 norms to accommodate standardization.
			// If data were pre-standardized, then lambda=lambda2 and Ups=Ups2.
			// If standardization is on-the-fly and incorporated in penalty loadings,
			// then lambdas and penalty loadings for L1 and L2 norms are different.
			lambda=lvec[1,k]
			lambda2=lvec2[1,k]
			m=0
			change = optTol*2 // starting value
			// optimization for a given lambda value. 
			while ((m < maxIter) & (change>optTol)) { 
			
				beta_old = beta
				for (j = 1;j<=p;j++) // keep all beta fixed, except beta[j]. update estimate.
				{
					S0 = quadcolsum(XX2[j,.]*beta) - XX2[j,j]*beta[j] - Xy2[j]

					if (alpha==1) {					//  lasso
							if (S0 > lambda*Ups[j])
							{
								beta[j] = (lambda*Ups[j] - S0)/(XX2[j,j])
							}
							else if (S0 < -lambda*Ups[j])	
							{
								beta[j] = (-lambda*Ups[j] - S0)/(XX2[j,j]) 
							}
							else 
							{
								beta[j] = 0
							}
					}								//  end lasso
					else if (alpha>0) {				//  elastic net
							if (S0 > lambda*Ups[j]*alpha)
							{
								beta[j] = (lambda*Ups[j]*alpha - S0)/(XX2[j,j]+ lambda2*Ups2[j]*(1-alpha))
							}
							else if (S0 < -lambda*Ups[j]*alpha)	
							{
								beta[j] = (-lambda*Ups[j]*alpha - S0)/(XX2[j,j]+ lambda2*Ups2[j]*(1-alpha)) 
							}
							else 
							{
								beta[j] = 0
							}
					} 								//  end elastic net
					else if (alpha==0) {			//  ridge  
					
							if (S0 > 0)
							{
								beta[j] = (-S0)/(XX2[j,j] + lambda2*Ups2[j])
							}
							else if (S0 < 0)	
							{
								beta[j] = (-S0)/(XX2[j,j] + lambda2*Ups2[j]) 
							}
							else 
							{
								beta[j] = 0
							}
					}								//  end ridge
				}									//  end j loop over components of beta

				m++
				change = quadcolsum(abs(beta-beta_old)) 
			
			}
			lpath[k,.]=beta'
		}	
		
		if (verbose==3) {
			printf("beta after shooting:\n")
			lpath[1,.]
		}	
		
		// following code should be the same for DoLassoPath() and DoSqrtLassoPath()
		
		lpath=edittozerotol(lpath, zeroTol)
		
		if (post) { 
			betasPL = J(lcount,p,0)
			nonzero0 = J(1,p,0)
			for (k = 1;k<=lcount;k++) { // loop over lambda points
				nonzero = lpath[k,.]:!=0  // 0-1 vector
				sk = sum(nonzero)			
				if ((0<sk) & (sk<n)) { // only if 0<s<n
					if ((nonzero0==nonzero) & (k>=2)) { // no change in active set
						betasPL[k,.] = betasPL[k-1,.]
					}
					else { 
						ix = selectindex(nonzero)	// index of non-zeros
						// obtain beta-hat
						if (d.dmflag) {
							// data are mean zero
							betak=qrsolve(select((*d.X),nonzero),*d.y)
						}
						else {
							betak=qrsolve(select((*d.X),nonzero):-select(d.mvec,nonzero),((*d.y):-d.ymvec))
						}
						betasPL[k,ix] = betak'
						nonzero0=nonzero
					}
				}
			}
			t.betas 		= betasPL	
		}
		else {
			t.betas 		= lpath
		}
		
		// use glmnet lambda
		if (lglmnet) {
			lvec=lvec/2/n
		}
		
		t.lambdalist	= lvec'
		t.Ups			= Ups
		if (d.dmflag) {
			// data are mean zero => intercept is zero
			t.intercept		= 0
		}
		else {
			t.intercept		= mean(*d.y):-mean((*d.X))*(t.betas')
		}
		t.cons 			= d.cons
		
		// degrees of freedom and dimension of the model (add constant)
		t.shat = quadrowsum(t.betas:!=0) :+ (d.cons | d.dmflag)
		if (!noic) {
			if (alpha==1) { 
				// lasso dof
				// need to add the constant / account for demeaning
				t.dof	= t.shat 
			}
			else if (alpha==0) {
				// ridge dof  
				df = J(lcount,1,.)
				if (d.cons) {
					Xt = (*d.X) :- d.mvec
				}
				else {
					Xt=(*d.X)
				}
				for (k=1;k<=lcount;k++) { // loop over lambda points
					df[k,1] =trace((Xt)*invsym(quadcross(Xt,Xt):+lvec2[1,k]/2*diag(Ups2))*(Xt)') 
				}
				// need to add the constant / account for demeaning
				df = df :+ (d.cons | d.dmflag)
				t.dof	= df  
			}
			else {
				// elastic net dof
				df = J(lcount,1,.)
				if (d.cons) {
					Xt = (*d.X) :- d.mvec
				}
				else {
					Xt=(*d.X)
				}				
				for (k=1;k<=lcount;k++) { // loop over lambda points
					nonzero = lpath[k,.]:!=0  // 0-1 vector
					XA=select((Xt),nonzero)
					Ups2A=select((Ups2),nonzero)
					df[k,1] =trace((XA)*invsym(quadcross(XA,XA):+(1-alpha)*lvec2[1,k]/2*diag(Ups2A))*(XA)')
				}
				// need to add the constant / account for demeaning
				df = df :+ (d.cons | d.dmflag)
				t.dof	= df  
			}
		}
				
		return(t)		

}
// end DoLassoPath

	
void ReturnResultsPath(		struct outputStructPath scalar t,	///  
							struct dataStruct scalar d,			///  
							string scalar Xnames,				///
							|									///
							real scalar sqrtflag,				///
							real scalar stdcoef)
{
		// default values
		if (args()<=3) sqrtflag = 0
		if (args()<=4) stdcoef = 0

		Xnamesall=tokens(Xnames)
		betas			= t.betas
		lambdalist		= t.lambdalist
		wl1norm			= rowsum(abs(betas :* t.Ups))

		// unstandardize unless overriden by stdcoef
		if (d.prestdflag & stdcoef==0) {
			betas		= betas			:/ d.sdvec * d.ysd
			if (sqrtflag==0) {
				// sqrt-lasso lambdas don't need unstandardizing (pivotal so doesn't depend on sigma/y)
				lambdalist	= lambdalist	* d.ysd
			}
			wl1norm		= wl1norm		* d.ysd
		}

		// "L1 norm" excluding constant and unpenalized vars
		l1norm			= rowsum(abs(betas :* (t.Ups :> 0)))

		//dof				= quadrowsum(betas:!=0)
		pall			= cols(Xnamesall)

		if (t.cons) {
			betas		= (betas , (t.intercept'))		//  no intercept if pre-standardized
			Xnamesall	= (Xnamesall, "_cons")
			pall		= pall+1
		}
		st_numscalar("r(lmax)", max(lambdalist))
		st_numscalar("r(lmin)", min(lambdalist))
		st_numscalar("r(lcount)",rows(lambdalist))
		st_matrix("r(lambdalist)",lambdalist)
		st_matrix("r(l1norm)",l1norm)
		st_matrix("r(wl1norm)",wl1norm)
		st_matrix("r(betas)",betas)
		st_matrix("r(Ups)",t.Ups)
		st_matrix("r(sUps)",t.sUps)
		st_matrix("r(shat)",t.shat)
		st_matrix("r(stdvec)",d.sdvec)
		st_matrixcolstripe("r(betas)",(J(pall,1,""),Xnamesall'))
		st_matrixcolstripe("r(lambdalist)",("","Lambdas"))
		st_matrixcolstripe("r(l1norm)",("","L1norm"))
		st_matrixcolstripe("r(wl1norm)",("","wL1norm"))
}
// end ReturnResultsPath
	
void ReturnCVResults(real rowvector Lam, ///
					real rowvector MSPE,
					real rowvector SD, 
					real scalar minid,
					real scalar minseid)
{
	
	lnum=cols(Lam)
	
	printf("{txt}%10s{c |} {space 3} {txt}%10s {space 3} {txt}%10s {space 3} {txt}%10s\n","","Lambda","MSPE","st. dev.")
	printf("{hline 10}{c +}{hline 45}\n")
			
	for (j = 1;j<=lnum;j++) 	{
		
		if (j==minid) {
			marker="*"
		}
		else {
			marker=""
		}
		if (j==minseid) {
			marker=marker+"^"
		}
		printf("{txt}%10.0g{c |} {space 3} {res}%10.0g {space 3} {res}%10.0g {space 3} {res}%10.0g  %s\n",j,Lam[1,j],MSPE[1,j],SD[1,j],marker)
	
	}
}
// end ReturnCVResults

real scalar lambdaCalc(struct dataStruct scalar d,		///
						real scalar pminus,				/// adjustment to number of Xs in model
						real scalar gamma,				///
						real scalar gammad,				///
						real scalar c,					/// 
						real scalar R,					///
						real scalar hetero,				///
						real scalar xdep,				///
						real colvector ehat,			///
						real scalar rmse,				///
						real scalar lalt,				///
						|								///
						real scalar newseed,			///
						real scalar dotsflag			///
						) 
{

	if (args()==11) {								//  default is not to set the seed
		newseed = -1
	}

	// model may have partialled-out var with coeffs estimated separately; subtract in formulae below
	p = d.p - pminus

	if (xdep==0) {									// X-independent
		if (lalt==0) {								// standard lambda
			lambda=2*c*sqrt(d.n)*invnormal(1-gamma/gammad/(2*p))
		}
		else {										//  alternative lambda
			lambda = 2 * c * sqrt(d.n) * sqrt(2 * log(2*p/(gamma/gammad)))
		}
	}												//  end X-independent block
	else {											//  X-dependent
		if (newseed>-1) {
			rseed(newseed)							//  set seed if requested
		}
		sim=J(R,1,.)								//  placeholder for simulated values of the distribution
		if (dotsflag) {
			dotscmd = "_dots 0 0, title(Estimating x-dependent lambda using " + strofreal(R) + " repetitions)"
			stata(dotscmd)
		}
		for (j = 1;j<=R;j++) 	{					//  simulate
			if (dotsflag) {
				dotscmd = "_dots " + strofreal(j) + " 0"
				stata(dotscmd)
			}
			if ((hetero==0) & (d.nclust==0)) {		//  X-dependent & homoskedastic
 				g=rnormal(d.n,1,0,1)
				// treatment of Xs depends on whether unpenalized Xs exist and/or Xs have been pre-standardized
				// Xs should be demeaned and standardized
				// centerpartial(.) returns Xp demeaned and with Xnp partialled-out.
				sim[j,1]=max(quadcolsum((centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat):/(d.sdvecp)):*g))
			}
			else {	 								//  X-dependent & (heteroskedastic or clustered)
				if (d.nclust==0) {					//  X-dependent & heteroskedastic
					g=rnormal(d.n,1,0,1)
				}
				else {								//  X-dependent & clustered
					// clustered case - this is provisional!!
					// g is iid by cluster and repeated within clusters
					g=J(d.n,1,.)
					info = panelsetup(*d.clustid, 1)
					forbound = rows(info)			//  faster
					for (i=1; i<=forbound; i++) {
						gi	=rnormal(1,1,0,1)
						// next line expands scalar gi into a vector of identical values
						gi	=gi*J(info[i,2]-info[i,1]+1,1,1)
						// now insert into g
						g[info[i,1]..info[i,2],1] = gi
					}
				}
				gehat=g:*ehat
				// treatment of Xs depends on whether unpenalized Xs exist and/or Xs have been pre-standardized
				// Xs should be demeaned and standardized
				// divide by rmse because we are multiplying x by g*ehat (rather than just ehat as in homosked case)
				// centerpartial(.) returns Xp demeaned and with Xnp partialled-out.
				sim[j,1]=max(quadcolsum((centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat):/(d.sdvecp)):*gehat)) * 1/rmse
			}
		}	// end simulation loop

		// use quantile from simulated values to get lambda
		lambda=2*c*mm_quantile(sim,1,1-(gamma/gammad))
	}

	return(lambda)
}
// end lambdaCalc

real colvector InitialResiduals(	struct dataStruct scalar d,
									real scalar corrnumber) 
{ 
	// applies OLS to a reduced set of regressors exhibiting highest correlation with y
	// size of restricted set = corrnumber

	// in case corrnumber < dim(X)
	if (d.pihat==0) {
		dimX=cols(*d.X)
	}
	else {
		dimX=cols(*d.Xp)
	}
	corrnumber=min((corrnumber,dimX))

	// just return if corrnum = 0
	if (corrnumber <= 0) {
		// centerpartial(.) returns y centered and with Xnp partialled out.
		return(centerpartial(d.y,d.ymvec,d.cons,d.Xnp,d.mvecnp,d.ypihat))
	}
	else {
		dimZ=dimX+1
		// corr matrix
		if ((d.pihat==0) & (d.dmflag)) {			//  no notpen Xs, zero-mean data
			Z=abs(correlation( (*d.y,*d.X)) )										//  join y and X
		}
		else if (d.pihat==0) {						//  no notpen Xs
			Z=abs(correlation(														///
								((*d.y):-d.ymvec,									/// join demeaned y and X
								 (*d.X):-d.mvec)									///
								))
		}
		else if (d.dmflag) {						//  notpen Xs, zero-mean data
			Z=abs(correlation(														///
								((*d.y)  :- (*d.Xnp)*d.ypihat,						/// join demeaned and
								 (*d.Xp) :- (*d.Xnp)*d.pihat)						/// partialled y and X
								))
		}
		else {										//  notpen Xs
			Z=abs(correlation(														///
								((*d.y):-d.ymvec  :- ((*d.Xnp):-d.mvecnp)*d.ypihat,	/// join demeaned and
								 (*d.Xp):-d.mvecp :- ((*d.Xnp):-d.mvecnp)*d.pihat)	/// partialled y and X
								))
		}
		z=Z[2..dimZ,1]
		ix=order(-z,1)
		ix=ix[1..corrnumber,1]

		if ((d.pihat==0) & (d.dmflag)) {				//  no notpen Xs, zero-mean data
			b	=qrsolve(								///
						(*d.X)[.,ix],					///
						 *d.y							///
						 )
			r	=(*d.y) :- ((*d.X)[.,ix])*b
		}
		else if (d.pihat==0) {							//  no notpen Xs
			b	=qrsolve(								///
						((*d.X)[.,ix]):-d.mvec[1,ix],	///
						(*d.y):-d.ymvec					///
						)
			r	=(*d.y) :- d.ymvec :- (((*d.X)[.,ix]):-d.mvec[1,ix])*b
		}
		else if (d.dmflag) {							//  notpen Xs, zero-mean data
			// ix has highest-correlation cols of Xp
			// replace with corresponding cols of X
			ix	= d.selindXp[1,ix]
			// and append np cols of X
			ix	= ix,d.selindXnp
			b	=qrsolve(								///
						 (*d.X)[.,ix],					/// selected cols of Xp plus all of Xnp
						 *d.y							///
						 )
			r	=(*d.y) :- ((*d.X)[.,ix])*b
		}
		else {											//  notpen Xs
			// ix has highest-correlation cols of Xp
			// replace with corresponding cols of X
			ix	= d.selindXp[1,ix]
			// and append np cols of X
			ix	= ix,d.selindXnp
			b	=qrsolve(								///
						 ((*d.X)[.,ix]):-d.mvec[1,ix],	/// selected cols of Xp plus all of Xnp
						 *d.y :- d.ymvec				///
						 )
			r	=(*d.y) :- d.ymvec :- (((*d.X)[.,ix]):-d.mvec[1,ix])*b
		}
		return(r)				// return residuals

	}
}
// end InitialResiduals

real rowvector MakeLassoWeights(struct dataStruct scalar d,	///
								real colvector v,			///
								real scalar rmse,			///
								real scalar hetero,			///
								|							///
								real scalar nclust1,		///
								real scalar center			///
								)
{

	if (args()==4) {
		nclust1 = 0											//  default value of 0 serves as boolean below
	}
	if (args()<=5) {
		center = 0
	}

	real rowvector Ups

	if (d.nclust==0) {										//  no cluster dependence
		if (hetero==0) {									//  homoskedastic
			Ups = d.sdvecpnp								//  this is the partialled-out version with zero penalties inserted.
															//  nb: rmse is in lambda
		}
		else if ((hetero==1) & (center)) {		
			// in lasso, rmse incuded in het Ups for standardization purposes only; also included in lambda
			// compare sqrt-lasso, where rmse is required in non-homosk Ups and does not appear in lambda (pivotal)
			// centerpartial(.) returns Xp demeaned and with Xnp partialled-out.
			centervec =	mean( centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat) :* v)
			St = J(1,cols(*d.X),0)
			St[1,(d.selindXp)] = quadcolsum( ( (centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat) :* v) :- centervec ):^2 )
			Ups = sqrt(St/d.n)*1/rmse						// relevant n is #obs
		}
		else {
			St = J(1,cols(*d.X),0)
			St[1,(d.selindXp)] = quadcolsum( ( (centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat) :* v)              ):^2 )
			Ups = sqrt(St/d.n)*1/rmse						// relevant n is #obs
		}
	}
	else {													//  cluster dependence
		//  this code is memory-intensive because St is a matrix that is the same dimension as X
		//  should recode to at least loop over clusters
		//  relevant n is #clusters
		info = panelsetup(*d.clustid, 1)
		//  St is nT x p; each element is x_it*e_it.
		// First create blank full-size matrix of zeros to be populated.
		St = J(rows(*d.X),cols(*d.X),0)
		// Now populate correct columns, leaving other as zeros (unpenalized).
		// centerpartial(.) returns Xp demeaned and with Xnp partialled-out.
		St[.,(d.selindXp)] = centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat):*(v*J(1,cols(*d.Xp),1))
		//  UpsTemp is n x p. Each row corresponds to a cluster.
		UpsTemp = J(d.nclust,d.p,0)
		forbound = d.nclust							//  faster
		for (i=1; i<=forbound; i++) {	
			Sti = panelsubmatrix(St, i, info)
			UpsTemp[i,.] = quadcolsum(Sti)			//  Populate each row with the cluster-sum of x_it*e_it for cluster i.
		}
		//  center UpsTemp
		if (center) {
			UpsTemp = UpsTemp - J(d.nclust,1,1)*quadcolsum(UpsTemp)/d.nclust
		}
		//  Default approach is simply to divide by nobs = nclust*T in a balanced panel.
		//  A finite-sample adjustment as in CBH's lassoCluster is to divide by (nclust-1)*T.
		//  In an unbalanced panel, we achieve this with 1/nobs * nclust/(nclust-1).
		//  In a balanced panel, 1/nobs * nclust/(nclust-1) = 1/(nclust*T) * nclust/(nclust-1) = 1/((nclust-1)*T).
		if (nclust1) {										//  override default: divide by (nclust-1)*T
			Ups = sqrt(quadcolsum(UpsTemp:^2)/(d.n)*(d.nclust)/(d.nclust-1))
		}
		else {												//  default: simply divide by n = nobs*T
			Ups = sqrt(quadcolsum(UpsTemp:^2)/(d.n))
		}
		//  in lasso, rmse incuded in homosked Ups for standardization purposes only; also included in lambda
		//  compare sqrt-lasso, where rmse is required in non-homosk Ups and does not appear in lambda (pivotal)
		Ups = Ups * 1/rmse
	}
	return(Ups)
}
// end MakeLassoWeights

void EstimateRLasso(							///  Complete Mata code for RLasso.
				string scalar nameY,			///
				string scalar nameX,			///
				string scalar nameX_o,			///
				string scalar pen,				///
				string scalar notpen,			///
				string scalar touse,			///
				string scalar nameclustid, 		///
				string scalar stdymat,			///
				string scalar stdxmat,			///
				real scalar sqrtflag,			/// lasso or sqrt-lasso?
				real scalar hetero,				/// homosk or heteroskedasticity?
				real scalar xdep,				/// X-dependent or independent?
				real scalar R,					/// number of simulations with xdep
				real scalar lassoUps,			/// use lasso or post-lasso residuals for estimating penalty loadings?
				real scalar optTol,				///
				real scalar maxIter,			///
				real scalar zeroTol,			///
				real scalar maxUpsIter,			///
				real scalar UpsTol,				///
				real scalar verbose,			///
				real scalar c,					///
				real scalar gamma,				///
				real scalar gammad,				///
				real scalar lambda0,			///
				real scalar lalt, 				///
				real scalar corrnumber,			///
				real scalar pminus,				///
				real scalar nclust1,			/// use #nclust-1 instead of #nclust in cluster-lasso
				real scalar center,				/// center x_i*e_i or x_ij*e_ij in cluster-lasso
				real scalar supscoreflag,		///
				real scalar ssnumsim,			///
				real scalar newseed,			/// rnd # seed; relevant for xdep and supscore
				real scalar dotsflag,			///
				real scalar cons,				///
				real scalar dmflag)
{
	struct dataStruct scalar d
	d = MakeData(nameY,nameX,nameX_o,touse,cons,dmflag,stdymat,stdxmat,nameclustid,pen,notpen)

	struct outputStruct scalar OUT
	if (sqrtflag) {
		OUT	= RSqrtLasso(d,hetero,xdep,R,lassoUps,optTol,maxIter,zeroTol,maxUpsIter,UpsTol,verbose,c,gamma,gammad,lambda0,lalt,corrnumber,pminus,nclust1,center,supscoreflag,ssnumsim,newseed,dotsflag)
	}
	else {
		OUT	= RLasso(d,hetero,xdep,R,lassoUps,optTol,maxIter,zeroTol,maxUpsIter,UpsTol,verbose,c,gamma,gammad,lambda0,lalt,corrnumber,pminus,nclust1,center,supscoreflag,ssnumsim,newseed,dotsflag)
	}
	
	ReturnResults(OUT,d,sqrtflag)		//  puts results into r(.) macros
}	// end EstimateRLasso

struct outputStruct scalar RLasso(							/// Mata code for BCH rlasso
							struct dataStruct scalar d,		/// data
							real scalar hetero,				/// homosk or heteroskedasticity?
							real scalar xdep,				/// X-dependent or independent?
							real scalar R,					/// number of simulations with xdep
							real scalar lassoUps,			/// use lasso or post-lasso residuals for estimating penalty loadings?
							real scalar optTol,				///
							real scalar maxIter,			///
							real scalar zeroTol,			///
							real scalar maxUpsIter,			///
							real scalar UpsTol,				///
							real scalar verbose,			///
							real scalar c,					///
							real scalar gamma,				///
							real scalar gammad,				///
							real scalar lambda0,			///
							real scalar lalt, 				///
							real scalar corrnumber,			///
							real scalar pminus,				///
							real scalar nclust1,			/// use #nclust-1 instead of #nclust in cluster-lasso
							real scalar center,				/// center x_i*e_i or x_ij*e_ij in cluster-lasso
							real scalar supscoreflag,		///
							real scalar ssnumsim,			///
							real scalar newseed,			///
							real scalar dotsflag			///
							)
{
	struct outputStruct scalar betas

	alpha=1 // elastic net parameter. always one.
	if (gammad<=0) {										//  not user-provided, so set here
		if (d.nclust==0)	gammad=log(d.n)					//  not cluster-lasso so #obs=n
		else				gammad=log(d.nclust)			//  cluster-lasso so #obs=nclust
	}
	
	// initial residuals, lambda, loadings
	v = InitialResiduals(d,corrnumber)
	s1 = sqrt(mean(v:^2))

	// calculate lambda0 unless provided by user
	if (lambda0==0) {
		lambda0=lambdaCalc(d,pminus,gamma,gammad,c,R,hetero,xdep,v,s1,lalt,newseed,dotsflag)
	}

	//  rmse is in lambda, not in lambda0
	lambda=lambda0*s1

	Ups1 = MakeLassoWeights(d,v,s1,hetero,nclust1,center)
	if (lassoUps==0) {
		Ups1=Ups1:*0.5	// initial penalty is multiplied by 0.5 as in rlasso/CBH's program [post-lasso only]
	}
	Delta=UpsTol*2														// initial value for "Delta"

	// updating
	MinIter=2		// need this to force minimum of 2 iterations and updating of Ups
	iter = 0
	do {																// start loop; always iterated at least once
		// increment counter
		iter = iter + 1
		Ups0=Ups1
		s0=s1
		if (verbose>=1) {
			printf("Estimation of penalty level/loadings: Step %f.\n",iter)
			printf("RMSE: %f\n",s0)
		}
		// lasso estimation for given penalty level & loadings; betas based on Ups0.

		// was: betas = DoLasso(d, Ups0, ., lambda, ., verbose, optTol, maxIter, zeroTol, alpha)
		betas = DoLasso(d, Ups0, Ups0, lambda, lambda, verbose, optTol, maxIter, zeroTol, alpha)
		// obtain residuals; based on betas(Ups0)
		if (lassoUps==1) {
			v = betas.v  //y - select(X,betas.index')*betas.beta					// lasso residuals
			s1 = betas.rmse 
		}
		else {
			v = betas.vPL 	//v = y - select(X,betas.index')*betas.betaPL	 				// post-lasso residuals
			s1 = betas.rmsePL
		}
		// s1 = rmse based on new residuals v(betas(Ups0)); use use large-sample definition of variance

		// change in RMSE				
		Delta = abs(s1-s0)

		// Reporting
		// unclear why cols need to be >1. Need c>1 even if r==0; c==0 and r==0 makes invtokens(.) crash.
		if (verbose>=1 & cols(betas.nameXSel)>1) {
			printf("Selected variables: %s\n",invtokens(betas.nameXSel'))
		}
		if (verbose>=1) {
			printf("Change in RMSE: %f\n\n",Delta)
		}

		// update Ups1 based on new residuals v, but only if loop is going to repeat
		// also update lambda (heterosk + X-dependent only), again only if loop is going to repeat
		if ( ((iter < maxUpsIter) & (Delta > UpsTol)) | (iter < MinIter) ) {

			// xdep lambda update
			if ((hetero==1) & (xdep==1)) {
				// see above; convention is that rmse is in lambda, but not in lambda0
				lambda0=lambdaCalc(d,pminus,gamma,gammad,c,R,hetero,xdep,v,s1,lalt,newseed,dotsflag)
				lambda=lambda0*s1
			}
			// Ups update
			Ups1 = MakeLassoWeights(d,v,s1,hetero,nclust1,center)
			// rmse is in lambda, not in lambda0
			lambda=lambda0*s1
		}

	} while ( ((iter < maxUpsIter) & (Delta > UpsTol)) | (iter < MinIter) )
	// When loop completes, we have: betas(Ups0); s1 based on v(betas(Ups0)); betas, rmse, Ups all consistent.
	// Note that when loop completes, we have avoided updating Ups1 based on v(betas(Ups0)).

	if (verbose>=1) {
		printf("Number of penalty loading iterations: %g\n",iter)
		if (iter == maxUpsIter) {
			printf("Warning: reached max penalty loading iterations w/o achieving convergence.\n")
		}
		else {
			printf("Penalty loadings (upsilon) convergence achieved.\n")
		}
	}

	// sup-score stat
	if (supscoreflag) {
		betas.supscore	= doSupScore(d, c, gamma, pminus, hetero, verbose, ssnumsim, newseed, dotsflag)
	}

	betas.sUps		= betas.Ups :/ d.sdvecpnp				//  should be =1 under homosk.
	betas.sUps		= editmissing(betas.sUps,0)				//  in case any unpenalized (div by 0)

	// Misc
	betas.n			= d.n
	betas.nclust	= d.nclust
	betas.nupsiter	= iter
	betas.lambda0	= lambda0
	betas.slambda	= lambda / d.ysdp

	betas.c			= c
	betas.gamma		= gamma
	betas.gammad	= gammad

	return(betas)
}
// end RLasso	

void EstimateSupScore(							///  Complete Mata code for RLasso.
				string scalar nameY,			///
				string scalar nameX,			///
				string scalar nameX_o,			///
				string scalar pen,				///
				string scalar notpen,			///
				string scalar touse,			///
				string scalar nameclustid, 		///
				string scalar stdymat,			///
				string scalar stdxmat,			///
				real scalar hetero,				/// homosk or heteroskedasticity?
				real scalar verbose,			///
				real scalar R,					///
				real scalar c,					///
				real scalar gamma,				///
				real scalar pminus,				///
				real scalar newseed,			///
				real scalar dotsflag,			///
				real scalar cons,				///
				real scalar dmflag)
{
	struct dataStruct scalar d
	d = MakeData(nameY,nameX,nameX_o,touse,cons,dmflag,stdymat,stdxmat,nameclustid,pen,notpen)

	struct outputStruct scalar OUT
	OUT.supscore	= doSupScore(d, c, gamma, pminus, hetero, verbose, R, newseed, dotsflag)
	// Misc
	OUT.n			= d.n
	OUT.nclust		= d.nclust
	OUT.c			= c
	OUT.gamma		= gamma

	ReturnResults(OUT,d)		//  puts results into r(.) macros

}	// end EstimateRLasso



real rowvector doSupScore(								///
							struct dataStruct scalar d,	///
							real scalar c,				///
							real scalar gamma,			///
							real scalar pminus,			///
							real scalar hetero,			///
							|							///
							real scalar verbose,		///
							real scalar R,				///
							real scalar newseed,		///
							real scalar dotsflag		///
							)
{

	// default is verbose=0
	if (args()<6) {
		verbose = 0
	}
	// default is reps=500
	if (args()<7) {
		R = 500
	}
	// default is no new seed
	if (args()<8) {
		newseed = -1
	}
	// default is no dots
	if (args()<9) {
		dotsflag = 0
	}

	// Null hypoth implies E(ynull*x)=0 where ynull and x are demeaned and standardized.
	if (d.dmflag) {
		// data are already demeaned
		ynull = (*d.y) * 1/d.ysd
	}
	else {
		ynull = ((*d.y):-(d.ymvec)) * 1/d.ysd
	}

// ******* sup-score test stat ********** //

// Follows Matlab code from BCH.
// Loop over jj values of aVec: assemble ynull (="eTemp") then calc SupScore for that jj ynull
// % Sup-Score Test
// aVec = (-.5:.001:1.5)';
// SupScore = zeros(size(aVec));
// for jj = 1:size(aVec,1)
//     aT = aVec(jj,1);
//     eTemp = My-Md*aT;
//     ScoreVec = eTemp'*Mz;
//     ScoreStd = sqrt((eTemp.^2)'*(Mz.^2));
//     ScaledScore = ScoreVec./(1.1*ScoreStd);
//     SupScore(jj,1) = max(abs(ScaledScore));
// end

// NB: can't use quadcross or mean with cde below because a column can have (all) missing (e.g. if std dev=0).
//     and quadcross and mean use rowwise deletion (all missing => all rows dropped!)
//     use 1/n * quadcolsum instead

	// centerpartial(.) returns Xp demeaned and with Xnp partialled-out.
	ScoreVec	= 1/(d.n) * quadcolsum( (centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat):/(d.sdvecp)) :* ynull )
	if ((d.nclust) | (hetero)) {
		// don't need this for the homoskedastic case
		ScoreStd	= sqrt( 1/(d.n) * quadcolsum( (((centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat):/(d.sdvecp)):^2) :* (ynull:^2) ) ) )
	}
// Question: if homoskedastic, isn't E(ScoreStd)=E(ynull^2)E(x^2)=1?
	if (d.nclust) {
		supscore	= (d.nclust)*max(abs(ScoreVec:/ScoreStd))
	}
	else if (hetero) {
		supscore	= (d.n)*max(abs(ScoreVec:/ScoreStd))
	}
	else {
		supscore	= (d.n)*max(abs(ScoreVec))
	}


// Question: BCH paper says (p. 20 of PDF) that gi is iid indep of controls and IVs but should use residuals
// after projecting gi on controls w  ... but if gi are iid independent of controls, why bother?

// Following rlasso R code:
//     object$supscore <- sqrt(n)*max(abs(colMeans(object$model*as.vector(object$dev))))

	// use 1/(d.n) * quadcolsum(.) instead of mean(.) to deal correctly with possible missings
	// centerpartial(.) returns Xp demeaned and with Xnp partialled-out.
	hdm_ss			= 1/(d.n) * quadcolsum( (centerpartial(d.X,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat):/(d.sdvecp)) :* ynull )
	hdm_ss			= sqrt(d.n)*max(abs(hdm_ss))

// NB: BCH says statistic is n*[something]/[denom]. rlasso doc says sqrt*[something].
// NB: BCH use denominator, rlasso doesn't.
	if (verbose>=1) {
		printf("\n")
		printf("hdm sup-score: %g\n",hdm_ss)
		printf("BCH sqrt(d.n)*max(abs(ScoreVec)) = hdm sup-score: %g\n",sqrt(d.n)*max(abs(ScoreVec)))
		printf("BCH (d.n)*max(abs(ScoreVec)): %g\n",(d.n)*max(abs(ScoreVec)))
		printf("BCH sup-score = sqrt(d.n)*max(abs(ScoreVec:/ScoreStd)): %g\n",supscore)
		printf("\n")
	}


// ******* sup-score crit values ********** //

// Following rlasso code:
//   object$supscore <- sqrt(n)*max(abs(colMeans(object$model*as.vector(object$dev))))
//    R <- 500
//    stat <- vector("numeric", length=R)
//    for (i in 1:R) {
//      g <- rnorm(n)
//      dev.g <- as.vector(g*object$dev)
//      mat <- object$model*dev.g
//      stat[i] <- sqrt(n)*max(abs(colMeans(mat)))
//    }
//    object$pvalue <- sum(stat>object$supscore)/R

	if (newseed>-1) {
		rseed(newseed)							//  set seed if requested
	}

	sim=J(R,1,.)
	for (j=1; j<=R; j++) {
		g=rnormal(d.n,1,0,1)
		// centerpartial(.) returns Xp demeaned and with Xnp partialled-out.
		sim[j,1]	= sqrt(d.n) * max(abs(mean((centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat):/(d.sdvecp)):*g:*ynull)))
	}
	hdm_pvalue		= sum(hdm_ss:<sim)/R

// Based on above but with BCH-type statistic instead.
	sim=J(R,1,.)
	if (dotsflag) {
		dotscmd = "_dots 0 0, title(Estimating sup-score p-value using " + strofreal(R) + " repetitions)"
		stata(dotscmd)
	}
	for (j=1; j<=R; j++) {
		if (dotsflag) {
			dotscmd = "_dots " + strofreal(j) + " 0"
			stata(dotscmd)
		}
		if (d.nclust==0) {
			// standard non-clustered case
			g=rnormal(d.n,1,0,1)
		}
		else {
			// clustered case - this is provisional!!
			// g is iid by cluster and repeated within clusters
			g=J(d.n,1,.)
			info = panelsetup(*d.clustid, 1)
			forbound = rows(info)				//  faster
			for (i=1; i<=forbound; i++) {
				gi	=rnormal(1,1,0,1)
				gi	=gi*J(info[i,2]-info[i,1]+1,1,1)
				g[info[i,1]..info[i,2],1] = gi
			}
		}
		// centerpartial(.) returns Xp demeaned and with Xnp partialled-out.
		ScoreVec		= 1/(d.n) * quadcolsum( (centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat):/(d.sdvecp)) :* (g:*ynull) )
		// ScoreStd is the same one generated earlier for the ss stat (i.e. not simulated)
		if (d.nclust) {
			// clustered case - this is provisional!!
			sim[j,1]		= (d.nclust)*max(abs(ScoreVec:/ScoreStd))
		}
		else if (hetero) {
			sim[j,1]		= (d.n)*max(abs(ScoreVec:/ScoreStd))
		}
		else {
			sim[j,1]		= (d.n)*max(abs(ScoreVec))
		}
	}
	supscore_pvalue		= sum(supscore:<sim)/R

	// BCH asymptotic bound, p. 20.
	// Also Matlab code (see above; note c=1.1 is in denom of SupScore:
	// ind = SupScore <= norminv(1-.05/(2*size(Mz,2)));
	if (d.nclust) {
		// clustered case - this is provisional!!
		supscore_critvalue	=c*sqrt(d.nclust)*invnormal(1-gamma/(2*((d.p)-pminus)))
	}
	else {
		supscore_critvalue	=c*sqrt(d.n)*invnormal(1-gamma/(2*((d.p)-pminus)))
	}

	if (verbose>=1) {
		printf("hdm p-value: %g\n",hdm_pvalue)
		printf("BCH p-value: %g\n",supscore_pvalue)
		printf("BCH %g percent critical value: %g\n",100*gamma,supscore_critvalue)
		printf("\n")
	}

	// BCH stat, p-value, crit-value, signif; rlasso stat, p-value
	res = (supscore, supscore_pvalue, supscore_critvalue, gamma, hdm_ss, hdm_pvalue)
	return(res)
}


real rowvector MakeSqrtLassoWeights(struct dataStruct scalar d,	///
									real colvector v,			///
									real scalar rmse,			///
									|							///
									real scalar nclust1,		///
									real scalar center			///
									)
{
	if (args()==3) {
		nclust1 = 0									//  default value of 0 serves as boolean below
		center = 0
	}

	// presumes heteroskedasticity or clustering; in homoskedastic case, Ups is just a vector of 1s

	real rowvector Ups
	
	if (d.nclust==0) {											// no cluster dependence
		// create blank full-size vector of zeros to be populated.
		St = J(1,cols(*d.X),0)
		if (center) {										// we have notpen vars, center
			// centerpartial(.) returns Xp demeaned and with Xnp partialled-out.
			centervec =	mean(centerpartial(d.Xp,d.mvecp,cons,d.Xnp,d.mvecnp,d.pihat) :* v)
			// populate correct columns, leaving other as zeros (unpenalized).
			St[1,(d.selindXp)] = quadcolsum( ((centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat) :* v) :- centervec ):^2 )/d.n
		}
		else {													// we have notpen vars, no center
			// populate correct columns, leaving other as zeros (unpenalized).
			St[1,(d.selindXp)] = quadcross((centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat)):^2,v:^2)'/d.n
		}
		// in sqrt-lasso, Ups must include rmse
		Ups = sqrt(St) * 1/rmse
	}
	else {														// cluster dependence
		// this code is memory-intensive because St is a matrix that is the same dimension as X
		// should recode to at least loop over clusters
		info = panelsetup(*d.clustid, 1)
		// First create blank full-size matrix of zeros to be populated.
		St = J(rows(*d.X),cols(*d.X),0)
		// Now populate correct columns, leaving other as zeros (unpenalized).
		// centerpartial(.) returns Xp demeaned and with Xnp partialled-out.
		St[.,(d.selindXp)] = centerpartial(d.Xp,d.mvecp,d.cons,d.Xnp,d.mvecnp,d.pihat):*v
		UpsTemp = J(d.nclust,d.p,0)
		forbound = d.nclust								//  faster
		for (i=1; i<=forbound; i++) {	
			Sti = panelsubmatrix(St, i, info)
			UpsTemp[i,.] = quadcolsum(Sti)
		}
		//  center UpsTemp
		if (center) {
			UpsTemp = UpsTemp - J(d.nclust,1,1)*quadcolsum(UpsTemp)/d.nclust
		}
		//  Default is simply to divide by nobs = nclust*T in a balanced panel.
		//  A finite-sample adjustment as in CBH's lassoCluster is to divide by (nclust-1)*T.
		//  In an unbalanced panel, we achieve this with 1/nobs * nclust/(nclust-1).
		//  In a balanced panel, 1/nobs * nclust/(nclust-1) = 1/(nclust*T) * nclust/(nclust-1) = 1/((nclust-1)*T).
		if (nclust1) {											//  override default: divide by (nclust-1)*T
			Ups = sqrt(quadcolsum(UpsTemp:^2)/(d.n)*(d.nclust)/(d.nclust-1))
		}
		else {													//  default: simply divide by n = nobs*T
			Ups = sqrt(quadcolsum(UpsTemp:^2)/(d.n))
		}
		//  in sqrt-lasso, Ups must include rmse
		Ups = Ups * 1/rmse
	}
	// =max(Ups,1), see Alg 1 in Ann of Stat
	// but since X isn't pre-standardized, need to compare to standardization vector instead of unit vector
	// nb: this the the std vector after partialling and with zeros in place of the zero-pen variables
	Ups = (Ups :< d.sdvecpnp):*d.sdvecpnp + (Ups :>= d.sdvecpnp):*Ups

	return(Ups)
}
// end MakeSqrtLassoWeights


struct outputStruct scalar RSqrtLasso(						/// Mata code for BCH sqrt rlasso
							struct dataStruct scalar d,		/// data
							real scalar hetero,				/// homosk or heteroskedasticity?
							real scalar xdep,				/// X-dependent or independent?
							real scalar R,					/// number of simulations with xdep
							real scalar lassoUps,			/// use lasso or post-lasso residuals for estimating penalty loadings?
							real scalar optTol,				///
							real scalar maxIter,			///
							real scalar zeroTol,			///
							real scalar maxUpsIter,			///
							real scalar UpsTol,				///
							real scalar verbose,			///
							real scalar c,					///
							real scalar gamma,				///
							real scalar gammad,				///
							real scalar lambda0,			///
							real scalar lalt, 				///
							real scalar corrnumber,			///
							real scalar pminus,				///
							real scalar nclust1,			/// use #nclust-1 instead of #nclust in cluster-lasso
							real scalar center,				/// center x_i*e_i or x_ij*e_ij in cluster-lasso
							real scalar supscoreflag,		///
							real scalar ssnumsim,			///
							real scalar newseed,			///
							real scalar dotsflag			///
							)
{
	struct outputStruct scalar betas

	if (gammad<=0) {										//  not user-provided, so set here
		if (d.nclust==0)	gammad=log(d.n)					//  not cluster-lasso so #obs=n
		else				gammad=log(d.nclust)			//  cluster-lasso so #obs=nclust
	}

	// in sqrt-lasso homoskedastic case, lambda does not incorporate the estimated error variance
	// so lambda0 returned by lambdaCalc is the sqrt-lasso lambda * 2 (lasso has factor of 2 built-in)
	// in heteroskedastic case, rmse is needed but only for penalty loadings;
	// interpretation is still that lambda0 is same thing as lambda
	
	// initial residuals needed for het, clustered or xdep cases
	if (xdep | hetero | d.nclust) {
		v = InitialResiduals(d,corrnumber)
		s1 = sqrt(mean(v:^2))
	}
	else {
		v = .
		s1 = .
	}
	
	if (lambda0==0) {
		// optimal lambda. divide by 2 since "sqrt-lambda"*2="lasso-lambda". 
	 	lambda=lambdaCalc(d,pminus,gamma,gammad,c,R,hetero,xdep,v,s1,lalt,newseed,dotsflag) / 2
	 }
	 else {
	 	// user-supplied lambda0
	 	lambda=lambda0
	 }

// MS:
// From Belloni web page; alpha=0.05.
// AA: Ann of Stat paper recommends 0.05/log(n). Maybe use that?
//	lambda = c*sqrt(n)*invnormal(1-0.05/(2*p))

	iter	= 0			// initialize here so that something is always saved in nupsiter

	if ((hetero==0) & (d.nclust==0)) {
		// homoskedasticity

		// Use sdvec rather than unit vector since raw data are not pre-standardized.
		// NB: this is the sdvec after partialling and with zeros for unpenalized vars.
		Ups0 = d.sdvecpnp
		betas = DoSqrtLasso(d, Ups0, lambda, verbose, optTol, maxIter, zeroTol)
	
	}
	else {
		// heteroskedasticity or clustering
		
		// initialize Ups
		Ups1 = MakeSqrtLassoWeights(d,v,s1,nclust1,center)
		// Ups1 = colmax(abs(X))  // Initial loadings suggested in Alg 1 in Ann of Stat (2014)

		do {								// start loop; always iterated at least once
			// increment counter
			iter = iter + 1
			Ups0=Ups1
			s0=s1
			if (verbose>=1) {
				printf("Estimation of penalty level/loadings: Step %f.\n",iter)
				printf("RMSE: %f\n",s0)
			}
			// lasso estimation for given penalty level & loadings; betas based on Ups0.
			betas = DoSqrtLasso(d, Ups0, lambda, verbose, optTol, maxIter, zeroTol)

			// obtain residuals; based on betas(Ups0)
			if (lassoUps==1) {
				v = betas.v		// v = (*d.y) - select(*d.X,betas.index')*betas.beta					// lasso residuals
				s1 = betas.rmse 
			}
			else {
				v = betas.vPL	// v = (*d.y) - select(*d.X,betas.index')*betas.betaPL	 				// post-lasso residuals
				s1 = betas.rmsePL
			}
			// s1 = rmse based on new residuals v(betas(Ups0)); use use large-sample definition of variance

			// change in RMSE				
			Delta = abs(s1-s0)

			// Reporting
			// unclear why cols need to be >1. Need c>1 even if r==0; c==0 and r==0 makes invtokens(.) crash.
			if (verbose>=1 & cols(betas.nameXSel)>1) {
				printf("Selected variables: %s\n",invtokens(betas.nameXSel'))
			}
			if (verbose>=1) {
				printf("Change in RMSE: %f\n\n",Delta)
			}

			// update Ups1 based on new residuals v, but only if loop is going to repeat
			// also update lambda (xdep+het case only), again only if loop is going to repeat
			if ((iter < maxUpsIter) & (Delta > UpsTol)) {
			
				// xdep lambda update
				if ((hetero==1) & (xdep==1)) {
					lambda=lambdaCalc(d,pminus,gamma,gammad,c,R,hetero,xdep,v,s1,lalt,newseed,dotsflag) / 2
				}
			
				// Ups update
				Ups1 = MakeSqrtLassoWeights(d,v,s1,nclust1,center)
			}
			
		} while ((iter < maxUpsIter) & (Delta > UpsTol))
		// When loop completes, we have: betas(Ups0); s1 based on v(betas(Ups0)); betas, rmse, Ups all consistent.
		// Note that when loop completes, we have avoided updating Ups1 based on v(betas(Ups0)).

		if (verbose>=1) {
			printf("Number of penalty loading iterations: %g\n",iter)
			if (iter == maxUpsIter) {
				printf("Warning: reached max penalty loading iterations w/o achieving convergence.\n")
			}
			else {
				printf("Penalty loadings (upsilon) convergence achieved.\n")
			}
		}

	}
	
	// sup-score stat
	if (supscoreflag) {
		betas.supscore	= doSupScore(d, c, gamma, pminus, hetero, verbose, ssnumsim, newseed, dotsflag)
	}

	// Misc
	betas.n			= d.n
	betas.nclust	= d.nclust
	betas.nupsiter	= iter
	betas.sUps		= betas.Ups :/ d.sdvecpnp		//  should be =1 under homosk.
	betas.sUps		= editmissing(betas.sUps,0)		//  in case any unpenalized (div by 0)
	betas.lambda0	= lambda						//  no diff between lambda and lambda0 with sqrt lasso
	betas.slambda	= lambda						//  no diff between lambda and std lambda with sqrt lasso
	betas.c			= c
	betas.gamma		= gamma
	betas.gammad	= gammad
	
	return(betas)
	
}
// end RSqrtLasso	

void EstimateLassoPath(							///  Complete Mata code for lassopath
				string scalar nameY,			///
				string scalar nameX,			///
				string scalar nameX_o,			///
				string scalar notpen_o,			///
				string scalar notpen_t,			///
				string scalar toest,			///
				string scalar holdout, 			/// validation data
				real scalar cons,				///
				real scalar dmflag,				///
				string scalar lambdamat,		/// single, list or missing (=> construct default list)
				real scalar lmax, 				///
				real scalar lcount,				///
				real scalar lminratio,			///
				real scalar lglmnet,			///
				string scalar UpsMat,			/// 
				string scalar stdymat,			/// 
				string scalar stdxmat,			///
				real scalar stdl, 				/// standardisation loadings?
				real scalar sqrtflag,			/// lasso or sqrt-lasso?
				real scalar alpha,				///
				real scalar post, 				///
				real scalar optTol,				///
				real scalar maxIter,			///
				real scalar zeroTol,			///
				real scalar verbose,			///
				real scalar stdcoef,			///
				real scalar noic 				///
				)
{

	struct dataStruct scalar d
	d = MakeData(nameY,nameX,nameX_o,toest,cons,dmflag,stdymat,stdxmat)
	p = cols(d.sdvec)
	
	// Estimation accommodates pre-standardized data and standardization on-the-fly.
	// Pre-standardized: lambdas and penalty loadings the same for L1 and L2 norms.
	// Standardization on-the-fly: standardization included in penalty loadings;
	//   L1 and L2 norm lambdas and penalties differ.
	
	if (UpsMat!="") {				//  overall pen loadings supplied
		Ups = st_matrix(UpsMat)
		Ups2 = Ups
	}
	else if (stdl) {				//  std loadings - standardize on the fly
		Ups = d.sdvec				//  L1 norm loadings - SDs
		Ups2 = d.varvec				//  L2 norm loadings - variances
	}
	else {
		Ups = J(1,p,1)				//  default is 1
		Ups2 = Ups
	}
	
	//  need to set loadings of notpen vars = 0
	if (notpen_o~="") {	
		npnames=tokens(notpen_o)
		forbound = cols(npnames)	//  faster
		for (i=1; i<=forbound; i++) {
				Ups		=Ups  :* (1:-(d.nameX_o:==npnames[1,i]))
				Ups2	=Ups2 :* (1:-(d.nameX_o:==npnames[1,i]))
		}
	}
	
	if (lambdamat=="") {
		if (lmax<=0) { // no lambda max given
			if (sqrtflag) {
				lmax = max(abs((d.Xy)):/((Ups)')) 
			}
			else {
				// see Friedman et al (J of Stats Software, 2010)  
				lmax = max(abs((d.Xy)):/((Ups)'))*2/max((0.001,alpha)) 
			}
		}
		lmin = lminratio*lmax
		lambda=exp(rangen(log(lmax),log(lmin),lcount))'
		lambda2=lambda
	} 
	else {
		lambda=st_matrix(lambdamat)
		lambda2=lambda
		if ((d.prestdflag) & (!sqrtflag)) {		//  data have been pre-standardized, so adjust lambdas accordingly
			lambda	=lambda * 1/(d.ysd)
		}
		if (lglmnet) {
			lambda=lambda*2*d.n
		}
	}


	if ((cols(lambda)==1) & (!hasmissing(lambda))) {				//  one lambda
		struct outputStruct scalar OUT
		if (sqrtflag) {
			OUT = DoSqrtLasso(d,Ups,lambda,verbose,optTol,maxIter,zeroTol)
		}
		else {
			OUT = DoLasso(d,Ups,Ups2,lambda,lambda2,verbose,optTol,maxIter,zeroTol,alpha,lglmnet)
		}
		ReturnResults(OUT,d,sqrtflag,stdcoef)
	}
	else if ((cols(lambda)>1) & (!hasmissing(lambda))) {		//  lambda is a vector or missing (=> default list)
		struct outputStructPath scalar OUTPATH
		if (sqrtflag) {
			OUTPATH = DoSqrtLassoPath(d,Ups,lambda,post,verbose,optTol,maxIter,zeroTol)
		}
		else {
			OUTPATH = DoLassoPath(d,Ups,Ups2,lambda,lambda2,post,verbose,optTol,maxIter,zeroTol,alpha,lglmnet,noic)
		}
		ReturnResultsPath(OUTPATH,d,nameX_o,sqrtflag,stdcoef)
		if (holdout!="") { // used for cross-validation
			getMSPE(OUTPATH,nameY,nameX,holdout,d.prestdflag,d.ysd)  
		}
		else if (!noic) { // don't calculate IC 
			getInfoCriteria(OUTPATH,d,sqrtflag)
		}
	}
}
// end EstimateLassoPath

struct outputStruct scalar DoLasso(								///
							struct dataStruct scalar d,			/// data (y,X)
							real rowvector Ups,					///	penalty loading vector (L1 norm)
							real rowvector Ups2,				///	penalty loading vector (L2 norm)
							real scalar lambda,					/// lambda (single value) (L1 norm)
							real scalar lambda2,				/// lambda (single value) (L2 norm)
							real scalar verbose,				/// reporting
							real scalar optTol,					/// convergence of beta estimates
							real scalar maxIter,				/// max number of shooting iterations
							real scalar zeroTol, 				/// tolerance to set coefficient estimates to zero
							real scalar alpha, 					/// elastic net parameter
							| real scalar lglmnet					/// 1= use lambda definition from glmnet
							)
{

	struct outputStruct scalar t

	if (args()<=10) lglmnet = 0

	p = cols(*d.X)
	n = rows(*d.X)
	
	XX = d.XX
	Xy = d.Xy

	if (verbose>=1) {
		avgUps=sum(abs(Ups))/sum(Ups:>0)
		printf("Lambda: %f\nAverage abs. loadings: %f\n", lambda,avgUps)
	}
	
	beta_init=lusolve(XX+lambda2/2*diag(Ups2),Xy)
	beta=beta_init

	if (verbose==2){
		w_old = beta
		printf("%8s %8s %10s %14s %14s\n","iter","shoots","n(w)","n(step)","f(w)")
		k=1
		wp = beta
	}

	m=0
	XX2=XX*2
	Xy2=Xy*2

	// Separate blocks for lasso, ridge, elastic net.
	// Separate lambdas and penalty loadings for L1 and L2 norms to accommodate standardization.
	// If data were pre-standardized, then lambda=lambda2 and Ups=Ups2.
	// If standardization is on-the-fly and incorporated in penalty loadings,
	// then lambdas and penalty loadings for L1 and L2 norms are different.
	while (m < maxIter)
	{
		beta_old = beta
		for (j = 1;j<=p;j++)
		{
			S0 = quadcolsum(XX2[j,.]*beta) - XX2[j,j]*beta[j] - Xy2[j]

			if (alpha==1) {					//  lasso
				if (S0 > lambda*Ups[j])
				{
					// beta[j,1] = (lambda*Ups[1,j] - S0)/(XX2[j,j])
					beta[j] = (lambda*Ups[j] - S0)/(XX2[j,j])
				}
				else if (S0 < -lambda*Ups[j])	
				{
					// beta[j,1] = (-lambda*Ups[1,j] - S0)/(XX2[j,j]) 
					beta[j] = (-lambda*Ups[j] - S0)/(XX2[j,j]) 
				}
				else 
				{
					// beta[j,1] = 0
					beta[j] = 0
				}
			}								//  end lasso
			else if (alpha>0) {				//  elastic net
				if (S0 > lambda*Ups[j]*alpha)
				{
					// beta[j,1] = (lambda*Ups[1,j]*alpha - S0)/(XX2[j,j] + lambda2*Ups2[1,j]*(1-alpha))
					beta[j] = (lambda*Ups[j]*alpha - S0)/(XX2[j,j] + lambda2*Ups2[j]*(1-alpha))
				}
				else if (S0 < -lambda*Ups[j]*alpha)
				{
					// beta[j,1] = (-lambda*Ups[1,j]*alpha - S0)/(XX2[j,j] + lambda2*Ups2[1,j]*(1-alpha))
					beta[j] = (-lambda*Ups[j]*alpha - S0)/(XX2[j,j] + lambda2*Ups2[j]*(1-alpha))
				}
				else
				{
					// beta[j,1] = 0
					beta[j] = 0
				}
			}								//  end elastic net
			else if ((alpha==0) & (1)) {	//  ridge (DISABLED)
			
				//		shooting not required for ridge since closed-formed solution exists
				// 		and is equal to initial beta.
	
				if (S0 > 0)
				{
					// beta[j,1] = (-S0)/(XX2[j,j] + lambda2*Ups2[1,j])
					beta[j] = (-S0)/(XX2[j,j] + lambda2*Ups2[j])
				}
				else if (S0 < 0)	
				{
					// beta[j,1] = (-S0)/(XX2[j,j] + lambda2*Ups2[1,j]) 
					beta[j] = (-S0)/(XX2[j,j] + lambda2*Ups2[j]) 
				}
				else 
				{
					// beta[j,1] = 0
					beta[j] = 0
				}
			}								//  end ridge
		}									//  end j loop over components of beta

       	m++

		if (verbose==2)
			{
				printf("%8.0g %8.0g %14.8e %14.8e %14.8e %14.8e\n",				///
					m,m*p,														///
					colsum(abs(beta)),											///
					colsum(abs(beta-w_old)),									///
					colsum(((*d.X)*beta-(*d.y)):^2)+lambda*colsum(abs(beta)),	///
					beta[1,1]	///
					)
				w_old = beta
				k=k+1
				wp =(wp, beta)
			}
    
		if (quadcolsum(abs(beta-beta_old))<optTol) break
	}
	
	if (verbose>=1)
	{
		printf("Number of iterations: %g\nTotal Shoots: %g\n",m,m*p)
		if (m == maxIter) {
			printf("Warning: reached max shooting iterations w/o achieving convergence.\n")
		}
		else {
			printf("Convergence achieved.\n")
		}
	}
	
	if (verbose==2) {
		printf("Initial beta and beta after estimation:\n")
		(beta_init'\beta')
	}	
	
	// convert lambda
	if (lglmnet) {
		lambda=lambda/2/n
	}
	
	// save results in t struct
	// following code should be the same in DoLasso() and DoSqrtLasso()
	
	t.niter = m
	
	// full vector
	t.betaAll = beta
	
	// compare initial beta vs estimated beta
	//(beta,beta_init)

	t.index = abs(beta) :> zeroTol
	s = sum(abs(beta) :> zeroTol) // number of selected vars, =0 if no var selected
	
	// reduce beta vector to only non-zero coeffs
	if (s>0) {
		t.beta = select(beta,t.index)
	}
	else {
		t.beta = .
	}
	
	// obtain post-OLS estimates
	if ((s>0) & (s<n) & (d.dmflag | (d.cons==0))) {
		// data are zero mean or there is no constant
		betaPL			= qrsolve(select(*d.X,t.index'),*d.y)
	}
	else if ((s>0) & (s<n)) {
		betaPL			= qrsolve((select(*d.X,t.index'):-select(d.mvec,t.index')),(*d.y:-d.ymvec))
	}
	else if (s>0) {
		betaPL			= J(s,1,.)		// set post-OLS vector = missing if s-hat > n. 
	}
	else if (s==0) {
		betaPL			= .
	}
	t.betaPL = betaPL

	// obtain intercept
	if ((d.dmflag) | (d.cons==0)) {
		// data are zero mean or there is no constant
		t.intercept 	= 0
		t.interceptPL	= 0
	}
	else if (s>0) {
		t.intercept		= mean(*d.y) - mean(select(*d.X,t.index'))*t.beta	//  something selected; obtain constant
		t.interceptPL	= mean(*d.y) - mean(select(*d.X,t.index'))*t.betaPL	// obtain constant
	}
	else {																	//  nothing selected; intercept is mean(y)
		t.intercept		= mean(*d.y)
		t.interceptPL	= mean(*d.y)
	}
	
	// "All" version of PL coefs
	t.betaAllPL = J(rows(beta),1,0)
	// need to look out for missing betaPL; if s=0 betaAllPL will just be a vector of zeros
	if (s>0) {
		t.betaAllPL[selectindex(t.index),1] = betaPL
	}

	// other objects
	if (s>0) {
		t.nameXSel	= select(d.nameX_o',t.index)
	}
	else {
		t.nameXSel	=""
	}
	t.Ups		= Ups
	t.lambda	= lambda
	t.cons		= d.cons
	t.n			= n
	t.beta_init	= beta_init
	t.s 		= s

	// obtain residuals
	if ((s>0) & (d.dmflag)) {
		// data are zero mean
		t.v		= *d.y - select(*d.X,t.index')*t.beta
		t.vPL	= *d.y - select(*d.X,t.index')*t.betaPL
	}
	else if (s>0) {
		t.v		= (*d.y:-d.ymvec) - (select(*d.X,t.index'):-select(d.mvec,t.index'))*t.beta
		t.vPL	= (*d.y:-d.ymvec) - (select(*d.X,t.index'):-select(d.mvec,t.index'))*t.betaPL
	}
	else if (d.dmflag) {
		// data are zero mean; nothing selected; residual is just y
		t.v		= *d.y
		t.vPL	= *d.y
	}
	else {
		// nothing selected; residual is demeaned y
		t.v		= (*d.y:-d.ymvec) 
		t.vPL	= (*d.y:-d.ymvec)
	}
	
	// RMSE
	t.rmse		=sqrt(mean(t.v:^2))
	t.rmsePL	=sqrt(mean(t.vPL:^2))
	
	return(t)
}
// end DoLasso

struct outputStruct scalar DoSqrtLasso(							///
							struct dataStruct scalar d,			/// data (y,X)
							real rowvector Ups,					///	penalty loading vector
							real scalar lambda,					/// lambda (single value)
							real scalar verbose,				/// reporting
							real scalar optTol,					/// convergence of beta estimates
							real scalar maxIter,				/// max number of shooting iterations
							real scalar zeroTol					/// tolerance to set coefficient estimates to zero
							)
{

	struct outputStruct scalar t

	p = cols(*d.X)
	n = rows(*d.X)
	
	XX = d.XX
	Xy = d.Xy
	
	if (verbose>=1) {
		avgUps=sum(abs(Ups))/sum(Ups:>0)
		printf("Lambda: %f\nAverage abs. loadings: %f\n", lambda,avgUps)
	}

	beta_init=lusolve(XX+lambda*diag(Ups),Xy)
	beta=beta_init

	if (verbose==2){
		w_old = beta
		printf("%8s %8s %10s %14s %14s\n","iter","shoots","n(w)","n(step)","f(w)")
		k=1
		wp = beta
	}

	m=0
	XX=XX/n
	Xy=Xy/n

	ERROR = ((*d.y):-(d.ymvec)) - ((*d.X):-(d.mvec))*beta
	Qhat = mean(ERROR:^2)

	while (m < maxIter)
	{
		beta_old = beta
		for (j = 1;j<=p;j++)
		{
			S0 = quadcolsum(XX[j,.]*beta) - XX[j,j]*beta[j] - Xy[j]

			if ( abs(beta[j])>0 ) {
				ERROR = ERROR + ((*d.X)[.,j]:-(d.mvec)[j])*beta[j]
				Qhat = mean(ERROR:^2)
			}
			
			if ( n^2 < (lambda * Ups[j])^2 / XX[j,j]) {
				beta[j] = 0
			}
			else if (S0 > lambda/n*Ups[j]*sqrt(Qhat))
			{
				beta[j]= ( ( lambda * Ups[j] / sqrt(n^2 - (lambda * Ups[j])^2 / XX[j,j] ) ) * sqrt(max((Qhat-(S0^2/XX[j,j]),0)))-S0)  / XX[j,j]
				ERROR = ERROR - ((*d.X)[.,j]:-(d.mvec)[j])*beta[j]
			}
			else if (S0 < -lambda/n*Ups[j]*sqrt(Qhat))	
       		{
				beta[j]= ( - ( lambda * Ups[j] / sqrt(n^2 - (lambda * Ups[j])^2 / XX[j,j] ) ) * sqrt(max((Qhat-(S0^2/XX[j,j]),0)))-S0)  / XX[j,j]
				ERROR = ERROR - ((*d.X)[.,j]:-(d.mvec)[j])*beta[j]
           	}
       		else 
			{
       			beta[j] = 0
			}
		}

		ERRnorm=sqrt( quadcolsum( (((*d.y):-(d.ymvec))-((*d.X):-(d.mvec))*beta):^2 ) )
		fobj = ERRnorm/sqrt(n) + (lambda/n)*Ups*abs(beta)
		
		if (ERRnorm>1e-10) {
			aaa = (sqrt(n)*ERROR/ERRnorm)
			dual = aaa'((*d.y):-(d.ymvec))/n  - abs(lambda/n*Ups' - abs(((*d.X):-(d.mvec))'aaa/n))'abs(beta)
		}
		else {
			dual = (lambda/n)*Ups*abs(beta)
		}
		
       	m++

		if (verbose==2) {
				printf("%8.0g %8.0g %14.8e %14.8e %14.8e\n",m,m*p,colsum(abs(beta)),colsum(abs(beta)),colsum((((*d.X):-(d.mvec))*beta-((*d.y):-(d.ymvec))):^2)+lambda*colsum(abs(beta)))
				w_old = beta
				k=k+1
				wp =(wp, beta)
		}

		if (quadcolsum(abs(beta-beta_old))<optTol) {
			if (fobj - dual < 1e-6) {
				break
			}
		}
	}
	
	if (verbose>=1)
	{
		printf("Number of iterations: %g\nTotal Shoots: %g\n",m,m*p)
		if (m == maxIter) {
			printf("Warning: reached max shooting iterations w/o achieving convergence.\n")
		}
		else {
			printf("Convergence achieved.\n")
		}
	}

	// save results in t structure
	// following code should be the same in DoLasso() and DoSqrtLasso()
	
	t.niter = m
	
	// full vector
	t.betaAll = beta

	t.index = abs(beta) :> zeroTol
	s = sum(abs(beta) :> zeroTol) // number of selected vars, =0 if no var selected
	
	// reduce beta vector to only non-zero coeffs
	if (s>0) {
		t.beta = select(beta,t.index)
	}
	else {
		t.beta = .
	}
	
	// obtain post-OLS estimates
	if ((s>0) & (s<n) & (d.dmflag | (d.cons==0))) {
		// data are zero mean or there is no constant
		betaPL			= qrsolve(select(*d.X,t.index'),*d.y)
	}
	else if ((s>0) & (s<n)) {
		betaPL			= qrsolve((select(*d.X,t.index'):-select(d.mvec,t.index')),(*d.y:-d.ymvec))
	}
	else if (s>0) {
		betaPL			= J(s,1,.)		// set post-OLS vector = missing if s-hat > n. 
	}
	else if (s==0) {
		betaPL			= .
	}
	t.betaPL = betaPL
	
	// obtain intercept
	if ((d.dmflag) | (d.cons==0)) {
		// data are zero mean or there is no constant
		t.intercept 	= 0
		t.interceptPL	= 0
	}
	else if (s>0) {
		t.intercept		= mean(*d.y) - mean(select(*d.X,t.index'))*t.beta	//  something selected; obtain constant
		t.interceptPL	= mean(*d.y) - mean(select(*d.X,t.index'))*t.betaPL	 // obtain constant
	}
	else {																	//  nothing selected; intercept is mean(y)
		t.intercept		= mean(*d.y)
		t.interceptPL	= mean(*d.y)
	}
	
	// "All" version of PL coefs
	t.betaAllPL = J(rows(beta),1,0)
	// need to look out for missing betaPL; if s=0 betaAllPL will just be a vector of zeros
	if (s>0) {
		t.betaAllPL[selectindex(t.index),1] = betaPL
	}

	// other objects
	if (s>0) {
		t.nameXSel	= select(d.nameX_o',t.index)
	}
	else {
		t.nameXSel	=""
	}
	t.Ups		= Ups
	t.lambda	= lambda
	t.cons		= d.cons
	t.n			= n
	t.beta_init	= beta_init
	t.s 		= s

	// obtain residuals
	if ((s>0) & (d.dmflag)) {
		// data are zero mean
		t.v		= *d.y - select(*d.X,t.index')*t.beta
		t.vPL	= *d.y - select(*d.X,t.index')*t.betaPL
	}
	else if (s>0) {
		t.v		= (*d.y:-d.ymvec) - (select(*d.X,t.index'):-select(d.mvec,t.index'))*t.beta
		t.vPL	= (*d.y:-d.ymvec) - (select(*d.X,t.index'):-select(d.mvec,t.index'))*t.betaPL
	}
	else if (d.dmflag) {
		// data are zero mean; nothing selected; residual is just y
		t.v		= *d.y
		t.vPL	= *d.y
	}
	else {
		// nothing selected; residual is demeaned y
		t.v		= (*d.y:-d.ymvec) 
		t.vPL	= (*d.y:-d.ymvec)
	}
	
	// RMSE
	t.rmse		=sqrt(mean(t.v:^2))
	t.rmsePL	=sqrt(mean(t.vPL:^2))	

	return(t)

}
// end DoSqrtLasso


struct outputStructPath scalar DoSqrtLassoPath(	struct dataStruct scalar d,
												real rowvector Ups,
												real rowvector lvec,
												real scalar post, 
												real scalar verbose, 
												real scalar optTol,
												real scalar maxIter,
												real scalar zeroTol)
{

		struct outputStructPath scalar t

		p = cols(*d.X)
		n = rows(*d.X)
		
		XX = d.XX
		Xy = d.Xy
		
		lmax=max(lvec)
		lcount=cols(lvec)
		beta=lusolve(XX+lmax*diag(Ups),Xy) // beta start	
		
		XX=XX/n
		Xy=Xy/n

		lpath = J(lcount,p,.) // create empty matrix which stores coef path
		for (k = 1;k<=lcount;k++) { // loop over lambda
		
			lambda=lvec[1,k]
			
			//beta=lusolve(XX*n+lmax/2*diag(Ups),Xy*n)	
			// more stable than only one initial beta at the top
				
			m=0
				
			ERROR = ((*d.y):-(d.ymvec)) - ((*d.X):-(d.mvec))*beta
			Qhat = mean(ERROR:^2)

			while (m < maxIter)
			{
				beta_old = beta
				for (j = 1;j<=p;j++)
				{
					S0 = quadcolsum(XX[j,.]*beta) - XX[j,j]*beta[j] - Xy[j]

					if ( abs(beta[j])>0 ) {
						ERROR = ERROR + ((*d.X)[.,j]:-(d.mvec)[j])*beta[j]
						Qhat = mean(ERROR:^2)
					}
					
					if ( n^2 < (lambda * Ups[j])^2 / XX[j,j]) {
						beta[j] = 0
					}
					else if (S0 > lambda/n*Ups[j]*sqrt(Qhat))
					{
						beta[j]= ( ( lambda * Ups[j] / sqrt(n^2 - (lambda * Ups[j])^2 / XX[j,j] ) ) * sqrt(max((Qhat-(S0^2/XX[j,j]),0)))-S0)  / XX[j,j]
						ERROR = ERROR - ((*d.X)[.,j]:-(d.mvec)[j])*beta[j]
					}
					else if (S0 < -lambda/n*Ups[j]*sqrt(Qhat))	
					{
						beta[j]= ( - ( lambda * Ups[j] / sqrt(n^2 - (lambda * Ups[j])^2 / XX[j,j] ) ) * sqrt(max((Qhat-(S0^2/XX[j,j]),0)))-S0)  / XX[j,j]
						ERROR = ERROR - ((*d.X)[.,j]:-(d.mvec)[j])*beta[j]
					}
					else 
					{
						beta[j] = 0
					}
				}

				ERRnorm=sqrt( quadcolsum( (((*d.y):-(d.ymvec))-((*d.X):-(d.mvec))*beta):^2 ) )
				fobj = ERRnorm/sqrt(n) + (lambda/n)*Ups*abs(beta)
				
				if (ERRnorm>1e-10) {
					aaa = (sqrt(n)*ERROR/ERRnorm)
					dual = aaa'((*d.y):-(d.ymvec))/n  - abs(lambda/n*Ups' - abs(((*d.X):-(d.mvec))'aaa/n))'abs(beta)
				}
				else {
					dual = (lambda/n)*Ups*abs(beta)
				}
				
				m++

				if (quadcolsum(abs(beta-beta_old))<optTol) { 
					if (fobj - dual < 1e-6) {
						break
					}
				}
				//terminate = (quadcolsum(abs(beta-beta_old))<optTol)*((fobj - dual) < 1e-6)
			}
			lpath[k,.]=beta'
		}
		
		// following code should be the same for DoLassoPath() and DoSqrtLassoPath()
		
		lpath=edittozerotol(lpath, zeroTol)
		
		if (post) { 
			betasPL = J(lcount,p,0)
			nonzero0 = J(1,p,0)
			for (k = 1;k<=lcount;k++) { // loop over lambda points
				nonzero = lpath[k,.]:!=0  // 0-1 vector
				sk = sum(nonzero)			
				if ((0<sk) & (sk<n)) { // only if 0<s<n
					if ((nonzero0==nonzero) & (k>=2)) { // no change in active set
						betasPL[k,.] = betasPL[k-1,.]
					}
					else { 
						ix = selectindex(nonzero)	// index of non-zeros
						// obtain beta-hat
						if (d.dmflag) {
							// data are mean zero
							betak=qrsolve(select((*d.X),nonzero),(*d.y))
						}
						else {
							betak=qrsolve(select((*d.X),nonzero):-select(d.mvec,nonzero),((*d.y):-d.ymvec))
						}
						betasPL[k,ix] = betak'
						nonzero0=nonzero
					}
				}
			}
			t.betas 		= betasPL	
		}
		else {
			t.betas 		= lpath
		}
		
		t.lambdalist	= lvec'
		t.Ups			= Ups
		if (d.dmflag) {
			t.intercept	= 0
		}
		else {
			t.intercept	= mean(*d.y):-mean((*d.X))*(t.betas')
		}
		t.cons 			= d.cons
		
		// sqrt-lasso df
		t.shat	= quadrowsum(t.betas:!=0) :+ (d.cons | d.dmflag)
		t.dof 	= t.shat
				
		return(t)		
		
}
// end DoSqrtLassoPath


void ReturnResults(		struct outputStruct scalar t,		///
						struct dataStruct scalar d,			///
						|									///
						real scalar sqrtflag,				///
						real scalar stdcoef					///
						)
{
	// default values
	if (args()<=2) sqrtflag = 0
	if (args()<=3) stdcoef = 0

	if (rows(t.betaAll)) {					// estimation results to insert

		// initialize from results struct
		s			= t.s					// all vars in param vector incl notpen but EXCL constant
		k			= t.s+t.cons			// all vars in param vector incl notpen PLUS constant
		betaAll		= t.betaAll
		betaAllPL	= t.betaAllPL
		Ups			= t.Ups
		eUps		= t.Ups					// rlasso only
		sUps		= t.sUps
		rmse		= t.rmse
		rmsePL		= t.rmsePL
		lambda		= t.lambda
		// initialize from data struct
		AllNames	= d.nameX_o'			// d.nameX_o is a row vector; AllNames is a col vector (to match coef vectors)
	
		// un-standardize unless overridden by stdcoef
		if (d.prestdflag & stdcoef==0) {
			betaAll		= betaAll		:/ d.sdvec' * d.ysd		//  beta is col vector, sdvec is row vector, ysd is scalar
			betaAllPL	= betaAllPL		:/ d.sdvec' * d.ysd		//  beta is col vector, sdvec is row vector, ysd is scalar
			eUps		= eUps			:* d.sdvec				//  rlasso only; eUps and sdvec both row vectors; Ups does not change
			rmse		= rmse			* d.ysd
			rmsePL		= rmsePL		* d.ysd
			if (sqrtflag==0) {									//  sqrt-lasso lambdas don't need unstandardizing
				lambda		= lambda		* d.ysd
			}
		}
		
		if (t.cons) {											//  pre-standardized means no constant
			intercept	= t.intercept
			interceptPL	= t.interceptPL
		}
		if (s>0) {												//  do here so that we select from std/unstd vector
			beta		= select(betaAll,t.index)
			betaPL		= select(betaAllPL,t.index)
			Names		= select(AllNames,t.index)
		}
		else {
			beta		= .
			betaPL		= .
			Names		= ""
		}
	
		if ((s>0) & (t.cons)) {									//  add constant to end of vectors
			beta		= (beta			\ intercept)		
			betaPL		= (betaPL		\ interceptPL)	
			betaAll		= (betaAll		\ intercept)		
			betaAllPL	= (betaAllPL	\ interceptPL)	
			NamesCons	= (Names		\ "_cons")
			AllNamesCons= (AllNames 	\ "_cons")
		}
		else if ((s>0) & (!t.cons)) {							//  coef vectors already OK, just need names with _cons
			NamesCons	= Names
			AllNamesCons= AllNames
		}
		else if ((s==0) & (t.cons)) {
			beta		= intercept						
			betaPL		= interceptPL					
			NamesCons	= "_cons"
			AllNamesCons= (AllNames 	\ "_cons")
			betaAll		= (betaAll		\ intercept)			//  will be all zeros + intercept at end
			betaAllPL	= (betaAllPL	\ interceptPL)	
		}
		else {
			beta		= .					
			betaPL		= .				
			NamesCons	= ""
			AllNamesCons= AllNames
			betaAll		= betaAll								//  will be all zeros
			betaAllPL	= betaAllPL
		}
		
		st_rclear() 
		if ((s > 0) | (t.cons)) {
			// matrix stripes
			coln		=(J(rows(Names),1,""),Names)			//  Names is a col vector so need #rows
			colnCons	=(J(rows(NamesCons),1,""),NamesCons)	//  Names is a col vector so need #rows
			// row vector of names
			st_global("r(sel)",invtokens(Names'))				//  selected vars exclude cons but here may include notpen
			// column vectors
			st_matrix("r(b)",beta)
			st_matrix("r(bOLS)",betaPL)
			st_matrixrowstripe("r(b)",colnCons)
			st_matrixrowstripe("r(bOLS)",colnCons)
		}
		
		// matrix stripe
		AllcolnCons=(J(rows(AllNamesCons),1,""),AllNamesCons)
		// column vectors
		st_matrix("r(bAll)",betaAll)
		st_matrix("r(bAllOLS)",betaAllPL)
		st_matrixrowstripe("r(bAll)",AllcolnCons)
		st_matrixrowstripe("r(bAllOLS)",AllcolnCons)
	
		// matrix stripe
		coln=(J(rows(AllNames),1,""),AllNames)
		// row vectors
		st_matrix("r(stdvec)",d.sdvec)
		st_matrix("r(Ups)",Ups)
		st_matrix("r(eUps)",eUps)						//  rlasso only
		st_matrixcolstripe("r(stdvec)",coln)
		st_matrixcolstripe("r(Ups)",coln)
		st_matrixcolstripe("r(eUps)",coln)				//  rlasso only
	
		if ((cols(sUps)>0) & (cols(sUps)<.)) {			//  "cols(t.sUps)<." may be unnecessary - if missing, cols(.) = 0.
			st_matrix("r(sUps)",sUps)
			st_matrix("r(stdvecpnp)",d.sdvecpnp)
			st_matrixcolstripe("r(sUps)",coln)
			st_matrixcolstripe("r(stdvecpnp)",coln)
		}
		
		st_matrix("r(beta_init)",t.beta_init)
	}
	
	// BCH stat, p-value, crit-value, signif; rlasso stat, p-value
	if (cols(t.supscore)) {
		st_matrix("r(supscore)",t.supscore)
		coln=(J(6,1,""),("BCH_ss" \ "BCH_p" \ "BCH_cv" \ "BCH_gamma" \ "hdm_ss" \ "hdm_p" ))
		st_matrixcolstripe("r(supscore)",coln)
	}
	
	// Can always return these; scalars will just be missing
	st_numscalar("r(lambda)",lambda)
	st_numscalar("r(rmse)",rmse)
	st_numscalar("r(rmsePL)",rmsePL)
	st_numscalar("r(s)",s)
	st_numscalar("r(k)",k)
	st_numscalar("r(lcount)",1)
	st_numscalar("r(lambda0)",t.lambda0)
	st_numscalar("r(slambda)",t.slambda)
	st_numscalar("r(c)",t.c)
	st_numscalar("r(gamma)",t.gamma)
	st_numscalar("r(gammad)",t.gammad)
	st_numscalar("r(N)",t.n)
	st_numscalar("r(N_clust)",t.nclust)
	st_numscalar("r(niter)",t.niter)
	st_numscalar("r(nupsiter)",t.nupsiter)
	
}
// end ReturnResults

real matrix getMinIC(real matrix IC)		//  =0 if nothing to partial out, =projection coefs if Zs.
{
		licid=.
		minindex(IC,1,licid,.)	// returns index of lambda that minimises ic
		if (rows(licid)>1) {    // no unique lopt 
			licid=licid[1,1] 	
			icmin=IC[1,licid]
			licunique=0
		}
		else {
			icmin=IC[1,licid]
			licunique=1
		}
	R = (licid,icmin,licunique)
	return(R)
}
// end get minimum IC

void getInfoCriteria(struct outputStructPath scalar t,
 			struct dataStruct scalar d,
			real scalar sqrtflag)
{		
		// t.betas is lcount by p	
		// t.dof is 1 by lcount

 		XB  = quadcross((*d.X)',(t.betas)') :+ t.intercept 	// n by lcount
		
 		TSS = quadcolsum(((*d.y):-(d.ymvec)):^2)	// 1 by lcount
 		ESS = quadcolsum(((XB) :-(d.ymvec)):^2)	 	
		RSS = quadcolsum(((*d.y) :-(XB)):^2)
		if ((d.prestdflag) & (!sqrtflag)) {	
			TSS=TSS*(d.ysd)^2
			ESS=ESS*(d.ysd)^2
			RSS=RSS*(d.ysd)^2
		}
		RSQ = ESS:/TSS
		//RSQ =1:-RSS:/TSS
		//(TSS:-RSS):-ESS

		// calculate aic and bic
		// For reference: in R the definitions would be
		// AIC = d.n + d.n*log(2*pi()) + d.n*log(RSS/d.n) + 2*(t.dof')
		// BIC = d.n + d.n*log(2*pi()) + d.n*log(RSS/d.n) + log(d.n)*(t.dof')
		AIC		= d.n*log(RSS/d.n) + (t.dof')*2 
		BIC 	= d.n*log(RSS/d.n) + (t.dof')*log(d.n) 
		EBIC 	= BIC :+ 2 * (t.dof') * log(d.p)
		AICC	= d.n*log(RSS/d.n) + (t.dof')*2:*((d.n):/(d.n:-t.dof'))

		// obtain minimum IC and obtimal lambda
		AICinfo = getMinIC(AIC)
		BICinfo = getMinIC(BIC)
		EBICinfo = getMinIC(EBIC)
		AICCinfo = getMinIC(AICC)
		laicid=AICinfo[1,1]
		aicmin=AICinfo[1,2]
		lbicid=BICinfo[1,1]
		bicmin=BICinfo[1,2]
		lebicid=EBICinfo[1,1]
		ebicmin=EBICinfo[1,2]
		laiccid=AICCinfo[1,1]
		aiccmin=AICCinfo[1,2]
		
		st_matrix("r(dof)",t.dof)
		st_matrixcolstripe("r(dof)",("","Df"))
		// return 
 		st_matrix("r(rsq)",RSQ')				
 		st_matrix("r(tss)",TSS')				
		st_matrix("r(ess)",ESS')				
		st_matrix("r(rss)",RSS')	
		/// aic
		st_matrix("r(aic)",AIC')
		st_numscalar("r(aicmin)",aicmin)
		st_numscalar("r(laicid)",laicid)
		/// bic
		st_matrix("r(bic)",BIC')
		st_numscalar("r(bicmin)",bicmin)
		st_numscalar("r(lbicid)",lbicid)	
		/// ebic
		st_matrix("r(ebic)",EBIC')
		st_numscalar("r(ebicmin)",ebicmin)
		st_numscalar("r(lebicid)",lebicid)	
		/// aicc
		st_matrix("r(aicc)",AICC')
		st_numscalar("r(aiccmin)",aiccmin)
		st_numscalar("r(laiccid)",laiccid)	
}
// end 


void getMSPE(struct outputStructPath scalar t,
								string scalar varY,
								string scalar varX, 
								string scalar holdout, // marks validation data set
								real scalar prestd,
								real scalar ysd)
{		
		// get beta matrix
		bhat=t.betas 		// lcount by p	
		
		// get validation data
		st_view(y0,.,varY,holdout)
		st_view(X0,.,varX,holdout) 	// n by p
		
		// predicted values
		X0B=quadcross(X0',bhat') 	// n by lcount
 		
		// add intercepts
		if (t.cons) { 	
			X0B=X0B :+ t.intercept 	// t.intercept is 1 by lcount
 		}
		
		// mean squared prediction error
		//X0
		MSPE= mean((y0:-X0B):^2) 	// 1 by lcount vector
		
		if (prestd) (
			MSPE = MSPE :* (ysd)^2
		)
		
		st_matrix("r(mspe)",MSPE)
}
// end getMSPE


// return the lambda id of the largest lambda at which the MSE
// is within one standard error of the minimal MSE.		
real scalar getOneSeLam(real rowvector mse, real rowvector sd, real scalar id) 
{
	minmse = mse[1,id] // minimal MSE
	minsd = sd[1,id]  // SE of minimal MSE
	criteria = mse[1,id]+sd[1,id] // max allowed MSE
	for (j=0; j<id; j++) {
		theid=id-j 
		thismspe= mse[1,theid]
		if (thismspe > criteria) { // if MSE is outside of interval, stop
				theid = id-j+1 // go back by one id and break
				break
		}
	} 
	return(theid)
}


// partial out program
void s_partial(	string scalar Ynames,
				string scalar Pnames,
				string scalar tYnames,
				string scalar touse,		//  indicates full sample
				string scalar toest,		//  indicates estimation subsample
				scalar dmflag,				//  indicates data are already mean zero; no constant
				scalar solver)

{

// All varnames should be basic form, no FV or TS operators etc.
// Y = structural variables
// P = variables to be partialled out
// touse = sample
// dmflag = 0 or 1
// solver = 0, 1 or 2
// Strategy is to demean (numerically more stable in case of scaling problems)
// and then use svqrsolve(.) Mata program:
//   svsolve if no collinearites (more accurate),
//   qrsolve if there are collinearities (drops columns/set coeffs to zero),
//   and svsolve if qrsolve can't find the collinearities that svsolve can.

	Ytokens=tokens(Ynames)
	Ptokens=tokens(Pnames)
	st_view(Y, ., Ytokens, touse)			//  full sample
	st_view(P, ., Ptokens, touse)
	st_view(Yest, ., Ytokens, toest)		//  estimation subsample
	st_view(Pest, ., Ptokens, toest)

	if (tYnames ~= "") {
		tYflag=1
		tYtokens=tokens(tYnames)
		st_view(tY, ., tYtokens, touse)		//  full sample
	}
	else {
		tYflag=0
	}

	L = cols(P)

	// means are based on estimation subsample
	// if dmflag=1, everything is already mean zero and there is no constant
	if ((!dmflag) & L>0) {					//  Vars to partial out including constant
		Ymeans = mean(Yest) 
		Pmeans = mean(Pest) 
	}
	else if (!dmflag) {						//  Only constant to partial out = demean
		Ymeans = mean(Yest) 
	}

	//	Partial-out coeffs.
	//	Not necessary if no vars other than constant. r=rank.
	//  coef vector b is based on estimation subsample
	if ((!dmflag) & L>0) {
		b = svqrsolve(Pest :- Pmeans, Yest :- Ymeans, r=., solver)	//  partial out P + cons
	}
	else if (L>0) {
		b = svqrsolve(Pest, Yest, r=., solver)						//  partial out P
	}
	else {
		r=1															//  partial out cons (demean)
	}

	//	Replace with residuals.
	//  Use full sample (Y, P) and not estimation subsample (Yest, Pest).
	if ((!dmflag) & L>0) {					//  Vars to partial out including constant
		if (tYflag) {
			tY[.,.] = (Y :- Ymeans) - (P :- Pmeans)*b
		}
		else {
			Y[.,.] = (Y :- Ymeans) - (P :- Pmeans)*b
		}
	}
	else if (L>0) {							//  Vars to partial out NOT including constant
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
	st_numscalar("r(rank)",r+(!dmflag))				//  Include constant in rank
}  
//end program s_partial


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


// Mata utility for standardizing; called from Stata.
void s_standard(	string scalar Xnames,				//  names of original Stata variables
					string scalar tXnames,				//  names of optional Stata temp vars to initialize/modify
					string scalar touse,				//  full sample
					string scalar toest,				//  sample on which standardization is based
					real scalar consmodel,				//  =1 if constant in model (=> data to be demeaned)
					real scalar dmflag,					//  =1 if data already demeaned
					real scalar transform				//  flag to indicate that data are to be transformed
					)
{

// All varnames should be basic form, no FV or TS operators etc.
// Can include tempvars.

	Xtokens=tokens(Xnames)
	st_view(X, ., Xtokens, touse)		// full sample
	st_view(Xest, ., Xtokens, toest)	// estimation subsample
	
	if (tXnames ~= "") {
		tXflag=1
		tXtokens=tokens(tXnames)
		st_view(tX, ., tXtokens, touse)
	}
	else {
		tXflag=0
	}

	if (dmflag) {
		// data already demeaned
		mvec			= J(1,cols(Xest),0)
		s				= sqrt(mean((Xest):^2))
		// if SD=0 (constant var), convention is to set SD to 1
		if (sum(s:==0)) {
			_editvalue(s,0,1)
		}
		if (transform) {
			if (tXflag) {
				tX[.,.]=X:/s
			}
			else {
				X[.,.]=X:/s
			}
		}
	}
	else {
		// mean and SD based on estimation subsample
		mvec			= mean(Xest)
		s				= sqrt(mean((Xest:-mvec):^2))

		// if SD=0 (constant var), convention is to set SD to 1
		if (sum(s:==0)) {
			_editvalue(s,0,1)
		}
		if (transform) {
			if ((tXflag) & (consmodel)) {
				tX[.,.]=(X:-mvec):/s			//  consmodel => demean temp vars before standardizing
			}
			else if (consmodel) {
				X[.,.]=(X:-mvec):/s				//  consmodel => demean orig vars before standardizing
			}
			else if (tXflag) {
				tX[.,.]=X:/s					//  no consmodel => just standardize temp vars
			}
			else {
				X[.,.]=X:/s						//  no consmodel => just standardize orig vars
			}
		}
	}
	
	st_matrix("r(stdvec)",s)
	st_matrix("r(mvec)",mvec)
}  
//end program s_standard

// utility for rlasso
// returns X after centering (if cons present)
// and partialling-out (if any partialled-out vars present)
real matrix centerpartial(	pointer matrix X,		///
							real rowvector Xmvec,	///
							real scalar cons,		///
							pointer matrix Z,		///
							real matrix Zmvec,		///
							real matrix pihat)		//  =0 if nothing to partial out, =projection coefs if Zs.
{
	if ((pihat==0) & (cons==0)) {
		return(*X)
	}
	else if ((pihat==0) & (cons==1)) {
		return((*X) :- Xmvec)
	}
	else {
		return(((*X):-Xmvec) - (((*Z):-Zmvec)*pihat))
	}
}
// end program centerpartial

//********************* FROM MOREMATA BY BEN JANN *********************//
// Based on:
// mm_quantile.mata
// version 1.0.8  20dec2007  Ben Jann

real matrix mm_quantile(real matrix X, | real colvector w,
 real matrix P, real scalar altdef)
{
    real rowvector result
    real scalar c, cX, cP, r, i

    if (args()<2) w = 1
    if (args()<3) P = (0, .25, .50, .75, 1)'
    if (args()<4) altdef = 0
    if (cols(X)==1 & cols(P)!=1 & rows(P)==1)
     return(mm_quantile(X, w, P', altdef)')
    if (missing(P) | missing(X) | missing(w)) _error(3351)
    if (rows(w)!=1 & rows(w)!=rows(X)) _error(3200)
    r = rows(P)
    c = max(((cX=cols(X)), (cP=cols(P))))
    if (cX!=1 & cX<c) _error(3200)
    if (cP!=1 & cP<c) _error(3200)
    if (rows(X)==0 | r==0 | c==0) return(J(r,c,.))
    if (c==1) return(_mm_quantile(X, w, P, altdef))
    result = J(r, c, .)
    if (cP==1) for (i=1; i<=c; i++)
     result[,i] = _mm_quantile(X[,i], w, P, altdef)
    else if (cX==1) for (i=1; i<=c; i++)
     result[,i] = _mm_quantile(X, w, P[,i], altdef)
    else for (i=1; i<=c; i++)
     result[,i] = _mm_quantile(X[,i], w, P[,i], altdef)
    return(result)
}

real colvector _mm_quantile(
 real colvector X,
 real colvector w,
 real colvector P,
 real scalar altdef)
{
    real colvector g, j, j1, p
    real scalar N

    if (w!=1) return(_mm_quantilew(X, w, P, altdef))
    N = rows(X)
    p = order(X,1)
    if (altdef) g = P*N + P
    else g = P*N
    j = floor(g)
    if (altdef) g = g - j
    else g = 0.5 :+ 0.5*((g - j):>0)
    j1 = j:+1
    j = j :* (j:>=1)
    _editvalue(j, 0, 1)
    j = j :* (j:<=N)
    _editvalue(j, 0, N)
    j1 = j1 :* (j1:>=1)
    _editvalue(j1, 0, 1)
    j1 = j1 :* (j1:<=N)
    _editvalue(j1, 0, N)
    return((1:-g):*X[p[j]] + g:*X[p[j1]])
}

real colvector _mm_quantilew(
 real colvector X,
 real colvector w,
 real colvector P,
 real scalar altdef)
{
    real colvector Q, pi, pj
    real scalar i, I, j, jj, J, rsum, W
    pointer scalar ww

    I  = rows(X)
    ww = (rows(w)==1 ? &J(I,1,w) : &w)
    if (altdef) return(_mm_quantilewalt(X, *ww, P))
    W  = quadsum(*ww)
    pi = order(X, 1)
    if (anyof(*ww, 0)) {
        pi = select(pi,(*ww)[pi]:!=0)
        I = rows(pi)
    }
    pj = order(P, 1)
    J  = rows(P)
    Q  = J(J, 1, .)
    j  = 1
    jj = pj[1]
    rsum = 0
    for (i=1; i<=I; i++) {
        rsum = rsum + (*ww)[pi[i]]
        if (i<I) {
            if (rsum<P[jj]*W) continue
            if (X[pi[i]]==X[pi[i+1]]) continue
        }
        while (1) {
            if (rsum>P[jj]*W | i==I) Q[jj] = X[pi[i]]
            else Q[jj] = (X[pi[i]] + X[pi[i+1]])/2
            j++
            if (j>J) break
            jj = pj[j]
            if (i<I & rsum<P[jj]*W) break
        }
        if (j>J) break
    }
    return(Q)
}

real colvector _mm_quantilewalt(
 real colvector X,
 real colvector w,
 real colvector P)
{
    real colvector Q, pi, pj
    real scalar i, I, j, jj, J, rsum, rsum0, W, ub, g

    W  = quadsum(w) + 1
    pi = order(X, 1)
    if (anyof(w, 0)) pi = select(pi, w[pi]:!=0)
    I  = rows(pi)
    pj = order(P, 1)
    J  = rows(P)
    Q  = J(J, 1, .)
    rsum = w[pi[1]]
    for (j=1; j<=J; j++) {
        jj = pj[j]
        if (P[jj]*W <= rsum) Q[jj] = X[pi[1]]
        else break
    }
    for (i=2; i<=I; i++) {
        rsum0 = rsum
        rsum = rsum + w[pi[i]]
        if (i<I & rsum < P[jj]*W) continue
        while (1) {
            ub = rsum0+1
            if (P[jj]*W>=ub | X[pi[i]]==X[pi[i-1]]) Q[jj] = X[pi[i]]
            else {
                g = (ub - P[jj]*W) / (ub - rsum0)
                Q[jj] = X[pi[i-1]]*g + X[pi[i]]*(1-g)
            }
            j++
            if (j>J) break
            jj = pj[j]
            if (i<I & rsum < P[jj]*W) break
        }
        if (j>J) break
    }
    return(Q)
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

// ftools can flip matastrict to on, causing program to fail to load.
local mstrict `c(matastrict)'
cap ftools, check
if _rc {
	// fails check, likely not installed, so do not compile
	// _fe Stata program will use (slower) Stata code
	exit
}
else {
	// temporarily set matastrict off
	qui set matastrict off
}

// Compile Mata function s_fe.
version 13
mata:

// FE transformation.
// Uses Sergio Correia's FTOOLS package - faster and does not require the data to be sorted.
void s_fe(		string scalar Xnames,
				string scalar tXnames,
				string scalar fe,
				string scalar touse,
				string scalar toest)
{

	class Factor scalar F
	F = factor(fe, touse)
	F.panelsetup()

	Xtokens=tokens(Xnames)
	X = st_data( ., Xtokens, touse)
	tXtokens=tokens(tXnames)
	st_view(tX, ., tXtokens, touse)

	w = F.sort(st_data(., toest, touse))				//  extract toest variable
	counts = panelsum(w, F.info)
	means = editmissing(panelsum(F.sort(X), w, F.info) :/ counts, 0)
	tX[.,.] = X - means[F.levels, .]

	N_g = F.num_levels
	st_numscalar("r(N_g)",N_g)

}

// End Mata section for s_fe
end

// reset matastrict to what it was prior to loading this file
qui set matastrict `mstrict'

// END CONDITIONAL COMPILATION SECTION

******************** END ALL PROGRAM CODE *******************
exit



*************************** NOTES ***************************

// Conditional coding section above works just once, for a single
// possibly uninstalled Mata package.
// To implement for multiple possibly uninstalled Mata packages,
// use the following comment-out trick (hat-tip to Sergio Correia).

// FTOOLS
// Set comment local to "*" before trying to complile Mata function
// s_fe if required FTOOLS package is not installed.
// Note this runs only once, on loading, and if ftools has not
// been compiled, the check option will trigger compilation.

// ftools can flip matastrict to on, causing program to fail to load.
local mstrict `c(matastrict)'
cap ftools, check
if _rc {
	// fails check, likely not installed, so do not complile
	// _fe Stata program will use (slower) Stata code
	loc c *
}
else {
	// temporarily set matastrict off
	qui set matastrict off
}

// Compile Mata function s_fe.
`c' version 13
`c' mata:

// FE transformation.
// Uses Sergio Correia's FTOOLS package - faster and does not require the data to be sorted.
`c' void s_fe(		string scalar Xnames,
`c' 				string scalar tXnames,
`c' 				string scalar fe,
`c' 				string scalar touse,
`c' 				string scalar toest)
`c' {
`c' 
`c' 	class Factor scalar F
`c' 	F = factor(fe, touse)
`c' 	F.panelsetup()
`c' 
`c' 	Xtokens=tokens(Xnames)
`c' 	X = st_data( ., Xtokens, touse)
`c' 	tXtokens=tokens(tXnames)
`c' 	st_view(tX, ., tXtokens, touse)
`c' 
`c' 	w = F.sort(st_data(., toest, touse))				//  extract toest variable
`c' 	counts = panelsum(w, F.info)
`c' 	means = editmissing(panelsum(F.sort(X), w, F.info) :/ counts, 0)
`c' 	tX[.,.] = X - means[F.levels, .]
`c' 
`c' 	N_g = F.num_levels
`c' 	st_numscalar("r(N_g)",N_g)
`c' 
`c' }
`c' 
`c' // End Mata section for s_fe
`c' end

// reset matastrict
qui set matastrict `mstrict'

// Reset comment local for next package.
loc c

// ... further possibly uninstalled Mata packages ...

exit

