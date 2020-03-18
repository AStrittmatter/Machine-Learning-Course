*! pdslasso 1.0.01 30jan2018
*! authors aa/cbh/ms
*  wrapper for ivlasso

program define pdslasso, eclass sortpreserve
	syntax [anything] [if] [in] ,		///
		[								///
		OLSOPTions(string)				/// options passed to IV or OLS estimation
		* ]

	version 13
	ivlasso `anything' `if' `in', `options' cmdname(pdslasso) ivoptions(`olsoptions')
	
	ereturn local cmd 		pdslasso
	
end

