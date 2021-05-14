clear

// Data is from Acemoglu, Robinson, and Johnson (2001) "The Colonial Origins of Comparative Development: An Empirical Investigation"
use https://statalasso.github.io/dta/AJR.dta

// We estimate the effect of institutions (avexpr) on income (logpgp95)
// logpgp95 - log of GDP per capita in 1995
// avexpr - average protection against exprorition risk, 1985-1995

* Unconditional OLS estimate
reg logpgp95 avexpr, robust

* Conditional OLS estimate
// We have 24 control variables (latitude, temperature, humidity, ethnical diversity, soil, commodities, etc.)
// The data contains only have 64 country-level observations
reg logpgp95 avexpr lat_abst edes1975 avelf temp* humid* steplow-oilres, robust

* Post-Lasso Double Selection Procedure
// Let the data decide which control variables are important
pdslasso logpgp95 avexpr (lat_abst edes1975 avelf temp* humid* steplow-oilres), robust nois


* Useful Links:
// https://statalasso.github.io/
// https://github.com/aahrens1
// https://economics.mit.edu/files/4123
