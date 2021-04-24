uDvine_simulation folder contains codes for uDvine
	- AR1000.R: code for using uDvine to approximate AR(p) processes. (can be time consuming)
	- Garch1000.R: code for using uDvine to approximate GARCH(1,1) processes.
	- GarchGJR1000.R: code for using uDvine to approximate GJR-GARCH(1,1) processes.

mDvine_simulation folder contains codes for CuDvine
	- Selection1000.R: code for selecting bivariate copulas for each uDvine of a CuDvine.
	- Estimation1000.R: code for estimating component uDvines + cross-sectional copula for a CuDvine.

Realdata folder contains codes for Australian electricity data
	- Result(general5dim)DCCt.R: code for estimating the CuDvine-DCC model and for out-of-sample forecasting. 
				     (can be time consuming)
	- Result(general5dim)Plots.R: code for estimating the CuDvine-DCC model and for in-sample analysis.


Note that the estimation and model selection procedure for uDvine and CuDvine are computationally fast. However, the out-of-sample forecast can be computationlly slow, as it is a bootstrap-based procedure and each bootstrap run requires finding the zero-point of a one-dimensional conditional CDF function based on the uDvine.