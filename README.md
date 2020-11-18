# RStatisticalSoftware
A repository to hold statistical software.
Better documentation for this software is coming.

DistCorScreen:
A package that implements distance correlation based screening. It supports both backwards and forward selection approaches. 
Note that this package uses the partial distance correlation to screen unlike other implementations of distance correlation screening.
The package was designed for screening large datasets of 10000+ variables. To handle these big datasets in reasonable time,
I implemented forward batch selection and clustered backward selection. There are a number of parameters that can be tuned to improve the performance
and computational speed of the program. Note this package is in development and has not been documented yet.

Halgene:
A (meta) machine learning package intended to allow using machine learning on huge datasets (e.g. methylation arrays with 500,000 covariates).
The general idea is to use the predictions of cheap baseline learners (e.g. glmnet) as dimension reductions of the dataset. 
Using these predictions as input, machine learning algorithms (e.g. monotone Highly Adaptive Lasso) can fit the best meta-function of the predictions.
By dividing the covariates into groups and fitting base line learners on each group, one can obtain multiple predictions which then can be combined optimally
into a single prediction via a machine learning algorithm. The package also supports prescreening and PCA. 
Additionally, by making some adjustments to the machine learning pipeline sl3, this package supports data too big to load into RAM which is stored on filebacked Big.Matrix 
objects. The machine learning and cross validation can be done without ever loading the data into RAM by employing BigLasso as baseline learner 
(and possibly screening with a BigMatrix supported screener). The Big.Matrix compatible sl3_Tasks and Big.Matrix compatible learners are completely general 
and can be used for any purpose, not just this package.

(OLD, see tmle3 fork for updated version) LTMLE:
Building upon tlverse/tmle3, I am building a general implementation of likelihood-based LTMLE. This package will support very general longitudinal data structures
(e.g. continuous-time, time dependent covariates, censoring, and counting process/monitoring processes).


OneStepSurvival:
A fast and scalable implementation of the one-step universal submodel TMLE estimator for the treatment-specific survival curves. Given baseline covariates, treatment assignment,
and censoring/failure times, the package obtains a nonparametric efficient estimator for the treatment-specific survival curves. By using the one-step universal submodel,
the entire survival curve (as a function of time) can be targeted simultaneously, giving a compatible fit and simultaneous confidence intervals.

Subsemble:
A general implementation of the Subsemble algorithm using the tlverse/sl3 machine learning pipeline. This implementation supports arbitrary baseline learners and meta learners,
and also accepts stacks/learners indexed by tuning parameters. The meta learner will only combine the predictions corresponding to a single learner in the stack/a single tuning parameter.
As a result, the Subsemble algorithm can be indexed by a number of learners/parameters and the optimal choice can be chosen with cross validation.
The subsemble implementation is also applied to hal9001 to create a very scalable Highly Adaptive Lasso estimator that scales to large datasets.

hal9001fast:
A fork of tlverse/hal9001 with additional screening, formula-based interface, and higher-order Hal implemented.
