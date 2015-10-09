
/**
 * FUNCTION: ARSLOGISTICSPLINE
 *
 * Sample beta from the posterior of (beta | y):
 *   y[i] ~ Bernoulli( sigmoid( inner_prod(beta, X[i,]) + offset[i]) )
 *   beta ~ Normal(mean = beta_mean, var = beta_var)
 *
 * Parameters:
 *   offset (nObs x 1): Offset vector
 *   beta_mean (nFactors x 1): Prior mean
 *   beta_var  (nFactors x 1): Prior variance
 *   X (nObs x nFactors): Feature matrix
 *   y (nObs x 1): Observation vector
 *   nFactors (1 x 1): Number of factors (features)
 *   nObs (1 x 1): Number of observations
 *   beta_sample (nFactors x 1):
 *      Input:  Previous sample
 *      Output: New sample
 *   qcent (ncent x 1): A number of centile points for ARS
 *   ncent (1 x 1): Number of centiles
 *   ninit (1 x 1): Number of initial points
 *   x_lower (1 x 1): Lower bound of any sample
 *   x_upper (1 x 1): Upper bound of any sample
 *   x_init (ninit x nFactors):
 *      Input:  Initial points for this sample
 *      Output: Initial points for the next sample
 *   alpha (1 x 1): A parameter used in ARS (default: 0.5)
 *   neval (nFactors x 1): Output: The number of function evaluations performed
 */
void ARSLOGISTICSPLINE(
  const double *offset, double *beta_mean, double *beta_var,
  double *X, double *y, const int *nFactors, const int *nObs,
  double *beta_sample, const double *qcent, const int *ncent,
  const int *ninit, const double *x_lower, const double *x_upper,
  double *x_init, const double *alpha, int *neval);

/**
 * The function to fit logistic regression using ARS sampler
 */
void fitLogistic(const double *offset, double *beta_mean, double *beta_var,
  double *X, double *y, const int *nFactors, const int *nObs, double* beta_sample,
  const double *qcent, const int *ncent, const int *ninit, const double *x_lower,
  const double *x_upper, const int* nSamples);
