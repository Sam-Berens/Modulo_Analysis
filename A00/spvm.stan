functions {
  // r2k: Converts r (resultant vector) to von Mises kappa
  real r2k(real r) {
    real cr = abs(r);
    real k;
    if (cr < 0.53) {
      k = 2 * cr + cr ^ 3 + 5 * pow(cr, 5) / 6;
    } else if (cr < 0.85) {
      k = -0.4 + 1.39 * cr + 0.43 / (1 - cr);
    } else {
      k = inv(pow(cr, 3) - 4 * pow(cr, 2) + 3 * cr);
    }
    // Clip to avoid numerical overflow.
    if (k > 700) 
      k = 700;
    if (k < -700) 
      k = -700;
    return k * (r >= 0 ? 1 : -1);
  }
  
  // findIdx: Finds arg2 in arg1 and returns the index
  int findIdx(array[] complex v, complex c) {
    real dist = positive_infinity();
    int idx = 0;
    for (i in 1 : size(v)) {
      real d = abs(v[i] - c);
      if (d < dist) {
        dist = d;
        idx = i;
      }
    }
    return idx;
  }
}
data {
  // n: Total number of trails.
  int<lower=0> n;
  
  // x: Predictor variable.
  vector<lower=0>[n] x;
  
  // c: Array of target responses for each trial.
  array[n] int<lower=0, upper=5> c;
  
  // Y: Array of responses for each trial (reals allow NaNs).
  array[n, 6] real Y;
}
transformed data {
  // maxX: Max value of the predictor variable x.
  real<lower=0> maxX = max(x);
  
  // theta: A vector to store all possible angular errors.
  // roots: A corresponding array of complex roots of unity.
  vector[6] theta;
  array[6] complex roots;
  for (i in 1 : 6) {
    real t = fmod((i - 1) * (pi() / 3), 2 * pi());
    theta[i] = t;
    roots[i] = exp(to_complex(0, t));
  }
  
  // Yidx: Array of response indices for theta.
  // nTry: Number of attempts per trial.
  array[n, 6] int Yidx;
  array[n] int<lower=1, upper=6> nTry;
  for (i in 1 : n) {
    int nk = 0;
    for (k in 1 : 6) {
      if (!is_nan(Y[i, k])) {
        nk += 1;
        real yt = fmod((Y[i, k] - c[i]) * (pi() / 3), 2 * pi());
        complex yc = exp(to_complex(0, yt));
        Yidx[i, k] = findIdx(roots, yc);
      } else {
        Yidx[i, k] = 0;
      }
    }
    nTry[i] = nk;
  }
}
parameters {
  // b0: Learning offset parameter.
  real<lower=0, upper=maxX> b0;
  
  // b1: Learning rate parameter.
  real b1;
}
model {
  // Uniform prior for b0
  target += -log(maxX);
  
  // Cauchy prior for b1
  target += cauchy_lpdf(b1 | 0, 0.05);
  
  // lq: A vector of log-likelihoods for all attempts per trial.
  vector[n] lq;
  
  // Compute the trial-wise predicted probabilities.
  vector[1] zero;
  zero[1] = 0;
  for (i in 1 : n) {
    // kappa: Model-predicted von-Mises κ.
    real xPred = log1p_exp(x[i] - b0) * b1;
    real tPred = tanh(xPred);
    real kappa = r2k(tPred);
    
    // pmf: An 6-vector of predicted response probabilities ...
    // ... unsorted, with entries corresponding to theta.
    vector[6] pmf;
    for (k in 1 : 6) {
      pmf[k] = exp(kappa * cos(theta[k]));
    }
    pmf = pmf / sum(pmf);
    
    // rpmf: A vector of predicted response probabilities ...
    // ... sorted, with entries corresponding to Y.
    vector[nTry[i]] rpmf;
    for (k in 1 : nTry[i]) {
      rpmf[k] = pmf[Yidx[i, k]];
    }
    
    // orpmf: rpmf vector offset by a single zero.
    vector[nTry[i]] orpmf = append_row(zero, rpmf[1 : (nTry[i] - 1)]);
    
    // corpmf: Cumulative offset response PMF ...
    // ... (sum probability ruled out before each attempt).
    vector[nTry[i]] corpmf = cumulative_sum(orpmf);
    
    // Update the log-likelihood for the current trial.
    lq[i] = 0;
    for (k in 1 : nTry[i]) {
      lq[i] += log(rpmf[k]) - log(1 - corpmf[k]);
    }
  }
  
  // Add the log-likelihood to the target.
  target += sum(lq);
}