functions {
  // findIdx: Finds arg2 in arg1 and returns the index
  // This is used to identify which angular error was made on each attempt
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
  // nTrials: Number of trials.
  int<lower=1> nTrials;
  
  // nPairs: Number of distinct pairs.
  int<lower=1> nPairs;
  
  // pairId: Identity of the pair tested on each trial.
  array[nTrials] int<lower=1, upper=nPairs> pairId;
  
  // nTypes: Number of pair types.
  int<lower=1> nTypes;
  
  // id2type: An array of integers encoding the relationship ...
  //          between pairId and pairType
  array[nPairs] int<lower=1, upper=nTypes> id2type;
  
  // x: Predictor variable.
  vector<lower=0>[nTrials] x;
  
  // c: Array of target responses for each trial.
  array[nTrials] int<lower=0, upper=5> c;
  
  // Y: Array of responses for each trial (reals allow NaNs).
  array[nTrials, 6] real Y;
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
  array[nTrials, 6] int Yidx;
  array[nTrials] int<lower=1, upper=6> nTry;
  for (iTrial in 1 : nTrials) {
    int nk = 0;
    for (k in 1 : 6) {
      if (!is_nan(Y[iTrial, k])) {
        nk += 1;
        real yt = fmod((Y[iTrial, k] - c[iTrial]) * (pi() / 3), 2 * pi());
        complex yc = exp(to_complex(0, yt));
        Yidx[iTrial, k] = findIdx(roots, yc);
      } else {
        Yidx[iTrial, k] = 0;
      }
    }
    nTry[iTrial] = nk;
  }
}
parameters {
  // alpha1 and alpha2:
  //    Distributional parameters for learning offsets of the 
  //    same type (mu and sigma for normal distribution).
  vector[nTypes] alpha1;
  vector<lower=0>[nTypes] alpha2;
  
  // beta1 and beta2:
  //    Distributional parameters for learning rates of the 
  //    same type (mu and sigma for normal distribution).
  vector[nTypes] beta1;
  vector<lower=0>[nTypes] beta2;
  
  // zau: Learning offset parameter (one per pair).
  //      (type-centered, unbounded, and unscaled).
  vector[nPairs] zau;
  
  // zb: Learning rate parameter (one per pair).
  //      (type-centered).
  vector[nPairs] zb;
}
transformed parameters {
  // au: Unbounded and unscaled learning offset parameter (one per pair).
  // b: Learning rate parameter (one per pair).
  vector[nPairs] au;
  vector[nPairs] b;
  for (iPair in 1 : nPairs) {
    int typeidx = id2type[iPair];
    au[iPair] = alpha1[typeidx] + (zau[iPair]*alpha2[typeidx]);
    b[iPair] = beta1[typeidx] + (zb[iPair]*beta2[typeidx]);
  }

  // a: Bounded and scaled learning offset parameter (one per pair).
   vector<lower=0, upper=maxX>[nPairs] a = inv_logit(au) * maxX;
}
model {
  // Logistic prior for alpha1.
  // Half-normal prior for alpha2.
  // Normal prior for beta1.
  // Half-normal prior for beta2.
  for (iType in 1 : nTypes) {
    alpha1[iType] ~ logistic(0, 1);
    alpha2[iType] ~ normal(1, 0.05) T[0,];
    beta1[iType] ~ normal(0, 0.1);
    beta2[iType] ~ normal(0.1, 0.005) T[0,];
  }
  
  // Parameter distributions for each pair type:
  // Normal(0,1) for zau;
  // Normal(0,1) for zb;
  // Yielding ...
  // Normal(alpha1,alpha2) for au;
  // Normal(beta1,beta2) for b;
  for (iPair in 1 : nPairs) {
    zau[iPair] ~ normal(0, 1);
    zb[iPair] ~ normal(0, 1);
  }
  
  // lq: A vector of log-likelihoods for all attempts per trial.
  vector[nTrials] lq;
  
  // Compute the trial-wise predicted probabilities.
  for (iTrial in 1 : nTrials) {
    // pid: ID of the pair tested on this trial
    int pid = pairId[iTrial];
    
    // Compute binomial probabilities: prCorr, prInco.
    // xPred ∈ [-∞, ∞], given by SoftPlus(x-a)*b.
    real xPred = log1p_exp(x[iTrial] - a[pid]) * b[pid];
    real prCorr = inv_logit(-log(6) + xPred);
    real prInco = 1 - prCorr;
    // Split prInco among all incorrect responses.
    real prBads = prInco / 5;

    // lpmf: A 6-vector of predicted response log-probabilities ...
    // ... unsorted, with entries corresponding to theta.
     vector[6] lpmf = log([prCorr,prBads,prBads,prBads,prBads,prBads]');

    // rlpmf: A vector of predicted response log-probabilities ...
    // ... sorted, with entries corresponding to Y.
    vector[nTry[iTrial]] rlpmf;
    for (k in 1 : nTry[iTrial]) {
      rlpmf[k] = lpmf[Yidx[iTrial, k]];
    }
    
    // Update the log-likelihood for the current trial.
    real cumsum_lq = negative_infinity();
    lq[iTrial] = 0;
    for (k in 1 : nTry[iTrial]) {
      cumsum_lq = fmin(cumsum_lq, -1e-15); // Protection
      lq[iTrial] += rlpmf[k] - log1m_exp(cumsum_lq);
      cumsum_lq = log_sum_exp(cumsum_lq, rlpmf[k]); 
    }
  }
  
  // Add the log-likelihood to the target.
  target += sum(lq);
}
