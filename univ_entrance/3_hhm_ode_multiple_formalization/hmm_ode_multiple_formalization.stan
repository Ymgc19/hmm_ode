data {
  int t;
  int<lower=1> K;  // 状態数（4を想定）
  real delta_t;
  vector[t] y;
  vector[t] salary_diff;
  vector[t] jobs_to_applicants_ratio;
  vector[t] unemployment_rate;
}

parameters {
  simplex[K] pi;
  array[K] simplex[K] A;

  ordered[3] i;
  vector[3] b;

  real<lower=0> lambda;

  vector<lower=1e-12, upper=0.1>[K] sigma;
}

model {
  matrix[t-1, K] log_alpha;

  // 事前分布
  i ~ normal(0, .5);

  b ~ normal(0, .5);

  lambda ~ normal(0, 1);
  sigma ~ gamma(0.5, 1);

  // 初期ステップ（t=2）
  for (k in 1:K) {
    real mu;
    if (k == 1) {
      mu = y[1] + delta_t * (i[1] + b[1] * salary_diff[1]) * y[1] * (1 - y[1]);
    } else if (k == 2) {
      mu = y[1] + delta_t * (i[2] + b[2] * jobs_to_applicants_ratio[1]) * y[1] * (1 - y[1]);
    } else if (k == 3) {
      mu = y[1] + delta_t * (i[3] + b[3] * unemployment_rate[1]) * y[1] * (1 - y[1]);
    } else {
      mu = y[1] + delta_t * lambda * y[1] * (1 - y[1]);
    }
    log_alpha[1, k] = log(pi[k]) + normal_lpdf(y[2] | mu, sigma[k]);
  }

  // 再帰ステップ（t ≥ 3）
  for (time in 3:t) {
    for (k in 1:K) {
      vector[K] temp;
      real mu;
      if (k == 1) {
        mu = y[time - 1] + delta_t * (i[1] + b[1] * salary_diff[time - 1]) * y[time - 1] * (1 - y[time - 1]);
      } else if (k == 2) {
        mu = y[time - 1] + delta_t * (i[2] + b[2] * jobs_to_applicants_ratio[time - 1]) * y[time - 1] * (1 - y[time - 1]);
      } else if (k == 3) {
        mu = y[time - 1] + delta_t * (i[3] + b[3] * unemployment_rate[time - 1]) * y[time - 1] * (1 - y[time - 1]);
      } else {
        mu = y[time - 1] + delta_t * lambda * y[time - 1] * (1 - y[time - 1]);
      }

      for (j in 1:K) {
        temp[j] = log_alpha[time - 2, j] + log(A[j, k]);
      }
      log_alpha[time - 1, k] = log_sum_exp(temp) + normal_lpdf(y[time] | mu, sigma[k]);
    }
  }

  // 対数尤度を加算
  target += log_sum_exp(log_alpha[t - 1]);
}
