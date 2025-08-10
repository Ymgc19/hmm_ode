// 現在は生活保護受給の利得を1に固定している


data {
  int t;                                  // 時系列の長さ
  int<lower=1> K;                         // 状態数
  real delta_t;                           // Δt 
  vector[t] y;                            // 大学進学率（目的変数）
  vector[t] p;                            // 完全失業率
  int N;                                  // 仕事に就くのに挑戦する回数
}

parameters {
  simplex[K] pi;                        // 初期状態分布
  array[K] simplex[K] A;                // 遷移確率行列
  vector<lower=0, upper=1>[K] lambda;                    // 働く利得のパラメタ
  vector<lower=0>[K] u_recipients;
  vector<lower=1e-12, upper = .1>[K] sigma;                      // 各状態の標準偏差
}


model {
  matrix[t-1, K] log_alpha;

//  lambda ~ normal(.25, .1);
  u_recipients ~ normal(1, .001);


  // ========== 初期ステップ（t=2） ========== //
  for (k in 1:K) {
    real mu;
    // ========== k = 1 の場合 ========== //
    if (k == 1) {
      real q = 1 - pow(p[1], 1.0 / N );
      real u = 0;
      for (m in 0:N) {
        real bin_prob = exp(binomial_lpmf(m | N, q));
        u += (pow(lambda[k] + 1, m) - 1) * bin_prob;
//        u += lambda[k]*log(m+1) * bin_prob;
      }
      // 次時点の予測
      mu = y[1] + delta_t * y[1] * (1 - y[1]) * (u_recipients[k] - (y[1] * u_recipients[k] + (1 - y[1]) * u));
    } 
    // ========== k == 2 の場合 ========== //
    else if (k == 2) {
      real q = 1 - pow(p[1], 1.0 / N );
      real u = 0;
      for (m in 0:N) {
        real bin_prob = exp(binomial_lpmf(m | N, q));
        u += lambda[k]*log(m+1) * bin_prob;
      }
      // 次時点の予測
      mu = y[1] + delta_t * y[1] * (1 - y[1]) * (u_recipients[k] - (y[1] * u_recipients[k] + (1 - y[1]) * u));
    }  
    // ========== k == 3 の場合 ========== //
    else if (k == 3) {
      real q = 1 - pow(p[1], 1.0 / N );
      real u = 0;
      for (m in 0:N) {
        real bin_prob = exp(binomial_lpmf(m | N, q));
//        u += lambda[k]*log(m+1) * bin_prob;
        u += lambda[k]*m * bin_prob;
      }
      // 次時点の予測
      mu = y[1] + delta_t * y[1] * (1 - y[1]) * (u_recipients[k] - (y[1] * u_recipients[k] + (1 - y[1]) * u));
    } 
    // ========== そのほかの場合 ========== //
    else {
      mu = y[1] + 100 * delta_t * y[1] * (1 - y[1]);
    }
    log_alpha[1, k] = log(pi[k]) + normal_lpdf(y[2] | mu, sigma[k]);
  }

  // ========== 再帰ステップ（t ≥ 3） ========== //
  for (time in 3:t) {
    for (k in 1:K) {
      vector[K] temp;
      real mu;
      if (k == 1) {
        real q = 1 - pow(p[time-1], 1.0 / N );
        real u = 0;
        for (m in 0:N) {
          real bin_prob = exp(binomial_lpmf(m | N, q));
          u += (pow(lambda[k] + 1, m) - 1) * bin_prob;
//          u += lambda[k]*log(m+1) * bin_prob;
        }
        // 次時点の予測
        mu = y[time-1] + delta_t * y[time-1] * (1 - y[time-1]) * (u_recipients[k] - (y[time-1] * u_recipients[k] + (1 - y[time-1]) * u));
      } else if (k == 2) {
        real q = 1 - pow(p[time-1], 1.0 / N );
        real u = 0;
        for (m in 0:N) {
          real bin_prob = exp(binomial_lpmf(m | N, q));
          u += lambda[k]*log(m+1) * bin_prob;
        }
        // 次時点の予測
        mu = y[time-1] + delta_t * y[time-1] * (1 - y[time-1]) * (u_recipients[k] - (y[time-1] * u_recipients[k] + (1 - y[time-1]) * u));
      } else if (k == 3) {
        real q = 1 - pow(p[time-1], 1.0 / N );
        real u = 0;
        for (m in 0:N) {
          real bin_prob = exp(binomial_lpmf(m | N, q));
          u += lambda[k]*m * bin_prob;
//          u += lambda[k]*log(m+1) * bin_prob;
        }
        // 次時点の予測
        mu = y[time-1] + delta_t * y[time-1] * (1 - y[time-1]) * (u_recipients[k] - (y[time-1] * u_recipients[k] + (1 - y[time-1]) * u));
      } else {
        mu = y[time - 1] + delta_t * 100 * y[time - 1] * (1 - y[time - 1]);
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
