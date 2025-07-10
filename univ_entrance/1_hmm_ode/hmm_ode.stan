data {
  int t;                                  // 時系列の長さ
  int<lower=1> K;                         // 状態数
  real delta_t;                           // Δt 
  vector[t] y;                            // 大学進学率（目的変数）
  vector[t] salary_diff;                  // 高卒と大卒の所得差
  vector[t] jobs_to_applicants_ratio;     // 有効求人倍率
  vector[t] unemployment_rate;            // 完全失業率
}

parameters {
  simplex[K] pi;                        // 初期状態分布
  array[K] simplex[K] A;                // 遷移確率行列
  ordered[K] b0;                        // 切片
  vector[K]  b1;  // 説明変数1の係数
  vector[K]  b2;  // 説明変数2の係数
  vector[K]  b3;  // 説明変数3の係数
  vector<lower=1e-12, upper = .1>[K] sigma;         // 各状態の標準偏差
}


model {
  b0 ~ normal(0, .5);
  b1 ~ normal(0, .5);
  b2 ~ normal(0, .5);
  b3 ~ normal(0, .5);
  sigma ~ gamma(.5, 1);
  
  matrix[t-1, K] log_alpha;
  
  // ここら辺もっと簡単にできないものかな
  // 初期ステップ（t=2を予測）
  for (k in 1:K) {
    real mu = y[1] + delta_t * (y[1]*(b0[k] + salary_diff[1]*b1[k] + jobs_to_applicants_ratio[1]*b2[k] + unemployment_rate[1]*b3[k] 
    - y[1]*(b0[k] + salary_diff[1]*b1[k] + jobs_to_applicants_ratio[1]*b2[k] + unemployment_rate[t-1]*b3[k]))); 
    log_alpha[1, k] = log(pi[k]) + normal_lpdf(y[2] | mu, sigma[k]);
  }

  // 再帰ステップ（t ≥ 3）
  for (time in 3:t) {
    for (k in 1:K) {
      vector[K] temp;
      // 次時点の平均値を計算
      real mu = y[time-1] + delta_t * (y[time-1]*(b0[k] + salary_diff[time-1]*b1[k] + jobs_to_applicants_ratio[time-1]*b2[k] + unemployment_rate[time-1]*b3[k] 
      - y[time-1]*(b0[k] + salary_diff[time-1]*b1[k] + jobs_to_applicants_ratio[time-1]*b2[k] + unemployment_rate[time-1]*b3[k])));
      // t時点に状態jにいる確率を計算（対数）
      for (j in 1:K) {
        temp[j] = log_alpha[time - 2, j] + log(A[j, k]);
      }
      log_alpha[time-1, k] = log_sum_exp(temp) + normal_lpdf(y[time] | mu, sigma[k]);
    }
  }

  target += log_sum_exp(log_alpha[t-1]);
}
