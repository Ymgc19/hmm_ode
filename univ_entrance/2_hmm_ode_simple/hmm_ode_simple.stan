data {
  int t;                                  // 時系列の長さ
  int<lower=1> K;                         // 状態数
  real delta_t;                           // Δt 時間幅
  vector[t] y;                            // 大学進学率
}

parameters {
  simplex[K] pi;                      // 初期状態分布
  array[K] simplex[K] A;              // 遷移確率行列
  ordered[K] H;                       // 相対利得
  vector<lower=1e-12, upper = .1>[K] sigma;       // 各状態の標準偏差
}


model {
  H ~ normal(0, .5);
  sigma ~ gamma(.5, 1);
  
  // HMMは以下
  matrix[t-1, K] log_alpha;
  // 初期ステップ（t=2を予測）
  for (k in 1:K) {
    real mu = y[1] + delta_t * (y[1]*(H[k] - y[1]*H[k]));     // 最初の更新式を定義 t=1からt=2への変化
    log_alpha[1, k] = log(pi[k] + 1e-12) + normal_lpdf(y[2] | mu, sigma[k]);
  }

  // 再帰ステップ（t ≥ 3）
  for (time in 3:t) {
    for (k in 1:K) {
      vector[K] temp;
      // 次時点の平均値を計算
      real mu = y[time-1] + delta_t * (y[time-1]*(H[k] - y[time-1]*H[k]));
      // t時点に状態jにいる確率を計算（対数）
      for (j in 1:K) {
        temp[j] = log_alpha[time - 2, j] + log(A[j, k] + 1e-12);
      }
      log_alpha[time-1, k] = log_sum_exp(temp) + normal_lpdf(y[time] | mu, sigma[k]);
    }
  }

  target += log_sum_exp(log_alpha[t-1]);
}
