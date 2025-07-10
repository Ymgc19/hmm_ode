# これは，隠れマルコフ微分方程式をを用いて，大学進学率を扱うファイル
# 説明変数を使っています

# データの読み込みと前処理
source("univ_entrance/data_preprocessing.R")

# ===== MCMCのためのデータ ===== #
data = list(
  t = nrow(df),
  K = 2,
  delta_t = 1,
  y = df$enter_univ_rate,
  salary_diff = df$salary_diff,
  jobs_to_applicants_ratio = df$jobs_to_applicants_ratio,
  unemployment_rate = df$unemployment_rate
)


# ========== モデルを読み込む ========== #
model <- cmdstan_model("univ_entrance/1_hmm_ode/hmm_ode.stan")







# ========== サンプリング ========== #
# サンプリングが不安定なので注意
fit <- model$sample(
  data = data,
  chains = 2, # 4とかにした方がいい
  parallel_chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  refresh = 1000
)
# MCMCの診断は悪い

# 結果を見る
options(scipen = 4)
fit$cmdstan_summary()



# お絵描き
#draws <- fit$draws(
#  variables = c("b0", "b1", "b2", "b3", "sigma"),
#  format = "draws_array"
#)

# traceplot
#color_scheme_set("mix-blue-red")  # 色設定（任意）
#mcmc_trace(draws)

# posterior
#mcmc_dens_overlay(draws)











# ========== 結果を詳しく描画したりしていく ========== #
compute_posterior_states <- function(y, x1, x2, x3, delta_t, pi, A, b0, b1, b2, b3, sigma) {
  T_graph <- length(y) - 1 
  K <- length(pi) 
  log_alpha <- matrix(NA, nrow = T_graph, ncol = K) 
  log_beta <- matrix(NA, nrow = T_graph, ncol = K) 
  # forward
  for (k in 1:K) {
    u_u <- b0[k] + b1[k]*x1[1] + b2[k]*x2[1] + b3[k]*x3[1]
    u_mean <- y[1]*u_u
    mu = y[1] + delta_t * y[1]*(u_u - u_mean)
    log_alpha[1, k] <- log(pi[k]) + dnorm(y[2], mean = mu, sd = sigma[k], log = TRUE)
  }
  for (t_r in 2:T_graph) { 
    for (k in 1:K) { 
      u_u <- b0[k] + b1[k]*x1[t_r] + b2[k]*x2[t_r] + b3[k]*x3[t_r]
      u_mean <- y[t_r]*u_u
      mu <- y[t_r] + delta_t * y[t_r]*(u_u - u_mean)
      
      temp <- sapply(1:K, function(j) { 
        log_alpha[t_r - 1, j] + log(A[j, k]) 
      })
      
      log_alpha[t_r, k] <- log_sum_exp(temp) + dnorm(y[t_r + 1], mu, sigma[k], log = TRUE)
    }
  }
  # backward
  for (k in 1:K) {
    u_u <- b0[k] + b1[k]*x1[T_graph] + b2[k]*x2[T_graph] + b3[k]*x3[T_graph]
    u_mean <- y[T_graph]*u_u
    mu_final_beta <- y[T_graph] + delta_t * y[T_graph] * (u_u - u_mean)
    log_beta[T_graph, k] <- dnorm(y[T_graph + 1], mu_final_beta, sigma[k], log = TRUE)
  }
  for (t_r in (T_graph - 1):1) { 
    for (k in 1:K) { 
      temp <- sapply(1:K, function(j) { 
        u_u <- b0[k] + b1[k]*x1[t_r] + b2[k]*x2[t_r] + b3[k]*x3[t_r]
        u_mean <- y[t_r]*u_u
        
        mu_inner <- y[t_r] + delta_t * y[t_r]*(u_u - u_mean)
        
        log(A[k, j]) + dnorm(y[t_r + 1], mu_inner, sigma[j], log = TRUE) + log_beta[t_r + 1, j]
      })
      log_beta[t_r, k] <- log_sum_exp(temp)
    }
  }
  # posterior copmputatation
  gamma <- matrix(NA, nrow = T_graph, ncol = K)
  for (t_r in 1:T_graph) {
    log_gamma_t <- log_alpha[t_r, ] + log_beta[t_r, ]
    log_gamma_t <- log_gamma_t - log_sum_exp(log_gamma_t) # 正規化
    gamma[t_r, ] <- exp(log_gamma_t)
  }
  return(gamma)
}

log_sum_exp <- function(x) {
  m <- max(x)
  return(m + log(sum(exp(x - m))))
}


# 潜在変数の事後分布を計算
# 変数たち
y <- df$enter_univ_rate
x1 <- df$salary_diff
x2 <- df$jobs_to_applicants_ratio
x3 <- df$unemployment_rate

# パラメタたち
# pi
pi_draws <- fit$draws(variables = "pi", format = "draws_matrix")
pi <- colMeans(pi_draws) %>% as.vector()
# A
library(posterior)
A_mat <- as_draws_matrix(fit$draws("A"))
A_names <- colnames(A_mat)
K <- sqrt(length(grep("^A\\[", A_names)))
A <- matrix(colMeans(A_mat), nrow = K, byrow = TRUE)
# b
b0 <- fit$draws(variables = c("b0"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b1 <- fit$draws(variables = c("b1"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b2 <- fit$draws(variables = c("b2"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b3 <- fit$draws(variables = c("b3"), format = "draws_matrix") %>% colMeans() %>% as.vector()
# sigma
sigma <- fit$draws(variables = c("sigma"), format = "draws_matrix") %>% colMeans() %>% as.vector()

# 潜在変数
gamma <- compute_posterior_states(
  y = y,
  x1 = x1,
  x2 = x2,
  x3 = x3,
  delta_t = 1,
  pi = pi,
  A = A,
  b0 = b0,
  b1 = b1,
  b2 = b2,
  b3 = b3,
  sigma = sigma
)
gamma





# ========== 予測値の計算 ========== #
# 予測値に関してはそんなに重要視しなくてもいい
de <- function(y, x1, x2, x3, b0, b1, b2, b3, delta_t = 1){
  u <- b0 + x1*b1 + x2*b2 + x3*b3
  u_mean <- y*u
  return(
    y + delta_t*y*(u - u_mean)
  )
}

# 潜在状態の決定
s <- c()
for (i in 1:nrow(gamma)){
  s <- c(s, gamma[i,] %>% which.max)
}
s
s %>% table # 潜在変数のtable


# 予測値を入れるリスト
y_pred <- c(df$enter_univ_rate[1])
y_low <- c()
y_high <- c()
# 予測値を計算
# 平均ではなく中央値とか？
for (i in 1:length(s)){
  y_pred <- c(y_pred, de(y_pred[i], x1[i], x2[i], x3[i], b0[s[i]], b1[s[i]], b2[s[i]], b3[s[i]]))
  y_low <- c(y_low, y_pred[i+1] - 2*sigma[s[i]])
  y_high <- c(y_high, y_pred[i+1] + 2*sigma[s[i]])
}

ggplot() +
  # 実測値
  geom_line(
    data = df,
    aes(x = year, y = enter_univ_rate),
    color = "black"
  ) +
  geom_point(
    data = df,
    aes(x = year, y = enter_univ_rate),
    color = "black", shape = 15
  ) +
  # 微分方程式
  geom_ribbon(
    aes(x = df$year[2:nrow(df)],
        ymin = y_low, ymax = y_high),
    #    fill = "skyblue1",
    fill = "black", alpha = .125
  ) +
  geom_line(
    aes(
      1976:(1976 + length(y_pred) - 1), y_pred,
      color = c(s, NA)
    )
  ) +
  geom_point(
    aes(
      1976:(1976 + length(y_pred) - 1), y_pred,
      color = c(s, NA)
    ), shape = 3
  ) +
  labs(x = "year", y = "university enrollment rate") +
  scale_colour_gradientn(colours = c("orange", "turquoise", "purple")) +
  theme_minimal()



