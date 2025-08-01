library(rstan)
library(cmdstanr)

# データの読み込み
source("heater/heater_data_handling.R")

# モデルの読み込み
model <- cmdstan_model("heater/0_heater_simple/heater_simple.stan")

# ======================================================== #
# ===================== K = 4の場合 ======================= #
# ======================================================== #
# データの用意
data <- list(
  t = nrow(df),
  K = 4,
  delta_t = 1,
  y = df$heater_preveil*.01
)
data

# sMCMC回す
fit <- model$sample(
  data = data,
  chains = 2,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  refresh = 1
)


# 結果を見る
fit$cmdstan_summary()


# 潜在状態を計算する
compute_posterior_states <- function(y, delta_t, pi, A, b, sigma) {
  T_graph <- length(y) - 1  # y[1] to y[T+1]
  K <- length(pi)
  
  log_alpha <- matrix(NA, nrow = T_graph, ncol = K)
  log_beta <- matrix(NA, nrow = T_graph, ncol = K)
  
  # ----- Forward pass -----
  for (k in 1:K) {
    mu <- y[1] + delta_t * b[k]
    log_alpha[1, k] <- log(pi[k]) + dnorm(y[2], mean = mu, sd = sigma[k], log = TRUE)
  }
  
  for (t_r in 2:T_graph) {
    for (k in 1:K) {
      mu <- y[t_r] + delta_t * b[k]
      temp <- sapply(1:K, function(j) {
        log_alpha[t_r - 1, j] + log(A[j, k])
      })
      log_alpha[t_r, k] <- log_sum_exp(temp) + dnorm(y[t_r + 1], mean = mu, sd = sigma[k], log = TRUE)
    }
  }
  
  # ----- Backward pass -----
  for (k in 1:K) {
    mu <- y[T_graph] + delta_t * b[k]
    log_beta[T_graph, k] <- dnorm(y[T_graph + 1], mean = mu, sd = sigma[k], log = TRUE)
  }
  
  for (t_r in (T_graph - 1):1) {
    for (k in 1:K) {
      temp <- sapply(1:K, function(j) {
        mu <- y[t_r + 1] + delta_t * b[j]
        log(A[k, j]) + dnorm(y[t_r + 2], mean = mu, sd = sigma[j], log = TRUE) + log_beta[t_r + 1, j]
      })
      log_beta[t_r, k] <- log_sum_exp(temp)
    }
  }
  
  # ----- Posterior state probabilities -----
  gamma <- matrix(NA, nrow = T_graph, ncol = K)
  for (t_r in 1:T_graph) {
    log_gamma_t <- log_alpha[t_r, ] + log_beta[t_r, ]
    log_gamma_t <- log_gamma_t - log_sum_exp(log_gamma_t)
    gamma[t_r, ] <- exp(log_gamma_t)
  }
  
  return(gamma)
}

# 安定な log-sum-exp
log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}





# パラメタを取り出す
pi_draws <- fit$draws(variables = "pi", format = "draws_matrix")
pi <- colMeans(pi_draws) %>% as.vector()
pi

library(posterior)
A_mat <- as_draws_matrix(fit$draws("A"))
A_names <- colnames(A_mat)
K <- sqrt(length(grep("^A\\[", A_names)))
A <- matrix(colMeans(A_mat), nrow = K, byrow = TRUE)
A

b <- fit$draws(variables = c("b"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b

sigma <- fit$draws(variables = c("sigma"), format = "draws_matrix") %>% colMeans() %>% as.vector()
sigma

gamma <- compute_posterior_states(df$heater_preveil*.01, delta_t = 1, pi, A, b, sigma)

# 潜在状態の決定
s <- c()
for (i in 1:nrow(gamma)){
  s <- c(s, gamma[i,] %>% which.max)
}
s
s %>% table






# 予測値を入れるリスト
y <- df$heater_preveil*.01
y_pred <- c(y[1])
y_pred_mu <- c(y[1])
y_low <- c(NA)
y_high <- c(NA)
# 予測値を計算
# 基礎方程式
de <- function(y, b, delta_t = 1) {
  return(y + delta_t*b)
}
# 平均ではなく中央値とか？
for (i in 1:length(s)){
  y_pred <- c(y_pred, de(y[i], b[s[i]]))
  y_low <- c(y_low, y_pred_mu[i+1] - 2*sigma[s[i]])
  y_high <- c(y_high, y_pred_mu[i+1] + 2*sigma[s[i]])
}

plot(df$year, y_pred)

df %>% 
  ggplot() +
  geom_line(
    aes(year, heater_preveil, color = c(s, NA))
  ) +
  geom_point(
    aes(year, heater_preveil, color = c(s, NA))
  ) +
  scale_color_gradientn(colours = c("tomato", "gold", "turquoise")) +
  theme_minimal() +
  theme(legend.position = "none")

# 例: gamma_matrix が T行4列の数値行列
gamma_df <- as.data.frame(gamma)
colnames(gamma_df) <- c("state0", "state1", "state2", "state3")
gamma_df$time <- 1985:2023

# long形式に変換
gamma_long <- gamma_df %>%
  pivot_longer(cols = starts_with("state"), names_to = "state", values_to = "probability") %>%
  mutate(state = factor(state, levels = c("state0", "state1", "state2", "state3")))

ggplot() +
  geom_area(
    data = gamma_long, 
    aes(x = time, y = probability, fill = state), 
    alpha = .5
  ) +
  labs(y = "Probability / Penetration rate of fan heater", x = "Time") +
  scale_fill_manual(values = c("tomato", "gold", "turquoise", "purple"), 
                    labels = c("0", "1", "2", "3")
  ) +
  geom_line(
    data = df, color = "black", lwd = .5,
    aes(year, heater_preveil*.01, color = c(s, NA))
  ) +
  geom_point(
    aes(df$year, y_pred),
    shape = 4
  ) +
  geom_point(
    data = df, 
    aes(year, heater_preveil*.01)
  ) +
  theme_minimal() +
  theme(legend.position = "none")








# ======================================================== #
# ===================== K = 3の場合 ======================= #
# ======================================================== #
# データの用意
data <- list(
  t = nrow(df),
  K = 3,
  delta_t = 1,
  y = df$heater_preveil*.01
)
data

# sMCMC回す
fit <- model$sample(
  data = data,
  chains = 2,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  refresh = 1
)


# 結果を見る
fit$cmdstan_summary()


# 潜在状態を計算する
compute_posterior_states <- function(y, delta_t, pi, A, b, sigma) {
  T_graph <- length(y) - 1  # y[1] to y[T+1]
  K <- length(pi)
  
  log_alpha <- matrix(NA, nrow = T_graph, ncol = K)
  log_beta <- matrix(NA, nrow = T_graph, ncol = K)
  
  # ----- Forward pass -----
  for (k in 1:K) {
    mu <- y[1] + delta_t * b[k]
    log_alpha[1, k] <- log(pi[k]) + dnorm(y[2], mean = mu, sd = sigma[k], log = TRUE)
  }
  
  for (t_r in 2:T_graph) {
    for (k in 1:K) {
      mu <- y[t_r] + delta_t * b[k]
      temp <- sapply(1:K, function(j) {
        log_alpha[t_r - 1, j] + log(A[j, k])
      })
      log_alpha[t_r, k] <- log_sum_exp(temp) + dnorm(y[t_r + 1], mean = mu, sd = sigma[k], log = TRUE)
    }
  }
  
  # ----- Backward pass -----
  for (k in 1:K) {
    mu <- y[T_graph] + delta_t * b[k]
    log_beta[T_graph, k] <- dnorm(y[T_graph + 1], mean = mu, sd = sigma[k], log = TRUE)
  }
  
  for (t_r in (T_graph - 1):1) {
    for (k in 1:K) {
      temp <- sapply(1:K, function(j) {
        mu <- y[t_r + 1] + delta_t * b[j]
        log(A[k, j]) + dnorm(y[t_r + 2], mean = mu, sd = sigma[j], log = TRUE) + log_beta[t_r + 1, j]
      })
      log_beta[t_r, k] <- log_sum_exp(temp)
    }
  }
  
  # ----- Posterior state probabilities -----
  gamma <- matrix(NA, nrow = T_graph, ncol = K)
  for (t_r in 1:T_graph) {
    log_gamma_t <- log_alpha[t_r, ] + log_beta[t_r, ]
    log_gamma_t <- log_gamma_t - log_sum_exp(log_gamma_t)
    gamma[t_r, ] <- exp(log_gamma_t)
  }
  
  return(gamma)
}

# 安定な log-sum-exp
log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# パラメタを取り出す
pi_draws <- fit$draws(variables = "pi", format = "draws_matrix")
pi <- colMeans(pi_draws) %>% as.vector()
pi

library(posterior)
A_mat <- as_draws_matrix(fit$draws("A"))
A_names <- colnames(A_mat)
K <- sqrt(length(grep("^A\\[", A_names)))
A <- matrix(colMeans(A_mat), nrow = K, byrow = TRUE)
A

b <- fit$draws(variables = c("b"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b

sigma <- fit$draws(variables = c("sigma"), format = "draws_matrix") %>% colMeans() %>% as.vector()
sigma

gamma <- compute_posterior_states(df$heater_preveil*.01, delta_t = 1, pi, A, b, sigma)

# 潜在状態の決定
s <- c()
for (i in 1:nrow(gamma)){
  s <- c(s, gamma[i,] %>% which.max)
}
s
s %>% table


# 予測値を入れるリスト
y <- df$heater_preveil*.01
y_pred <- c(y[1])
y_pred_mu <- c(y[1])
y_low <- c(NA)
y_high <- c(NA)
# 予測値を計算
# 基礎方程式
de <- function(y, b, delta_t = 1) {
  return(y + delta_t*b)
}
# 平均ではなく中央値とか？
for (i in 1:length(s)){
  y_pred <- c(y_pred, de(y[i], b[s[i]]))
  y_low <- c(y_low, y_pred_mu[i+1] - 2*sigma[s[i]])
  y_high <- c(y_high, y_pred_mu[i+1] + 2*sigma[s[i]])
}

plot(df$year, y_pred)

df %>% 
  ggplot() +
  geom_line(
    aes(year, heater_preveil, color = c(s, NA))
  ) +
  geom_point(
    aes(year, heater_preveil, color = c(s, NA))
  ) +
  scale_color_gradientn(colours = c("tomato", "gold", "turquoise")) +
  theme_minimal() +
  theme(legend.position = "none")

# 例: gamma_matrix が T行4列の数値行列
gamma_df <- as.data.frame(gamma)
colnames(gamma_df) <- c("state0", "state1", "state2")
gamma_df$time <- 1985:2023

# long形式に変換
gamma_long <- gamma_df %>%
  pivot_longer(cols = starts_with("state"), names_to = "state", values_to = "probability") %>%
  mutate(state = factor(state, levels = c("state0", "state1", "state2")))

ggplot() +
  geom_area(
    data = gamma_long, 
    aes(x = time, y = probability, fill = state), 
    alpha = .5
  ) +
  labs(y = "Probability / Penetration rate of fan heater", x = "Time") +
  scale_fill_manual(values = c("tomato", "gold", "turquoise"), 
                    labels = c("0", "1", "2", "3")
  ) +
  geom_line(
    data = df, color = "black", lwd = .5,
    aes(year, heater_preveil*.01, color = c(s, NA))
  ) +
  geom_point(
    aes(df$year, y_pred),
    shape = 4
  ) +
  geom_point(
    data = df, 
    aes(year, heater_preveil*.01)
  ) +
  theme_minimal() +
  theme(legend.position = "none")

