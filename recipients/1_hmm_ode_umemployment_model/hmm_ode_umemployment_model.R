source("recipients/data_handling.R")
library(rstan)
library(cmdstanr)

model <- cmdstan_model("recipients/1_hmm_ode_umemployment_model/hmm_ode_umemployment_model.stan")

# データの用意
data <- list(
  t = nrow(df),
  K = 5,
  delta_t = 1,
  y = df$recipients_density,
  p = df$umemploy
)
data


fit <- model$sample(
  data = data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 3000,
  iter_sampling = 3000,
  refresh = 100
)

fit$cmdstan_summary()














# 潜在状態を計算する
compute_posterior_states <- function(y, p, delta_t = 1, pi, A, lambda, beta, sigma) {
  T_graph <- length(y) - 1  # y[1] to y[T+1]
  K <- length(pi)
  
  log_alpha <- matrix(NA, nrow = T_graph, ncol = K)
  log_beta <- matrix(NA, nrow = T_graph, ncol = K)
  
  # ----- Forward pass -----
  for (k in 1:K) {
#    mu <- y[1] + delta_t * (b1[k]*kero[1]) * y[1] * (1-y[1])
    mu <- y[1] + delta_t * y[1]*(1-y[1])*(lambda[k] + (1-p[1])*beta[k])
    log_alpha[1, k] <- log(pi[k]) + dnorm(y[2], mean = mu, sd = sigma[k], log = TRUE)
  }
  
  for (t_r in 2:T_graph) {
    for (k in 1:K) {
#      mu <- y[t_r] + delta_t * (b1[k]*kero[t_r]) * y[t_r] * (1-y[t_r])
      mu <- y[t_r] + delta_t * y[t_r]*(1-y[t_r])*(lambda[k] + (1-p[t_r])*beta[k])
      temp <- sapply(1:K, function(j) {
        log_alpha[t_r - 1, j] + log(A[j, k])
      })
      log_alpha[t_r, k] <- log_sum_exp(temp) + dnorm(y[t_r + 1], mean = mu, sd = sigma[k], log = TRUE)
    }
  }
  
  # ----- Backward pass -----
  for (k in 1:K) {
#    mu <- y[T_graph] + delta_t * b[k]
    mu <- y[T_graph] + delta_t * y[T_graph]*(1-y[T_graph])*(lambda[k] + (1-p[T_graph])*beta[k])
    log_beta[T_graph, k] <- dnorm(y[T_graph + 1], mean = mu, sd = sigma[k], log = TRUE)
  }
  
  for (t_r in (T_graph - 1):1) {
    for (k in 1:K) {
      temp <- sapply(1:K, function(j) {
#        mu <- y[t_r + 1] + delta_t * b[j]
        mu <- y[t_r + 1] + delta_t * y[t_r + 1]*(1-y[t_r + 1])*(lambda[k] + (1-p[t_r + 1])*beta[k]);
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

lambda <- fit$draws(variables = c("lambda"), format = "draws_matrix") %>% colMeans() %>% as.vector()
lambda
beta <- fit$draws(variables = c("beta"), format = "draws_matrix") %>% colMeans() %>% as.vector()
beta

sigma <- fit$draws(variables = c("sigma"), format = "draws_matrix") %>% colMeans() %>% as.vector()
sigma


y <- df$recipients_density
p <- df$umemploy
gamma <- compute_posterior_states(y, p, delta_t = 1, pi, A, lambda, beta, sigma)
gamma


# 潜在状態の決定
s <- c()
for (i in 1:nrow(gamma)){
  s <- c(s, gamma[i,] %>% which.max)
}
s
s %>% table








# 予測値を入れるリスト
y_pred <- c(y[1])
y_pred_mu <- c(y[1])
y_lower <- c(NA)
y_upper <- c(NA)
# 予測値を計算
# 基礎方程式
de <- function(y, p, lambda, beta, delta_t = 1) {
  return(y + y*(1-y)*delta_t*(lambda+(1-p)*beta))
}
# 平均ではなく中央値とか？
for (i in 1:length(s)){
  y_pred <- c(y_pred, de(y_pred[i], p[i], lambda[s[i]], beta[s[i]]))
  y_lower <- c(y_lower, y_pred[i+1] - 2*sigma[s[i]])
  y_upper <- c(y_upper, y_pred[i+1] + 2*sigma[s[i]])
}

plot(df$year, y_pred)



# 例: gamma_matrix が T行4列の数値行列
gamma_df <- as.data.frame(gamma)
colnames(gamma_df) <- c("state0", "state1", "state2")
gamma_df$year <- 1960:2011

# long形式に変換
gamma_long <- gamma_df %>%
  pivot_longer(cols = starts_with("state"), names_to = "state", values_to = "probability") %>%
  mutate(state = factor(state, levels = c("state0", "state1", "state2")))

ggplot() +
  geom_area(
    data = gamma_long, 
    aes(x = year, y = probability, fill = state), 
    alpha = .5
  ) +
  labs(y = "Probability / Penetration rate of fan heater", x = "Year") +
  scale_fill_manual(values = c("royalblue3", "yellow", "turquoise"), 
                    labels = c("0", "1", "2")
  ) +
  geom_ribbon(
    aes(x = df$year, ymin = c(y_lower), ymax = c(y_upper)), 
    alpha = .25
  ) +
  geom_line(
    data = df, color = "black", lwd = .5,
    aes(year, recipients_density*50, color = c(s, NA))
  ) +
  geom_point(
    aes(df$year, y_pred),
    shape = 4
  ) +
  geom_point(
    data = df, 
    aes(year, recipients_density)
  ) +
  theme_minimal() +
  theme(legend.position = "none")


ggplot() +
  geom_point(
    aes(year, recipients_density), color = "royalblue",
    data = df
  ) + 
  geom_point(
    aes(df$year, y_pred, color = c(s, NA)),
  ) +
  theme_minimal() +
  scale_color_gradientn(colours = c("tomato", "cyan", "black"))





lambda; beta
s


