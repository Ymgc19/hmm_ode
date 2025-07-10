library(cmdstanr)
library(tidyverse)
library(bayesplot)
library(patchwork)

# ========== データを用意 ========== #
enter_univ <- read_csv("univ_entrance/data/enter_univ.csv")
salary <- read_csv("univ_entrance/data/salary.csv")
jobs_to_applicants_ratio <- read_csv("univ_entrance/data/jobs_to_applicants_ratio.csv")

# データを結合
df <- inner_join(enter_univ, salary, by = "year") %>% 
  inner_join(., jobs_to_applicants_ratio, by = "year") %>% 
  dplyr::select(year, enter_univ_rate, sarary_univirsity, salary_high_scool,
                unemployment_rate, jobs_to_applicants_ratio) %>% 
  # 変数を作成
  mutate(
    salary_diff = (sarary_univirsity - salary_high_scool) * .1, # 大卒初任給と高卒初任給の差
    enter_univ_rate = enter_univ_rate / 100, # 大学進学率
    unemployment_rate = unemployment_rate # 完全失業率
  ) %>% 
  arrange(year) # 西暦で並び替え
df %>% glimpse()
df %>% summary

# 大学進学率の推移を描画
df$enter_univ_rate %>% plot


