library(tidyverse)
library(slider)
library(cmdstanr)

# ====================================================================================================== #
# =========================================== データの読み込み ========================================= #
# ====================================================================================================== #
df2001 <- read_csv("covid_self_restraint/data/covid_self_restraint/Gender_age_01.csv")
df2002 <- read_csv("covid_self_restraint/data/covid_self_restraint/Gender_age_02.csv")
df2003 <- read_csv("covid_self_restraint/data/covid_self_restraint/Gender_age_03.csv")
df2004 <- read_csv("covid_self_restraint/data/covid_self_restraint/Gender_age_04.csv")
df2005 <- read_csv("covid_self_restraint/data/covid_self_restraint/Gender_age_05.csv")
df2006 <- read_csv("covid_self_restraint/data/covid_self_restraint/Gender_age_06.csv")
df2007 <- read_csv("covid_self_restraint/data/covid_self_restraint/Gender_age_07.csv")
df2008 <- read_csv("covid_self_restraint/data/covid_self_restraint/Gender_age_08.csv")
df2009 <- read_csv("covid_self_restraint/data/covid_self_restraint/Gender_age_09.csv")
df2010 <- read_csv("covid_self_restraint/data/covid_self_restraint/Gender_age_10.csv")
df2011 <- read_csv("covid_self_restraint/data/covid_self_restraint/Gender_age_11.csv")
df2012 <- read_csv("covid_self_restraint/data/covid_self_restraint/Gender_age_12.csv")
df2101 <- read_csv("covid_self_restraint/data/covid_self_restraint/Gender_age_2021_01.csv")
df2102 <- read_csv("covid_self_restraint/data/covid_self_restraint/Gender_age_2021_02.csv")
df2103 <- read_csv("covid_self_restraint/data/covid_self_restraint/Gender_age_2021_03.csv")

df <- bind_rows(df2001, df2002, df2003, df2004, df2005, df2006, 
                df2007, df2008, df2009, df2010, df2011, df2012, 
                df2101, df2102, df2103)
df %>% glimpse()

# 新規感染者数
df_infectious <- read_csv("covid_self_restraint/data/covid_newly_infected/newly_confirmed_cases_daily.csv") %>% 
  rename(date = Date) %>% 
  mutate(date = as.character(date))
df_infectious %>% glimpse()

df <- left_join(df, df_infectious, by = "date")
df %>% glimpse()




# ====================================================================================================== #
# =========================================== 任意県のデータ =========================================== #
# ====================================================================================================== #
# ここでは東京
df_pref <- df %>%
  filter(prefecture == "東京都") %>%
  mutate(date = date %>% as_date(format = "%m/%d/%y"),
         # 緊急事態宣言下のダミー変数作成
         emergency = case_when(
           (date >= as.Date("2020-04-07")) ~ 1,
           TRUE ~ 0
         ),
         #緊急事態宣言解除後の経過日数
         pass_days = case_when(
           (date > as.Date("2020-05-25")) & (date < as.Date("2021-01-08")) ~ as.numeric(difftime(date, as.Date("2020-05-25"), units = "days")),
           date > as.Date("2021-03-21") ~ as.numeric(difftime(date, as.Date("2021-03-21"), units = "days")),
           T ~ 0
         )
  )

df_pref %>% glimpse()


# ====================================================================================================== #
# =========================================== 簡単に可視化！ =========================================== #
# ====================================================================================================== #
# 1週間の移動平均を計算
df_pref <- df_pref %>% 
  mutate(
    M30_smooth = slide_dbl(M30, mean, .before = 6, .complete = T),
    F30_smooth = slide_dbl(F30, mean, .before = 6, .complete = T)
  )
df_pref %>% glimpse()

# rawデータ
df_pref %>% 
  ggplot() +
  geom_point(aes(date, M30, color = emergency), shape = 3) +
  geom_point(aes(date, F30, color = emergency), shape = 5) +
  ylim(-.05, 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_gradient(low = "black", high = "tomato")

# 平滑化したデータ
df_pref %>% 
  ggplot() +
  geom_point(aes(date, M30_smooth, color = emergency), shape = 3) +
  geom_point(aes(date, F30_smooth, color = emergency), shape = 5) +
  ylim(-.05, 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_gradient(low = "black", high = "tomato")
# 女性のほうが自粛率が高い
# 女性のほうが接客業が多く，その結果として宣言の影響を強く受けるから？


