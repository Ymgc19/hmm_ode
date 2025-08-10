library(tidyverse)

df <- read_csv("recipients/recipients.csv") %>% 
  dplyr::select(year, umemploy, recipients, population, ratio) %>% 
  mutate(recipients_density = recipients / population)

plot(df$year, df$recipients_density)
plot(df$year, df$umemploy)
plot(df$umemploy, df$recipients_density)
plot(df$umemploy, df$ratio)

hist(df$ratio)
summary(df$ratio)

df %>% glimpse()




