library(dplyr)
library(tidyr)
library(psych)
library(readxl)
library(sjmisc)
library(ggplot2)
library(ggcorrplot)

library(estimatr)
library(lmtest)
library(FNN)
library(Hmisc)

library(plm)
library(modelr)
library(stargazer)

library(foreign)
library(gplots)


prep_week_data <- read.csv("prep_weekly_data.csv")
oblig_data <- read.csv("oblig_discrip_data.csv")
company_data <- read.csv("company_data.csv")

prep_week_data$turnover_norm <- prep_week_data$turnover_norm * 100

# Функция clse позволяет считать корректные стандартные ошибки 
# в случае использования моделей на панельных данных (FE)
# clustered SEs, clustered on "group"... could also cluster on "time" 
# compute Stata-like degrees of freedom adjustment for number of groups
# See http://www.richard-bluhm.com/clustered-ses-in-r-and-stata-2/

clse = function(reg) { 
  # index(reg, "id") returns the id or entity variable vector
  G = length(unique(index(reg, "id")))
  N = length(index(reg, "id"))
  dfa = (G/(G - 1))   # note Bluhm multiplies this by finite-sample df adjustment
  rob = sqrt(diag(dfa*vcovHC(reg, method = "arellano", type = "HC1", cluster = "group")))
  return(rob)
}

# names(prep_week_data)
# [1] "date"                   "bid"                    "ask"                    "bid_ask_bp"             "market_price"          
# [6] "admitted_price"         "close_price"            "indicative"             "turnover"               "num_of_deals"          
# [11] "turnover_in_securities" "current_coupon"         "ytm_ind_p"              "ytm_close_p"            "accrued_interest"      
# [16] "duration"               "mod_duration"           "ytm_oferta_ind_p"       "ytm_oferta_close_p"     "dur_to_p_c_date"       
# [21] "mod_dur_to_p_c_date"    "g_spread_bp"            "exchange"               "indicative_type"        "maturity"              
# [26] "put_call_date"          "bond"                   "cum_count"              "counted_vals"           "company"               
# [31] "sector"                 "year_num"               "quarter_num"            "month_num"              "week_num"              
# [36] "issuance"               "reg_date"               "days_left"              "days_passed"            "turnover_norm"         
# [41] "prop_bid_ask"           "yield2y"                "yield10y"               "moex"                   "moex_ret"              
# [46] "brent"                  "ref_rate"               "usdrub"                 "slope"                  "rating"                
# [51] "rating_num"             "rating_cut"             "liab_current"           "liab_noncurrent"        "default_point"         
# [56] "close"                  "daily_ret"              "trading_volume"         "variance"               "num_of_stocks"         
# [61] "cap"                    "dist_to_default"        "prob_of_default"        "illiquidity"            "year_month"            
# [66] "zero_coupon_bond"       "oferta"                 "sub_maturity"             

# plotmeans(g_spread_bp ~ date, main = "Гетерогенность по времени", data = prep_week_data)
# plotmeans(g_spread_bp ~ bond, main = "Гетерогенность по облигациям", data = prep_week_data)

#################################################### Панельные модели ####################################################

####### illiquidity #######

pool11 <- plm(g_spread_bp ~ illiquidity2 + dist_to_default + days_left + current_coupon + zero_coupon_bond, 
              index = c("bond", "date"), data = prep_week_data, model = "pooling")
pool12 <- plm(g_spread_bp ~ illiquidity2 + prob_of_default + days_left + current_coupon + zero_coupon_bond, 
              index = c("bond", "date"), data = prep_week_data, model = "pooling")

pcdtest(pool11, test = c("lm"))
pbgtest(pool11)

fixed11 <- plm(g_spread_bp ~ illiquidity2 + dist_to_default + days_left + current_coupon + zero_coupon_bond, 
               index = c("bond", "date"), data = prep_week_data, model = "within", effect = "twoways")
fixed12 <- plm(g_spread_bp ~ illiquidity2 + prob_of_default + days_left + current_coupon + zero_coupon_bond, 
               index = c("bond", "date"), data = prep_week_data, model = "within", effect = "twoways")

plmtest(fixed11, effect="twoways", type="ghm")
pcdtest(fixed11, test = c("lm"))
pbgtest(fixed11)

pooltest(pool11, fixed11)
pooltest(pool12, fixed12)

random11 <- plm(g_spread_bp ~ illiquidity2 + dist_to_default + days_left + current_coupon + zero_coupon_bond, 
                index = c("bond", "date"), data = prep_week_data, model = "random", effect = "twoways")
random12 <- plm(g_spread_bp ~ illiquidity2 + prob_of_default + days_left + current_coupon + zero_coupon_bond, 
                index = c("bond", "date"), data = prep_week_data, model = "random", effect = "twoways")

stargazer(random11, random12, type = "text")

phtest(fixed11, random11)
phtest(fixed12, random12)


stargazer(pool11, fixed11, random11, se = list(clse(pool11), clse(fixed11)), type = "html", out = "1_illiquidity_1.html")
stargazer(pool12, fixed12, random12, se = list(clse(pool12), clse(fixed12)), type = "html", out = "1_illiquidity_2.html")

stargazer(pool11, pool12, fixed11, fixed12, random11, random12, digits = 2,
          se = list(clse(pool11), clse(pool12), clse(fixed11), clse(fixed12)), type = "html", out = "1_illiquidity.html"
          )

####### turnover_norm #######

pool21 <- plm(g_spread_bp ~ turnover_norm + dist_to_default + days_left + current_coupon + zero_coupon_bond, 
              index = c("bond", "date"), data = prep_week_data, model = "pooling")
pool22 <- plm(g_spread_bp ~ turnover_norm + prob_of_default + days_left + current_coupon + zero_coupon_bond, 
              index = c("bond", "date"), data = prep_week_data, model = "pooling")

pcdtest(pool21, test = c("lm"))
pbgtest(pool21)

fixed21 <- plm(g_spread_bp ~ turnover_norm + dist_to_default + days_left + current_coupon + zero_coupon_bond, 
               index = c("bond", "date"), data = prep_week_data, model = "within", effect = "twoways")

fixed22 <- plm(g_spread_bp ~ turnover_norm + prob_of_default + days_left + current_coupon + zero_coupon_bond, 
               index = c("bond", "date"), data = prep_week_data, model = "within", effect = "twoways")

plmtest(fixed21, effect="twoways", type="ghm")
pcdtest(fixed21, test = c("lm"))
pbgtest(fixed21)

pooltest(pool21, fixed21)
pooltest(pool22, fixed22)

random21 <- plm(g_spread_bp ~ turnover_norm + dist_to_default + days_left + current_coupon + zero_coupon_bond, 
                index = c("bond", "date"), data = prep_week_data, model = "random", effect = "twoways")

random22 <- plm(g_spread_bp ~ turnover_norm + prob_of_default + days_left + current_coupon + zero_coupon_bond, 
                index = c("bond", "date"), data = prep_week_data, model = "random", effect = "twoways")

phtest(fixed21, random21)
phtest(fixed22, random22)


stargazer(pool21, fixed21, random21, se = list(clse(pool21), clse(fixed21)), type = "html", out = "2_turnover_norm_1.html")
stargazer(pool22, fixed22, random22, se = list(clse(pool22), clse(fixed22)), type = "html", out = "2_turnover_norm_2.html")

stargazer(pool21, pool22, fixed21, fixed22, random21, random22, digits = 2,
          se = list(clse(pool21), clse(pool22), clse(fixed21), clse(fixed22)), type = "html", out = "2_turnover_norm.html"
          )

####### prop_bid_ask ####### 

pool31 <- plm(g_spread_bp ~ prop_bid_ask + dist_to_default + days_left + current_coupon + zero_coupon_bond, 
              index = c("bond", "date"), data = prep_week_data, model = "pooling")
pool32 <- plm(g_spread_bp ~ prop_bid_ask + prob_of_default + days_left + current_coupon + zero_coupon_bond, 
              index = c("bond", "date"), data = prep_week_data, model = "pooling")

pcdtest(pool31, test = c("lm"))
pbgtest(pool31)

fixed31 <- plm(g_spread_bp ~ prop_bid_ask + dist_to_default + days_left + current_coupon + zero_coupon_bond, 
               index = c("bond", "date"), data = prep_week_data, model = "within", effect = "twoways")
fixed32 <- plm(g_spread_bp ~ prop_bid_ask + prob_of_default + days_left + current_coupon + zero_coupon_bond, 
               index = c("bond", "date"), data = prep_week_data, model = "within", effect = "twoways")

plmtest(fixed31, effect="twoways", type="ghm")
pcdtest(fixed31, test = c("lm"))
pbgtest(fixed31)

pooltest(pool31, fixed31)
pooltest(pool32, fixed32)

random31 <- plm(g_spread_bp ~ prop_bid_ask + dist_to_default + days_left + current_coupon + zero_coupon_bond, 
                index = c("bond", "date"), data = prep_week_data, model = "random", effect = "twoways")
random32 <- plm(g_spread_bp ~ prop_bid_ask + prob_of_default + days_left + current_coupon + zero_coupon_bond, 
                index = c("bond", "date"), data = prep_week_data, model = "random", effect = "twoways")

phtest(fixed31, random31)
phtest(fixed32, random32)


stargazer(pool31, fixed31, random31, se = list(clse(pool31), clse(fixed31)), type = "html", out = "3_prop_bid_ask_1.html")

stargazer(pool32, fixed32, random32, se = list(clse(pool32), clse(fixed32)), type = "html", out = "3_prop_bid_ask_2.html")

stargazer(pool31, pool32, fixed31, fixed32, random31, random32, digits = 2,
          se = list(clse(pool31), clse(pool32), clse(fixed31), clse(fixed32)), type = "html", out = "3_prop_bid_ask.html"
          )

####### Итоги #######

stargazer(random11, random21, random31, type = "html", out = "4_final_1.html", digits = 2)
stargazer(random12, random22, random32, type = "html", out = "4_final_2.html", digits = 2)

stargazer(random11, random21, random31, random12, random22, random32, type = "html", out = "4_final.html", digits = 2)

#################################################### Кластеризация #################################################### 

cluster_data <- prep_week_data %>% group_by(bond) %>% 
  summarise(mean = mean(prop_bid_ask, na.rm = TRUE), std = sd(prop_bid_ask, na.rm = TRUE)) %>%
  mutate(norm_mean = (mean - mean(mean)) / sd(mean), norm_std = (std - mean(std)) / sd(std)) %>% data.frame()

overall_mean_mean <- mean(cluster_data$mean)
overall_mean_std <- mean(cluster_data$std)

print(overall_mean_mean)
print(overall_mean_std)

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

rownames(cluster_data) <- cluster_data[, 1]
cluster_data <- cluster_data %>% select(-c(bond, mean, std))

set.seed(123)

fviz_nbclust(cluster_data, kmeans, method = "wss")
fviz_nbclust(cluster_data, kmeans, method = "silhouette")

k2 <- kmeans(cluster_data, centers = 2, nstart = 25)
k3 <- kmeans(cluster_data, centers = 3, nstart = 25)
k4 <- kmeans(cluster_data, centers = 4, nstart = 25)
k5 <- kmeans(cluster_data, centers = 5, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = cluster_data) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = cluster_data) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = cluster_data) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = cluster_data) + ggtitle("k = 5")

library(gridExtra)
grid.arrange(p1, p2, p3, p4, nrow = 2)

cluster_data["cluster"] <- k4$cluster
cluster_data %>% group_by(cluster) %>% summarise(n())

clusters <- cluster_data %>% distinct(cluster) %>% c()

all_models <- list()
i <- 1

for (k in clusters$cluster) {
  print(k)
  
  bonds <- cluster_data %>% filter(cluster == k) %>% rownames()
  
  mini_data <- prep_week_data %>% filter(bond %in% bonds)

  random011 <- plm(g_spread_bp ~ illiquidity2 + dist_to_default + days_left + current_coupon + zero_coupon_bond,
                   index = c("bond", "date"), data = mini_data, model = "random", effect = "twoways")
  all_models[[i]] = random011
  i <- i + 1
  
  random012 <- plm(g_spread_bp ~ illiquidity2 + prob_of_default + days_left + current_coupon + zero_coupon_bond,
                   index = c("bond", "date"), data = mini_data, model = "random", effect = "twoways")
  all_models[[i]] = random012
  i <- i + 1
  
  random021 <- plm(g_spread_bp ~ turnover_norm + dist_to_default + days_left + current_coupon + zero_coupon_bond,
                   index = c("bond", "date"), data = mini_data, model = "random", effect = "twoways")
  all_models[[i]] = random021
  i <- i + 1
  
  random022 <- plm(g_spread_bp ~ turnover_norm + prob_of_default + days_left + current_coupon + zero_coupon_bond,
                   index = c("bond", "date"), data = mini_data, model = "random", effect = "twoways")
  all_models[[i]] = random022
  i <- i + 1
  
  random031 <- plm(g_spread_bp ~ prop_bid_ask + dist_to_default + days_left + current_coupon + zero_coupon_bond,
                   index = c("bond", "date"), data = mini_data, model = "random", effect = "twoways")
  all_models[[i]] = random031
  i <- i + 1
  
  random032 <- plm(g_spread_bp ~ prop_bid_ask + prob_of_default + days_left + current_coupon + zero_coupon_bond,
                  index = c("bond", "date"), data = mini_data, model = "random", effect = "twoways")
  all_models[[i]] = random032
  i <- i + 1
}

stargazer(all_models[[7]], all_models[[8]], all_models[[9]], all_models[[10]], all_models[[11]], all_models[[12]],
          all_models[[13]], all_models[[14]], all_models[[15]], all_models[[16]], all_models[[17]], all_models[[18]],
          all_models[[1]], all_models[[2]], all_models[[3]], all_models[[4]], all_models[[5]], all_models[[6]],
          column.labels = c("1 cluster", "1 cluster", "1 cluster", "1 cluster", "1 cluster", "1 cluster",
                            "2 cluster", "2 cluster", "2 cluster", "2 cluster", "2 cluster", "2 cluster",
                            "3 cluster", "3 cluster", "3 cluster", "3 cluster", "3 cluster", "3 cluster"),
          omit = c("Constant", "days_left", "current_coupon", "zero_coupon_bond"),
          add.lines = list(c("Control", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes",
                                        "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes")),
          type = "html", out = "5_clusters.html", digits = 2)

stargazer(all_models[[1]], all_models[[2]], all_models[[3]], all_models[[4]], all_models[[5]], all_models[[6]],
          title = "3 cluster",
          # omit = c("Constant", "days_left", "current_coupon", "zero_coupon_bond"),
          type = "html", out = "5_cluster_3.html", digits = 2)

stargazer(all_models[[7]], all_models[[8]], all_models[[9]], all_models[[10]], all_models[[11]], all_models[[12]],
          title = "1 cluster",
          # omit = c("Constant", "days_left", "current_coupon", "zero_coupon_bond"),
          type = "html", out = "5_cluster_1.html", digits = 2)

stargazer(all_models[[13]], all_models[[14]], all_models[[15]], all_models[[16]], all_models[[17]], all_models[[18]],
          title = "2 cluster",
          # omit = c("Constant", "days_left", "current_coupon", "zero_coupon_bond"),
          type = "html", out = "5_cluster_2.html", digits = 2)


cluster_another_data <- prep_week_data %>% group_by(bond) %>% 
  summarise(mean1 = mean(prop_bid_ask, na.rm = TRUE), std1 = sd(prop_bid_ask, na.rm = TRUE),
            mean2 = mean(g_spread_bp, na.rm = TRUE), std2 = sd(prop_bid_ask, na.rm = TRUE)) %>% data.frame()

cluster_another_data["cluster"] <- k4$cluster


hull_cyl1 <- cluster_another_data %>%
  group_by(cluster) %>%
  slice(chull(mean1, std1))

p5 <- ggplot(cluster_another_data, aes(x=mean1, y=std1, color=factor(cluster), fill = factor(cluster))) + 
  geom_point() + geom_polygon(data = hull_cyl1, alpha = 0.5) + labs(fill="Кластеры", color="Кластеры") + 
  labs(x = "Среднее значение пропорц. бид-аск спреда", 
       y = "Стандартное отклонение пропорционального бид-аск спреда") + 
  theme(text = element_text(size = 20))

hull_cyl2 <- cluster_another_data %>%
  group_by(cluster) %>%
  slice(chull(mean2, std2))

p6 <- ggplot(cluster_another_data, aes(x=mean2, y=std2, color=factor(cluster), fill = factor(cluster))) + 
  geom_point() + geom_polygon(data = hull_cyl2, alpha = 0.5) + labs(fill="Кластеры", color="Кластеры") + 
  labs(x = "Среднее значение G-спреда, б.п.", 
       y = "Стандартное отклонение G-спреда, б.п.") + 
  theme(text = element_text(size = 20))

grid.arrange(p5, p6, nrow = 1)

