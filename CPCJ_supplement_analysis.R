library(zeallot)
library(tidyverse)
library(rstanarm)
library(cowplot)
library(rethinking)
library(lme4)
library(optimx)
library(bayesplot)
library(knitr)
library(kableExtra)
library(parallel)
library(latex2exp)
library(data.table)
library(matrixStats)
library(ggrepel)
library(gganimate)
library(paletteer)

#This directory contains all the models that will be visualized with this code
setwd("C:/Users/annam/OneDrive/Current Projects/GC Methods")

source("CF_functions_GC.R")

filenames <- list.files(path = "C:/Users/annam/OneDrive/Current Projects/GC Methods",
                        pattern = "population.Rds",
                        full.names = TRUE,
                        recursive = TRUE)

args <- c("Female_ans_pns",
          "Female_gonion_pogonion",
          "Female_nasion_basion",
          "Female_nasion_menton",
          "Female_sella_gonion",
          "Male_ans_pns",
          "Male_gonion_pogonion",
          "Male_nasion_basion",
          "Male_nasion_menton",
          "Male_sella_gonion")

names <- c("Female ANS to PNS",
           "Female Gonion to Pogonion",
           "Female Nasion to Basion",
           "Female Nasion to Menton",
           "Female Sella to Gonion",
           "Male ANS to PNS",
           "Male Gonion to Pogonion",
           "Male Nasion to Basion",
           "Male Nasion to Menton",
           "Male Sella to Gonion")

traits <- c("ANS-PNS",
            "Go-Pog",
            "N-Ba",
            "N-Me",
            "S-Go",
            "ANS-PNS",
            "Go-Pog",
            "N-Ba",
            "N-Me",
            "S-Go")

err_val <- c(0.77,
             0.79,
             0.82,
             1.05,
             0.86,
             0.86,
             0.86,
             0.89,
             1.19,
             0.96)

iter <- 10000L
lower_age <- 0L
upper_age <- 100L
pred_lo <- 1L
pred_hi <- 50L
drop_1 <- TRUE
ncores <- 4
pred_len <- 200
max_lo <- 10
max_hi <- 16
pct_of_max <- 0.98 # sets the percentage of the f asymptote that we use for the growth cessation milestone
DL_pars <- c("f", "a1", "b1", "c1", "b2", "c2")
GC1 <- 0.1

summ_all <- NULL
max_all <- NULL
cess1_all <- NULL
cess05_all <- NULL
asym_all <- NULL
abs_all <- NULL
pred_all <- NULL
err_all <- NULL

#looping through 5 traits in 2 sexes
for (i in 1:10) {
  
  path <- args[i]
  c(trait, sex) %<-% parse_path(path)
  print(trait)
  print(sex)
  
  model <- readRDS(filenames[i])
  
  DL <- model$DL_longitudinal
  
  c(M, dat) %<-% load_CF_data(variable = trait,
                              which_sex = sex,
                              lower_age = lower_age,
                              upper_age = upper_age,
                              drop_1 = drop_1)
  
  summ <- M %>%
    group_by(individual) %>%
    mutate(count = n_distinct(age)) %>%
    ungroup() %>%
    summarise(age_min = min(age),
              age_max = max(age),
              n_ID = n_distinct(individual),
              n_image = n_distinct(y),
              med_count = median(count)) %>%
    mutate(trait = trait,
           sex = sex)
  
  summ_all <- rbind(summ_all, summ)
  
  post_DL <- extract.samples(DL)
  a_ID_sim <- matrix(0, nrow = nrow(post_DL$a_ID), 
                     ncol = ncol(post_DL$a_ID))
  
  pred_new <- crossing(age = seq(pred_lo, pred_hi, 
                                 length.out = pred_len),
                       ID = 3L) %>%
    mutate(
      mean = colMeans(link(DL, data = ., refresh = 0,
                           replace = list(a_ID = a_ID_sim))),
      min1 = apply(link(DL, data = ., refresh = 0,
                        replace = list(a_ID = a_ID_sim)), 
                   2, quantile, probs = c(0.01)),
      max99 = apply(link(DL, data = ., refresh = 0,
                         replace = list(a_ID = a_ID_sim)), 
                    2, quantile, probs = c(0.99)) ) %>% 
    dplyr::select(-ID)
  
  rate <- derivative(pred_new) %>%
    rename(rate = "mean") %>%
    select(age, rate)
  
  pred_full <- merge(pred_new, rate, 
                     by = c("age"), all.x = TRUE) %>%
    rename(y = mean) %>%
    mutate(trait = trait,
           sex = sex)
    
    fvalue <- as.numeric(coef(DL)[1])

#max is the peak growth velocity (PGV)
    max <- pred_full %>%
      filter(age > max_lo, age < max_hi) %>%
      mutate(rank = rank(desc(rate)),
             trait = trait,
             sex = sex,
             perc_f = (y / fvalue)) %>%
      filter(rank == 1) %>%
      rename(max_age = "age",
             max_y = "y",
             max_min1 = "min1",
             max_max99 = "max99",
             max_rate = "rate",
             max_percf = "perc_f") %>%
      select(-rank)

#cess1 is GC10% (growth cessation occurs at 10% of the PGV)
    cess1 <- pred_full %>%
      filter(rate < (GC1*max$max_rate)) %>%
      mutate(rank = rank(age),
             trait = trait,
             sex = sex,
             perc_f = (y / fvalue),
             f = fvalue) %>%
      filter(rank == 1) %>%
      rename(cess1_age = "age",
             cess1_y = "y",
             cess1_min1 = "min1",
             cess1_max99 = "max99",
             cess1_rate = "rate",
             cess1_percf = "perc_f") %>%
      select(-rank)

#asym is GCasym (growth cessation occurs at 98% of the f asymptote)
    asym <- pred_full %>%
      mutate(asym = (y / fvalue)) %>%
      filter(asym >= pct_of_max) %>%
      mutate(rank = rank(age),
             trait = trait,
             sex = sex) %>%
      filter(rank == 1) %>%
      rename(asym_age = "age",
             asym_y = "y",
             asym_min1 = "min1",
             asym_max99 = "max99",
             asym_rate = "rate",
             asym_percf = "asym") %>%
      select(-rank)

#abs is GCabs (growth cessation occurs when growth velocity is 0.1 mm/year)
    abs <- pred_full %>%
      filter(rate <= 0.1) %>%
      mutate(rank = rank(age),
             trait = trait,
             sex = sex,
             perc_f = (y /fvalue)) %>%
      filter(rank == 1) %>%
      rename(abs_age = "age",
             abs_y = "y",
             abs_min1 = "min1",
             abs_max99 = "max99",
             abs_rate = "rate",
             abs_percf = "perc_f") %>%
      select(-rank)

#err is GCerr (growth cessation occurs when trait value is within measurement error of f asymptote)
    err <- pred_full %>%
      filter(y >= fvalue - err_val[i]) %>%
      mutate(rank = rank(age),
             trait = trait,
             sex = sex,
             perc_f = (y/fvalue)) %>%
      filter(rank == 1) %>%
      rename(err_age = "age",
             err_y = "y",
             err_min1 = "min1",
             err_max99 = "max99",
             err_rate = "rate",
             err_percf = "perc_f") %>%
      select(-rank)

#combine results from previous loops together    
    max_all <- rbind(max_all, max)
    cess1_all <- rbind(cess1_all, cess1)
    asym_all <- rbind(asym_all, asym)
    abs_all <- rbind(abs_all, abs)
    pred_all <- rbind(pred_all, pred_full)
    err_all <- rbind(err_all, err)
    
}

#reformatting results
est1 <- merge(max_all,
              cess1_all,
              by = c("trait", "sex"))

est2 <- merge(est1,
              asym_all,
              by = c("trait", "sex"))

est3 <- merge(abs_all, err_all,
              by = c("trait", "sex"))

ests <- merge(est2, est3,
              by = c("trait", "sex"))

summ <- ests %>%
  group_by(sex) %>%
  select(max_age, max_y,
         cess1_age, cess1_y,
         asym_age, asym_y,
         abs_age, abs_y,
         err_age, err_y) %>%
  summarise_if(is.numeric, 
               list(~min(.),
                    ~max(.)))

ests2 <- ests %>%
  pivot_longer(cols = ends_with("age"),
               names_to = "age_var",
               values_to = "age_val")

ests3 <- ests2 %>%
  pivot_longer(cols = ends_with("y"),
               names_to = "y_var",
               values_to = "y_val") %>%
  select(trait, sex, age_var, age_val, y_var, y_val) %>%
  mutate(age_name = str_remove(age_var, "_age"),
         y_name = str_remove(y_var, "_y")) %>%
  filter(age_name == y_name)

trait.labs <- c("ANS-PNS", "Go-Po", "N-B",
                "N-M", "S-Go")
names(trait.labs) <- c("ans_pns","g_po", "n_b",
                       "n_m", "s_g")

sex.labs <- c("Female", "Male")
names(sex.labs) <- c("F", "M")

#code to make Figure 1
ggplot() +
  geom_line(data = (pred_all %>%
                      filter(age > 4, age < 30)),
            aes(x = age, y = y, color = sex),
            size = 2) +
  labs(x = "Age (years)",
       y = "Size (mm)",
       shape = "Method") +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  scale_x_continuous(limits = c(5, 30)) +
  facet_grid(rows = vars(trait), 
             cols = vars(sex), 
             scales = "free_y",
             labeller = labeller(trait = trait.labs,
                                 sex = sex.labs)) +
  geom_point(data = (ests3 %>%
               filter(age_name != "cess05")),
             aes(x = age_val,
                 y = y_val,
                 shape = age_name),
             size = 4, stroke = 1) +
  scale_shape_manual(values = c(1:2, 5:6, 3),
                     labels = c(expression(GC[abs]), expression(GC[asym]), 
                                expression(GC["10%"]),
                                expression(GC[err]), "PGV")) +
  guides(color = FALSE) +
  scale_color_manual(values = c("#7785AC", "#84DCC6"))
