library(tidyverse)
library(zeallot)
library(rstanarm)
library(rethinking)
library(rgenoud)
library(parallel)
library(data.table)
library(rlist)
library(readxl)

source("C:/Users/annam/Documents/GitHub/Craniofacial_Growth/IADR code/CF_functions.R")

#########################################################################

refresh <- 0  # Print interval for stan sampler

t0 <- Sys.time()
print(t0)

#args <- commandArgs(TRUE)

args <- c("Female_gonion_pogonion",
          "Female_sella_gonion",
          "Female_ans_pns",
          "Female_nasion_basion",
          "Female_nasion_menton",
          "Male_gonion_pogonion",
          "Male_sella_gonion",
          "Male_ans_pns",
          "Male_nasion_basion",
          "Male_nasion_menton")

path <- args[1]

models <- list()

print(path)

c(trait, sex) %<-% parse_path(path)

print(trait)
print(sex)

###

iter <- 10000L
lower_age <- 0L
upper_age <- 100L
pred_lo <- 6L
pred_hi <- 22L
drop_1 <- TRUE
ncores <- 2

Priors <- read_excel("C:/Users/annam/Documents/GitHub/Craniofacial_Growth/IADR code/DL_Priors.xlsx")

############
# read in data
c(M, dat) %<-% load_CF_data(variable = trait,
                            which_sex = sex,
                            lower_age = lower_age,
                            upper_age = upper_age,
                            drop_1 = drop_1)


## ------------------------------------------------------------------------
# Starting values

start_values_DL <- Priors[Priors$path == path, 2:7] %>% as.numeric()
print(start_values_DL)

constraints_DL <- list(
  b1 = if_else(start_values_DL[3] < 0, "upper=0", "lower=0"),
  c1 = if_else(start_values_DL[4] < 0, "upper=0", "lower=0"),
  c2 = if_else(start_values_DL[6] < 0, "upper=0", "lower=0"))
print(constraints_DL)

##  stan_models #########################################################

flist_s <- paste(
  "y ~ normal(mu, sigma)",
  "mu <- a1 / (1 + exp(-1 * b1 * (age + c1))) +
        (f - a1) / (1 + exp(-1 * b2 * (age + c2))) +
        a_ID[ID]",
  paste("f ~ normal(", start_values_DL[1], ", 0.5)"),
  paste("a1 ~ normal(", start_values_DL[2], ", 0.5)"),
  paste("b1 ~ cauchy(", start_values_DL[3], ", 1)"),
  paste("c1 ~ normal(", start_values_DL[4], ", 1)"),
  paste("b2 ~ normal(", start_values_DL[5], ", 0.2)"),
  paste("c2 ~ normal(", start_values_DL[6], ", 1)"),
  "sigma ~ cauchy(0, 0.5)",
  "a_ID[ID] ~ normal(0, sigma_ID)",
  "sigma_ID ~ cauchy(0, 1)",
  sep = ";"
)
flist <- lapply(strsplit(flist_s, ";")[[1]],
                function(x) parse(text = x)[[1]])
assign("DL_growth_Type",
       map2stan(
         flist,
         data = dat,
         iter = iter,
         control = list(adapt_delta = 0.99,
                        max_treedepth = 15),
         cores = ncores,
         chains = 4,
         refresh = refresh,
         constraints = constraints_DL,
         start = list(f = start_values_DL[1], a1 = start_values_DL[2],
                      b1 = start_values_DL[3], c1 = start_values_DL[4],
                      b2 = start_values_DL[5], c2 = start_values_DL[6])
       )
)
models <- list.append(models, DL_growth_Type = DL_growth_Type)

write_rds(models, path = paste0(path, ".Rds"), compress = "gz")
print(Sys.time() - t0)
