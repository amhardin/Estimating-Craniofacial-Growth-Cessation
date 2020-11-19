#' Parse path
#'
#' @param path 
#'
#' @return
#' @export
#'
#' @examples
parse_path <- function(path) {
  pp <- str_split(path, "_", simplify = TRUE)
  tr <- paste(pp[2], pp[3], sep = "_")
  sex <- if_else(pp[1] == "Female", "F", "M")
  trait <- case_when(
    tr == "sella_gonion" ~ "s_g",
    tr == "sella_nasion" ~ "s_n",
    tr == "sella_basion" ~ "s_b",
    tr == "nasion_ans" ~ "n_ans",
    tr == "menton_ans" ~ "ans_m",
    tr == "ans_pns" ~ "ans_pns",
    tr == "gonion_pogonion" ~ "g_po",
    tr == "nasion_basion" ~ "n_b",
    tr == "nasion_menton" ~ "n_m"
  )
  return(list(trait = trait, sex = sex))
}


#' Fix capitalization for datasets
#'
#' @param s (string): Stored name for growth study set
#'
#' @return
#' @export
#'
#' @examples
#cap_Type <- function(s) {
#  s <- case_when(
#    s == "Hyper" ~ "Hyper-divergent",
#    s == "Normo" ~ "Normo-divergent",
#    s == "Hypo" ~ "Hypo-divergent"
#  )
#  return(s)
#}


#' Load data for CF growth model
#'
#' @param variable (string): variable to select, will be renamed to y
#' @param which_sex (string): for sex ("F", "M")
#' @param lower_age (numeric): lower age of range to use for model fitting 
#' @param upper_age (numeric): upper age of range to use for model fitting
#' @param drop_1 
#'
#' @return
#' @export
#'
#' @examples
load_CF_data <- function(variable,
                         which_sex,
                         lower_age = 0,
                         upper_age = 100,
                         drop_1 = FALSE) {
  M <- data.table::fread("C:/Users/annam/Box Sync/Hardin Knigge Shared/Master dataset/master_0805_fixedDENVER_clean.csv") %>% 
    filter(sex == which_sex) %>% 
    filter(age.yrs > lower_age, age.yrs < upper_age) %>% 
    dplyr::select(age.yrs, individual, matches(variable)) %>% 
    rename(age = age.yrs) %>% 
    drop_na()
  
  # Rename variable to y
  names(M)[names(M) == variable] <- "y"
  
  if (drop_1) { # Drop those with 1 obs
    Obs_1 <- M %>% 
      group_by(individual) %>% 
      tally() %>% 
      filter(n == 1L)
    
    M <- M %>% 
      filter(!(individual %in% Obs_1$individual))
  }
  
  # Drop any y < 1 to get rid of 0 values
  M <- M %>% 
    filter(y > 1)
  
  M <- M %>% 
    mutate(individual = factor(individual),
           individual_int = as.integer(individual)) %>% 
    as.data.frame()
  
  dat <- data.frame(age = M$age,
                    y = M$y,
                    ID = M$individual_int)
  
  return(list(M = M, dat = dat))
}


# Calculate first derivative for different models
derivative <- function(P) {
  out <- tibble(age = numeric(),
                `Double Logistic` = numeric())
  
    d_age <- P %>% pull(age)
    d_age <- d_age[-1]
    
    P <- P %>% dplyr::select(-age)
    
    Pd <- apply(P, 2, diff) %>%
      as_tibble()
    Pd <- bind_cols(age = d_age,
                    Pd)
    out <- rbind(out, Pd)

  return(out)
}

nearest_value <- function(vec, target) {
  which(abs(vec - target) == min(abs(vec - target)))
}
