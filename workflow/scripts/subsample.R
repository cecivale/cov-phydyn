# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-euphydyn
#        V\ Y /V    Subsample sequences script
#    (\   / - \     29 July 2022
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------


subsample <- function(metadata_file, include_file, exclude_file, n, method,  seq_id, cases_file, seed, output_file) {
  
  set.seed(seed)
  df <- read_tsv(metadata_file) %>% 
    mutate(week = floor_date(date, "week", week_start = 1),
           id = genbankAccession)
  
  if (n != -1) {
    
    # Include and exclude specific sequences by id
    include <- if (include_file != "None") {
      read_lines(include_file)
      include_df <- df %>% filter(id %in% include)
    } else include_df <- tibble()
    exclude <- if (exclude_file != "None") {
      read_lines(exclude_file)
      df <- df %>%  
        filter(!id %in% exclude)
    }
    
    # Subsample
    if (method == "random") {
      # 1. Random in the full period
      subsample <- df %>%
        sample_n(size = min(n(), n), replace = F)
    
    } else if (method == "cases") {
        
      # 2. Proportional to cases in each week
      cases <- read.csv(cases_file)
      cases_filtered <- cases %>%
        mutate(week = floor_date(ymd(date), "week", week_start = 1)) %>%
        filter(ymd(date) > min(ymd(df$date)), ymd(date) < max(ymd(df$date)), 
               country %in% unique(df$country))
      p_week <- cases_filtered %>%
        group_by(week) %>%
        summarise(cases_week = sum(cases), .groups = "drop") %>%
        mutate(p_cases = cases_week/sum(cases_week)) %>% ungroup
      df_prop <- left_join(df, p_week, by = c("week")) %>%
        replace_na(list(p_cases = 0))
      print(p_week)
      print(head(cases))
      print(head(cases_filtered))
      subsample <- df_prop %>%
        sample_n(size = min(n(), n), replace = F, weight = p_cases)
      
      
    } else if (method == "uniform") {
        
      # 3. Uniform in each week
      n_weeks <- nrow(df %>% count(week = floor_date(date, "week", week_start = 1)))
      n_seq_week <- ceiling(n/n_weeks)
      df_seq_week <- df %>% count(week) %>%
        rowwise() %>%
        mutate(n_subsample = min(n, n_seq_week)) %>%
        right_join(df, by = "week")
      subsample <- df_seq_week %>%
        group_by(week) %>%
        sample_n(size = min(n(), n_subsample), replace = F) %>%
        ungroup()
    }
    
    subsample <- bind_rows(subsample, include_df) %>%
      distinct()
  }
  else subsample <- df
  
  # Plot
  gp <- ggplot(subsample %>% count(date, week)) +
    geom_bar(aes(date, n, fill = factor(week)), stat = "identity") +
    labs(title = paste(unique(df$deme), "subsampled sequences, method:", method))
  
  subsample_to_save <- subsample %>%
    mutate(seq_name = paste0(!!sym(seq_id), "|", deme, "|", date)) %>%
    select(seq_id, seq_name)
  
  # Save plot and subsample
  ggsave(gsub(".tsv", "_plot.pdf", output_file), gp, device = "pdf", width = 5)
  write_tsv(subsample_to_save, file = output_file)
  
}


# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(argparse)

# Parser -----------------------------------------------------------------------
parser <- argparse::ArgumentParser()
parser$add_argument("--metadata_file", type = "character", 
                    help = "metadata tsv file with country and date")
parser$add_argument("--include_file", type = "character", 
                    help = "include file with strain ids")
parser$add_argument("--exclude_file", type = "character", 
                    help = "exclude file with strain ids")
parser$add_argument("--n", type = "integer",
                    help = "number of seqs to subsample")
parser$add_argument("--method", type = "character",
                    help = "type of subsampling")
parser$add_argument("--seq_id", type = "character",
                    help = "sequence identifier")
parser$add_argument("--cases_file", type = "character",
                    help = "Cases file")
parser$add_argument("--seed", type = "integer")
parser$add_argument("--output_file", type = "character",
                    help = "Output file for subsample")

args <- parser$parse_args()

# Subsampling ------------------------------------------------------------------
subsample_output <- subsample(args$metadata_file, 
                              args$include_file,
                              args$exclude_file,
                              args$n, 
                              args$method, 
                              args$seq_id,
                              args$cases_file,
                              args$seed, 
                              args$output_file) 

#TODO set color for histogram


