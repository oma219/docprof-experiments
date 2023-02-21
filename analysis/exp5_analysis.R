#####################################################
# Name: exp5_analysis.R
# Description: Analyze experiment 5 where we 
#              compare the classification of SPUMONI
#              and document array profiles on genomic
#              datasets
#
# Date: October 13th, 2022
#####################################################

library(ggplot2)
library(ggpubr)
library(data.table)
library(tidyr)

########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
########################################################################
base_dir <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_5/trial_1/"
data_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_5/trial_1/ecoligenus/data/"
plots_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_5/trial_1/ecoligenus/plots/"

########################################################################
# Methods for plots ...
########################################################################

make_barplot_comparison <- function(total_df) {
  
  # Salmonella labels ...
  #dataset.labs <- c("0"="GCF_025398995.1", "1"="GCF_025399055.1", "2"="GCF_025399075.1", "3"="GCF_025399015.1")
  
  # E. coli labels ...
  #dataset.labs <- c("0"="GCF_025426235.1", "1"="GCF_025563475.1", "2"="GCF_025563515.1", "3"="GCF_025563435.1")
  
  # Mock Community labels ...
  #dataset.labs <- c("0"="Escherichia coli", "1"="Salmonella enterica", "2"="Listeria monocytogenes", "3"="Pseudomonas aeruginosa")
  
  # Bacillus labels ...
  #dataset.labs <- c("0"="Bacillus cereus", "1"="Bacillus anthracis", "2"="Bacillus thuringiensis", "3"="Bacillus myocides")
  
  # E. coli genus labels ...
  dataset.labs <- c("0"="Escherichia coli", "1"="Escherichia albertii", "2"="Escherichia fergusonii", "3"="Escherichia marmotae")
  
  bar.labels <- c("recall"="Recall", "precision"="Precision")
  tool.labels <- c("docprofile"="Doc. Array", "spumoni"="SPUMONI 2")
  
  plot <- ggplot(data=total_df, aes(x=metric, y=metric_val)) +
    geom_bar(position="dodge", stat="identity", width=0.5, aes(color=tool, fill=tool)) +
    facet_wrap(~dataset, ncol=4, nrow=1, labeller = labeller(dataset=dataset.labs)) +
    theme_bw() +
    theme(plot.title=element_text(hjust = 0.5, size=12, face="bold"),
          strip.text=element_text(face="bold", size=12),
          axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.position = "bottom", 
          #legend.position = c(0.6, 0.8),
          legend.text=element_text(size=12),
          legend.box="horizontal",
          #legend.background = element_rect(size=0.5, colour ="black"),
          legend.title=element_text(size=12, face="bold"),
          axis.text=element_text(size=10, color="black")) +
    #scale_color_discrete(name="Number of Classes:", labels=c("3", "5", "8")) +
    #scale_shape_discrete(name="Approach:", labels=c("Doc. Array", "Doc. Array (optimized)", "r-index")) +
    #coord_cartesian(ylim = c(0.75,1.00)) +
    labs(x="Metric", y="Value") +
    scale_x_discrete(labels=bar.labels) +
    scale_y_continuous(breaks=seq(0, 1.0, 0.10)) +
    geom_hline(yintercept=0.25, linetype="dashed", color = "red") +
    scale_color_discrete(name="Tools", labels=tool.labels) +
    scale_fill_discrete(name="Tools", labels=tool.labels) 
  return (plot)
}

make_grid_barplot_comparison <- function(total_df) {
  dataset.labs <- c("0"="Class 1", "1"="Class 2", "2"="Class 3", "3"="Class 4")
  #ecoli_strain_label <- subtitute(paste(italic("S. enterica"), " strains"))
  task.labs <- c("Different Genera"="Different Genera", 
                 "Same Genus"="Same Genus", 
                 "E. coli strains"="ecoli_strain_label", 
                 "S. enterica strains"= paste("hello", "world"))
  bar.labels <- c("recall"="Recall", "precision"="Precision")
  tool.labels <- c("docprofile"="Doc. Array", "spumoni"="SPUMONI 2")
  
  #total_df$task <- as.factor(total_df$task)
  #levels(total_df$task) <- c("Different Genera"="task1", "Same Genus"="task2", "E. coli strains"="task3", "S. enterica strains"="task4")
  
  total_df$task <- factor(total_df$task,    # Change factor labels
                          levels = c("Different Genera", "Same Genus", "E. coli strains", "S. enterica strains"),
                          labels = c("paste('Different', ' Genera')", "paste('Same', ' Genus')", "paste(italic('E. coli'), ' strains')", "paste(italic('S. enterica'), ' strains')"))

  
  plot <- ggplot(data=total_df, aes(x=metric, y=metric_val)) +
    geom_bar(position="dodge", stat="identity", width=0.5, aes(fill=tool), color="black") +
    facet_grid(task ~ dataset, labeller=labeller(dataset=dataset.labs, task=label_parsed)) +
    #facet_grid(task ~ dataset, labeller=labeller(dataset=dataset.labs, task=task.labs)) +
    theme_bw() +
    theme(plot.title=element_text(hjust = 0.5, size=12, face="bold"),
          strip.text=element_text(face="bold", size=12),
          axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.position = "bottom", 
          #legend.position = c(0.6, 0.8),
          legend.text=element_text(size=12),
          legend.box="horizontal",
          #legend.background = element_rect(size=0.5, colour ="black"),
          legend.title=element_text(size=12, face="bold"),
          axis.text=element_text(size=10, color="black")) +
    #scale_color_discrete(name="Number of Classes:", labels=c("3", "5", "8")) +
    #scale_shape_discrete(name="Approach:", labels=c("Doc. Array", "Doc. Array (optimized)", "r-index")) +
    #coord_cartesian(ylim = c(0.75,1.00)) +
    labs(x="Metric", y="Value") +
    scale_x_discrete(labels=bar.labels) +
    scale_y_continuous(breaks=seq(0, 1.0, 0.10)) +
    geom_hline(yintercept=0.25, linetype="dashed", color = "red") +
    scale_color_discrete(name="Tools", labels=tool.labels) +
    scale_fill_discrete(name="Tools", labels=tool.labels) 
  return (plot)
}


########################################################################
# Main Method Code
########################################################################

data_files <- list.files(path=data_folder, pattern = "\\.csv$")
df_list = list()

# Build a total dataframe with all data
for(i in seq_along(data_files)) {
  # Read in csv and annotate with name and dataset
  temp_df <- read.csv(file=paste(data_folder, data_files[i], sep = ""), header=TRUE)
  df_list[[data_files[i]]] <- temp_df
}
merge_df <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)

# Generate classification metrics based on raw data ...
merge_df["recall"] <- merge_df["tp"]/(merge_df["tp"] + merge_df["fn"])
merge_df["precision"] <- merge_df["tp"]/(merge_df["tp"] + merge_df["fp"])

# Rearrange the dataset for plotting purposes ...
wide_df<- merge_df %>% pivot_longer(cols=c('recall', 'precision'),
                                    names_to='metric',
                                    values_to='metric_val') 

# Make the plot and save the plot ...
plot <- make_barplot_comparison(wide_df)
plot

output_name <- paste(plots_folder, "comparison_plot.jpeg", sep="")
ggsave(output_name, plot=plot, dpi=800, device="jpeg", width=10, height=4)


########################################################################
# Main Method Code - Make a grid plot ...
########################################################################

# Collection of data-frames ...
total_df_list = list()

###############################################################
# Dataset 1: Mock Community
###############################################################

data_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_5/trial_1/mock/data/"
plots_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_5/trial_1/mock/plots/"

data_files <- list.files(path=data_folder, pattern = "\\.csv$")
df_list = list()

# Build a total dataframe with all data
for(i in seq_along(data_files)) {
  # Read in csv and annotate with name and dataset
  temp_df <- read.csv(file=paste(data_folder, data_files[i], sep = ""), header=TRUE)
  df_list[[data_files[i]]] <- temp_df
}
merge_df <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)

# Generate classification metrics based on raw data ...
merge_df["recall"] <- merge_df["tp"]/(merge_df["tp"] + merge_df["fn"])
merge_df["precision"] <- merge_df["tp"]/(merge_df["tp"] + merge_df["fp"])

# Rearrange the dataset for plotting purposes ...
wide_df<- merge_df %>% pivot_longer(cols=c('recall', 'precision'),
                                    names_to='metric',
                                    values_to='metric_val') 

# Replace dataset names ...
# dataset.labs <- c("0"="Escherichia coli", "1"="Salmonella enterica", "2"="Listeria monocytogenes", "3"="Pseudomonas aeruginosa")
# 
# for (x in unique(wide_df$dataset)) {
#   wide_df$dataset[wide_df$dataset == x] <- dataset.labs[[as.character(x)]]
#   print(x)
#   print(dataset.labs[[as.character(x)]])
# }

wide_df["task"] <- "Different Genera"
total_df_list[[data_folder]] <- wide_df


#################################################################
# Dataset 2: Escherichia genus
#################################################################

data_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_5/trial_1/ecoligenus/data/"
plots_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_5/trial_1/ecoligenus/plots/"

data_files <- list.files(path=data_folder, pattern = "\\.csv$")
df_list = list()

# Build a total dataframe with all data
for(i in seq_along(data_files)) {
  # Read in csv and annotate with name and dataset
  temp_df <- read.csv(file=paste(data_folder, data_files[i], sep = ""), header=TRUE)
  df_list[[data_files[i]]] <- temp_df
}
merge_df <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)

# Generate classification metrics based on raw data ...
merge_df["recall"] <- merge_df["tp"]/(merge_df["tp"] + merge_df["fn"])
merge_df["precision"] <- merge_df["tp"]/(merge_df["tp"] + merge_df["fp"])

# Rearrange the dataset for plotting purposes ...
wide_df<- merge_df %>% pivot_longer(cols=c('recall', 'precision'),
                                    names_to='metric',
                                    values_to='metric_val') 

# Replace dataset names ...
# dataset.labs <- c("0"="Escherichia coli", "1"="Escherichia albertii", "2"="Escherichia fergusonii", "3"="Escherichia marmotae")
# 
# for (x in unique(wide_df$dataset)) {
#   wide_df$dataset[wide_df$dataset == x] <- dataset.labs[[as.character(x)]]
#   print(x)
#   print(dataset.labs[[as.character(x)]])
# }

wide_df["task"] <- "Same Genus"
total_df_list[[data_folder]] <- wide_df

######################################################################
# Dataset 3: Escherichia coli strains
######################################################################

data_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_5/trial_1/ecoli/data/"
plots_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_5/trial_1/ecoli/plots/"

data_files <- list.files(path=data_folder, pattern = "\\.csv$")
df_list = list()

# Build a total dataframe with all data
for(i in seq_along(data_files)) {
  # Read in csv and annotate with name and dataset
  temp_df <- read.csv(file=paste(data_folder, data_files[i], sep = ""), header=TRUE)
  df_list[[data_files[i]]] <- temp_df
}
merge_df <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)

# Generate classification metrics based on raw data ...
merge_df["recall"] <- merge_df["tp"]/(merge_df["tp"] + merge_df["fn"])
merge_df["precision"] <- merge_df["tp"]/(merge_df["tp"] + merge_df["fp"])

# Rearrange the dataset for plotting purposes ...
wide_df<- merge_df %>% pivot_longer(cols=c('recall', 'precision'),
                                    names_to='metric',
                                    values_to='metric_val') 

# Replace dataset names ...
# dataset.labs <- c("0"="GCF_025426235.1", "1"="GCF_025563475.1", "2"="GCF_025563515.1", "3"="GCF_025563435.1")
# 
# for (x in unique(wide_df$dataset)) {
#   wide_df$dataset[wide_df$dataset == x] <- dataset.labs[[as.character(x)]]
#   print(x)
#   print(dataset.labs[[as.character(x)]])
# }

wide_df["task"] <- "E. coli strains"
total_df_list[[data_folder]] <- wide_df

##########################################################################
# Dataset 4: Salmonella enterica strains
##########################################################################

data_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_5/trial_1/salmonella/data/"
plots_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_5/trial_1/salmonella/plots/"

data_files <- list.files(path=data_folder, pattern = "\\.csv$")
df_list = list()

# Build a total dataframe with all data
for(i in seq_along(data_files)) {
  # Read in csv and annotate with name and dataset
  temp_df <- read.csv(file=paste(data_folder, data_files[i], sep = ""), header=TRUE)
  df_list[[data_files[i]]] <- temp_df
}
merge_df <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)

# Generate classification metrics based on raw data ...
merge_df["recall"] <- merge_df["tp"]/(merge_df["tp"] + merge_df["fn"])
merge_df["precision"] <- merge_df["tp"]/(merge_df["tp"] + merge_df["fp"])

# Rearrange the dataset for plotting purposes ...
wide_df<- merge_df %>% pivot_longer(cols=c('recall', 'precision'),
                                    names_to='metric',
                                    values_to='metric_val') 

# Replace dataset names ...
# dataset.labs <- c("0"="GCF_025426235.1", "1"="GCF_025563475.1", "2"="GCF_025563515.1", "3"="GCF_025563435.1")
# 
# for (x in unique(wide_df$dataset)) {
#   wide_df$dataset[wide_df$dataset == x] <- dataset.labs[[as.character(x)]]
#   print(x)
#   print(dataset.labs[[as.character(x)]])
# }

wide_df["task"] <- "S. enterica strains"
total_df_list[[data_folder]] <- wide_df

###############################################
# Merge the plots, and make the grid plot
###############################################

merge_df <- Reduce(function(x, y) merge(x, y, all=TRUE), total_df_list)

plot <- make_grid_barplot_comparison(merge_df)
plot

output_name <- paste("/Users/omarahmed/downloads/current_research/docarray_exps/exp_5/", "sim_classification_task.pdf", sep="")
ggsave(output_name, plot=plot, dpi=800, device="pdf", width=8, height=8)

