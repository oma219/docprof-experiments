#####################################################
# Name: exp7_analysis.R
# Description: Analyze experiment7 where we 
#              assess sequence similarity of datasets
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

data_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_7/trial_1/data/"
plots_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_7/trial_1/plots/"

########################################################################
# Methods for plots ...
########################################################################

make_distribution_plot <- function(total_df) {
  plot <- ggplot(data=total_df, aes(x=name, y=ANI)) +
          geom_boxplot(aes(fill=name), color="black") +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                strip.text=element_text(face="bold", size=12),
                axis.title.x=element_text(size=12),
                axis.title.y=element_text(size=12),
                legend.position = "none", 
                legend.text=element_text(size=12),
                legend.box="horizontal",
                legend.title=element_text(size=12, face="bold"),
                axis.text=element_text(size=12, color="black")) +
          scale_x_discrete(limits=c("salmonstrain","ecolistrain","egenus","MOCK"), 
                           labels=c("S. enterica strains", "E. coli strains", "Same Genus", "Different Genera")) +
          coord_flip() +
          labs(x="Dataset", y="Average Nucleotide Identity (ANI)") 
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

# Make the violin plot ...
plot <- make_distribution_plot(merge_df)
plot

output_name <- paste(plots_folder, "ani_result.pdf", sep="")
ggsave(output_name, plot=plot, dpi=800, device="pdf", width=5, height=3)
