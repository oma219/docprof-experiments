#####################################################
# Name: exp6_analysis.R
# Description: Analyze experiment65 where we 
#              compare the classification of SPUMONI
#              and document array profiles on real
#              nanopore data
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

data_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_6/trial_1/data/"
plots_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_6/trial_1/plots/"

########################################################################
# Methods for plots ...
########################################################################

make_grouped_bar_chart <- function(total_df) {
  bar.labels <- c("0"="S. aureus", "1"="S. enterica", "2"="E. coli", "3"="P. aeruginosa", "4"="L. monocytogenes", "5"="E. faecalis", "6"="B. subtilis")
  tool.labels <- c("docprofile"="Doc. Array", "spumoni"="SPUMONI 2")
  
  plot <- ggplot(data=total_df, aes(x=as.factor(dataset), y=recall)) +
          geom_bar(position="dodge", stat="identity", width=0.5, aes(fill=tool), color="black") +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
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
          #coord_cartesian(ylim = c(0.75,1.00)) +
          labs(x="Strain Classification Dataset", y="Recall") +
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

# Make a grouped bar chart ...
plot <- make_grouped_bar_chart(merge_df)
plot

output_name <- paste(plots_folder, "real_classification_test.pdf", sep="")
ggsave(output_name, plot=plot, dpi=800, device="pdf", width=9, height=4)




