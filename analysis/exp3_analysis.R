#####################################################
# Name: exp3_analysis.R
# Description: Analyze experiment 3 where we 
#              compare the classification of SPUMONI
#              and document array profiles on AMR 
#              genes
#
# Date: October 5th, 2022
#####################################################

library(ggplot2)
library(ggpubr)
library(data.table)
library(tidyr)

########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
########################################################################

data_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_3/trial_2/data/"
plots_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_3/trial_2/plots/"

########################################################################
# Methods for plots ...
########################################################################

make_barplot_comparison <- function(total_df) {
  
  dataset.labs <- c("0"="Class A", "1"="Class B", "2"="Class C", "3"="Class D")
  bar.labels <- c("sensitivity"="Sensitivity", "specificity"="Specificity")
  tool.labels <- c("docprofile"="Doc. Array", "spumoni"="SPUMONI 2")
  
  plot <- ggplot(data=total_df, aes(x=metric, y=metric_val)) +
          geom_bar(position="dodge", stat="identity", width=0.5, aes(color=tool, fill=tool)) +
          facet_wrap(~dataset, ncol=2, nrow=2, labeller = labeller(dataset=dataset.labs)) +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                strip.text=element_text(face="bold"),
                axis.title.x=element_text(size=10),
                axis.title.y=element_text(size=10),
                legend.position = "right", 
                #legend.position = c(0.6, 0.8),
                legend.text=element_text(size=10),
                legend.box="vertical",
                legend.background = element_rect(size=0.5, linetype="solid", colour ="black"),
                legend.title=element_text(size=10, face="bold"),
                axis.text=element_text(size=10, color="black")) +
          #scale_color_discrete(name="Number of Classes:", labels=c("3", "5", "8")) +
          #scale_shape_discrete(name="Approach:", labels=c("Doc. Array", "Doc. Array (optimized)", "r-index")) +
          coord_cartesian(ylim = c(0.75,1.00)) +
          labs(x="Metric", y="Value") +
          scale_x_discrete(labels=bar.labels) +
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
merge_df["sensitivity"] <- merge_df["tp"]/(merge_df["tp"] + merge_df["fn"])
merge_df["specificity"] <- merge_df["tn"]/(merge_df["tn"] + merge_df["fp"])

# Rearrange the dataset for plotting purposes ...
wide_df<- merge_df %>% pivot_longer(cols=c('sensitivity', 'specificity'),
                          names_to='metric',
                          values_to='metric_val') 

# Make the plot and save the plot ...
plot <- make_barplot_comparison(wide_df)
plot

output_name <- paste(plots_folder, "comparison_plot.jpeg", sep="")
ggsave(output_name, plot=plot, dpi=800, device="jpeg", width=6, height=4)
