#####################################################
# Name: exp2_analysis.R
# Description: Analyze experiment 2 where we run queries
#              through both the r-index and the document
#              arrays and compare the index sizes.
# Date: October 3rd, 2022
#####################################################

library(ggplot2)
library(ggpubr)
library(data.table)

########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
########################################################################

data_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_2/trial_3/data/"
plots_folder <- "/Users/omarahmed/downloads/current_research/docarray_exps/exp_2/trial_3/plots/"

########################################################################
# Methods for plots ...
########################################################################

make_comparison_plot <- function(total_df) {
  facet.labs <- c("100"="100+", "300"="300+", "30"="30+")
  plot <- ggplot(data=total_df, aes(x=indexsize/1073741824, y= time)) +
          geom_point(aes(color=as.factor(numclasses), shape=approach), size=3.0) + 
          theme_bw() +
          facet_wrap(~numgenomes, labeller=labeller(numgenomes=facet.labs)) +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                strip.text=element_text(face="bold", size=12),
                axis.title.x=element_text(size=12),
                axis.title.y=element_text(size=12),
                legend.position = "right", 
                #legend.position = c(0.6, 0.8),
                legend.text=element_text(size=12),
                legend.box="vertical",
                legend.background = element_rect(size=0.5, linetype="solid", colour ="black"),
                legend.title=element_text(size=12, face="bold"),
                axis.text=element_text(size=10, color="black")) +
          #guides(shape=guide_legend(override.aes=list(size = 0.5))) +
          scale_y_continuous(breaks=seq(0, 400, 50), limits=c(0, 400)) +
          scale_x_continuous(breaks=seq(0, 8, 1)) +
          scale_color_discrete(name="Number of Classes:", labels=c("3", "5", "8")) +
          scale_shape_discrete(name="Approach:", labels=c("Doc. Array", "Doc. Array (optimized)", "r-index")) +
          labs(x="Index Size (GB)",
               y="Time (sec)",
               title="") 
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


# Make a plot
plot_1 <- make_comparison_plot(merge_df)
plot_1

output_name <- paste(plots_folder, "query_comparison_plot.pdf", sep="")
ggsave(output_name, plot=plot_1, dpi=800, device="pdf", width=10, height=4)
