#####################################################
# Name: exp9_analysis.R
# Description: Analyze experiment7 where we 
#              build taxonomic compressed docprofiles
#
# Date: April 4th, 2022
#####################################################

library(ggplot2)
library(ggpubr)
library(data.table)
library(tidyr)

data_folder <- "/Users/omarahmed/Downloads/test.csv"

df <- read.csv(data_folder, header=FALSE)
colnames(df) <- c("numcols", "table", "overflow", "ktable")
df["total"] <- df["table"] + df["overflow"]

df_pivot <- df %>% pivot_longer(cols=!numcols, names_to="component", values_to="size")

ggplot(df_pivot, aes(x=numcols, y=size, color=component)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  #scale_y_continuous(trans='log2') +
  geom_hline(yintercept=160217920, linetype="dashed", color = "red") +
  scale_color_discrete(name="Index Component", labels=c("Top-k Table", "Overflow", "Table", "Table+Overflow")) +
  labs(x="Number of Columns in Table", y="Size (bytes)", title="Analysis for leftmost 100 genera in SILVA") +
  scale_x_continuous(breaks=seq(2, 16, 1))