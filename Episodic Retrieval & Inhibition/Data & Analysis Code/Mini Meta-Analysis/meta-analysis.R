#------------------------------------------------------------------------------
### set up

# clear workspace
rm(list = ls())

# set working directory
setwd("~/Git/lab-book/Agi's PhD/Episodic Retrieval/Mini Meta-Analysis")

# load data
exp1 <- read.csv("experiment_1_respRep.csv")
exp1 <- exp1[, 1]

exp2_arrows <- read.csv("experiment_2_respRep_arrows.csv")
exp2_arrows <- exp2_arrows[, 1]

exp2_shapes <- read.csv("experiment_2_respRep_shapes.csv")
exp2_shapes <- exp2_shapes[, 1]

exp3_match <- read.csv("experiment_3_respRep_match.csv")
exp3_match <- exp3_match[, 1]

exp3_mismatch <- read.csv("experiment_3_respRep_mismatch.csv")
exp3_mismatch <- exp3_mismatch[, 1]
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### get means, SDs, and N for each condition and plug into ESCI Excel file
### meta-analyses can be done in R, but I haven't learned how yet
mean(exp1)
sd(exp1)
length(exp1)

mean(exp2_arrows)
sd(exp2_arrows)
length(exp2_arrows)

mean(exp2_shapes)
sd(exp2_shapes)
length(exp2_shapes)

mean(exp3_match)
sd(exp3_match)
length(exp3_match)

mean(exp3_mismatch)
sd(exp3_mismatch)
length(exp3_mismatch)
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### construct the forest plot. Thanks to this resource:
### https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html
library(forestplot)

# enter the data from the meta-analysis
meta_data <- structure(list(
  
  mean = c(NA, 48.13, 11.83, 26.85, 39.62, 16.33, NA, 28.42), 
  lower = c(NA, 25.26, -5.11, -11.12, -1.72, -40.17, NA, 10.76), 
  upper = c(NA, 71.00, 28.77, 64.84, 80.96, 72.83, NA, 46.09)), 
  .Names = c("mean", "lower", "upper"), 
  row.names = c(NA, -6L),
  class = "data.frame")

# enter the details of the plot
table_text <- cbind(
  
  c("Experiment", "1", "2: Arrows", "2: Shapes", "3: Match", "3: Mismatch", 
    NA, "Meta Estimate"), 

  c("Weight (%)", "27.2", "34.1", "15.4", "14.4", "8.9", NA, NA), 
  
  c("Cost (ms)", "48", "12", "27", "40", "16", NA, "28")
)

# generate the plot
pdf("meta_analysis.pdf", width = 8, height = 6)
forestplot(table_text, 
           txt_gp = fpTxtGp(ticks = gpar(fontfamily = "", 
                                         cex = 1), 
                            xlab = gpar(fontfamily = "", 
                                        cex = 1)),
           hrzl_lines = list("2" = gpar(lty = 1), 
                             "8" = gpar(lwd=1, columns=1:3, col = "#000044")),
           meta_data,
           new_page = FALSE,
           is.summary = c(TRUE, rep(FALSE, 7)), 
           col = fpColors(box = "royalblue", line = "darkblue", 
                          summary = "royalblue", zero = "darkblue"), 
           xlab = "N-2 Task Repetition Cost (ms)",
           lwd.zero = 1.5,
           vertices = T)
dev.off()
#------------------------------------------------------------------------------