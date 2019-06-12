#!/usr/local/apps/R/R-3.4.4/bin/ Rscript

library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(optparse)



opt_list <- list(
  make_option(c("-d", "--data")),
  make_option(c("-i", "--info")),
  make_option(c("-nql", "--no_qstat_lines"), default = TRUE, action = "store_false"),
  make_option(c("--name"), default = NA, action="store", type="character"),
  make_option(c("-rs", "--rdata_save"), default=NA, action="store_true")
)
parser <- OptionParser(option_list = opt_list)
parsed_args <- parse_args(parser)

# First argument is the data file
df <- parsed_args$d
data <- read.csv(df, stringsAsFactors = FALSE)
data <- as_tibble(data)
data$datetime <- as.POSIXct(data$datetime, format = "%Y-%m-%d %H:%M:%S", tz="UTC")
data2 <- gather(data, Memory_Type, Memory, resources_used.mem, resources_used.vmem)

info_file <- parsed_args$i
info<-read.csv(info_file, stringsAsFactors = FALSE)
info <- as_tibble(info)
info2 <- gather(info, Time_Type, T, ctime, qtime, mtime, stime, etime)
info2$T <- as.POSIXct(info2$T, format = "%a %b  %d %H:%M:%S %Y", tz="UTC")

center_of_lim <- max(as.numeric(stringr::str_replace(data2$Memory, "kb", "")))/2
center_of_lim <- paste0(center_of_lim, "kb")

qstat_line_flag <- parsed_args$no_qstat_lines
if(qstat_line_flag){
  left_limit <- with_tz(min(c(info2$T, data2$datetime)), "UTC")
  right_limit <- with_tz(max(c(info2$T, data2$datetime)), "UTC")
p <- ggplot(data2, aes(datetime, Memory, color=Memory_Type)) + geom_line() + scale_x_datetime(date_labels="%Y-%m-%d %H:%M:%S", limits = c(left_limit, right_limit)) +
  theme(axis.text.x = element_text(angle=60, hjust=1)) +
  geom_vline(data=info2, aes(xintercept = info2$T), show.legend = TRUE, linetype=4) +
  geom_label_repel(data=info2, aes(x=info2$T, y=center_of_lim, label=info2$Time_Type), inherit.aes = FALSE,  size=4)
} else{
  left_limit <- with_tz(min(c(data2$datetime)), "UTC")
  right_limit <- with_tz(max(c(data2$datetime)), "UTC")
  p <- ggplot(data2, aes(datetime, Memory, color=Memory_Type)) + geom_line() + scale_x_datetime(date_labels="%Y-%m-%d %H:%M:%S", limits = c(left_limit, right_limit)) +
    theme(axis.text.x = element_text(angle=60, hjust=1))
}

rdata_save <- parsed_args$rdata_save
new_name <- parsed_args$name
if(!is.na(new_name)){
  if(!is.na(rdata_save)){
    save(p, file = sprintf("%s.RData", new_name))
  } else{
  ggsave(plot=p,filename = sprintf("%s.png", new_name))
  }
  } else {
    if(!is.na(rdata_save)){
      save(p, file = sprintf("%s.RData", df))
    } else {
  ggsave(plot=p,filename = sprintf("%s.png", df))
  }
  }

