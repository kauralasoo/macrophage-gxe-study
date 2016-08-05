library("readr")
library("dplyr")
library("ggplot2")

#Import cell counts from disk
colnames = c("line_id", "start_date", "date", "harvest", "cell_count", "comment")
counts = readr::read_delim("macrophage-gxe-study/data/sample_lists/macrophage_counts_tidy.txt", 
                           delim = "\t", col_names = colnames, col_types = "ccccdc") %>%
  dplyr::mutate(start_date = as.Date(start_date,"%d/%m/%y")) %>%
  dplyr::mutate(date = as.Date(date, "%d/%m/%y")) %>%
  tbl_df() %>%
  dplyr::mutate(harvest = ifelse(is.na(harvest),"No",harvest)) %>%
  dplyr::mutate(cell_count = ifelse(is.na(cell_count),0,cell_count)) %>%
  dplyr::mutate(diff_day = as.numeric(date - start_date))

#Calculate some count summaries
count_summaries = dplyr::filter(counts, cell_count > 0) %>% 
  dplyr::group_by(line_id, start_date) %>% 
  dplyr::summarise(measurement_count = length(cell_count), mean_count = mean(cell_count), 
                   total_count = sum(cell_count), max_count = max(cell_count),
                   median_count = median(cell_count)) %>% 
  ungroup()

#Make a plot of the distribution of mean cell count per line
median_count = median(count_summaries$mean_count)/1e6
mean_plot = ggplot(count_summaries, aes(x = mean_count/1e6)) + 
  geom_histogram(binwidth = 1) +
  theme_light() +
  xlab("Mean cells per harvest (millions)") + 
  ylab("Number of lines") +
  geom_vline(xintercept = median_count, color = "red")
ggsave("figures/supplementary/diff_mean_cell_count.pdf", plot = mean_plot, width = 4.5, height = 4.5)


#Explore the trajectories of cell counts
multiple_measurements = dplyr::filter(count_summaries, measurement_count >= 3)

#Remove intermediate zero counts
cumulative_counts = dplyr::semi_join(counts, multiple_measurements, by = c("line_id", "start_date")) %>%
dplyr::group_by(line_id, start_date) %>% 
  dplyr::mutate(cum_count = cumsum(cell_count)) %>%
  dplyr::mutate(max_count = max(cell_count)) %>%
  dplyr::filter(!(cell_count == 0 & cum_count > 0)) %>% 
  dplyr::filter(max_count > 0) %>%
  dplyr::mutate(std_count = (cell_count - mean(cell_count))/sd(cell_count)) %>%
  dplyr::mutate(std_count = std_count - min(std_count)) %>%
  ungroup() %>%
  dplyr::mutate(rel_count = cell_count/max_count)

diff_trajectories = ggplot(cumulative_counts, aes(x = diff_day, y = rel_count, group = paste(line_id, start_date))) +
  geom_point() + geom_line() + 
  geom_smooth(group = 1) +
  theme_light() + 
  ylab("Standardised cell count") +
  xlab("Days from differentiation start")
ggsave("figures/supplementary/diff_trajectories.pdf", plot = diff_trajectories, width = 5.5, height = 4.5)


#How long does it take to harvest first batch of cells
first_harvest_day = dplyr::group_by(counts,line_id, start_date) %>% 
  dplyr::mutate(cum_count = cumsum(cell_count)) %>%
  dplyr::mutate(max_count = max(cell_count)) %>%
  dplyr::filter(!(cell_count == 0 & cum_count > 0)) %>% 
  dplyr::filter(max_count > 0) %>%
  dplyr::filter(cell_count > 0) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::mutate(relative_count = cell_count/max_count) %>%
  dplyr::select(line_id, start_date, cell_count, max_count, relative_count, diff_day) %>%
  ungroup()

#Make a histogram
median_day = median(first_harvest_day$diff_day)
first_harvest_hist = ggplot(first_harvest_day, aes(x = diff_day)) + 
  geom_histogram(binwidth = 3) +
  xlab("Days until first harvest") +
  ylab("Number of lines") +
  theme_light() +
  geom_vline(xintercept = median_day, color = "red")
ggsave("figures/supplementary/diff_first_harvest.pdf", plot = first_harvest_hist, width = 4.5, height = 4.5)


