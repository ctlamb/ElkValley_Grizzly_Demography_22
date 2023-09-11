library(tidyverse)
library(lubridate)


# Create a "date bin" data.frame with columns c('start','end')
# start_date, and last_date are dates
# time_bin_size should be string like '1 year','6 month', etc. see ?seq.Date for more examples
bin_dates = function(first_date, last_date, time_bin_size='1 year'){
  if(!(lubridate::is.Date(first_date) & lubridate::is.Date(last_date))) stop('first_date and last_date must be dates')
  if(first_date>=last_date) stop('first_date must come before last date')
  
  full_range = seq.Date(first_date,last_date, by=time_bin_size)
  starts = full_range[-length(full_range)]
  ends   = full_range[-1]
  
  return(data.frame(start=starts, end=ends))
}

# Takes a data.frame 'df' with columns c('id','start','end','dead')
# time_bin_size should be string like '1 year','6 month', etc. see ?seq.Date for more
# returns a data.frame where each unique ID
stretch_survival_data = function(df, time_bin_size){
  dead_animals = df %>%
    filter(event==1) %>%
    pull(id_period)
  
  # apply the appropriate binning
  # here group_by() %>% summarize() will submit each unique ID to the bin_dates() function
  # and build the appropriate data.frame with all the new rows.
  stretched_df = df %>% 
    group_by(id_period, id, sex) %>% 
    summarise(bin_dates(first_date = start, last_date=end, time_bin_size)) %>%
    ungroup() 
  
  # mark dead animals as such (1) in their last respective timestep, othwerwise mark as
  # alive (0)
  stretched_df = stretched_df%>%
    group_by(id_period, id, sex) %>%
    mutate(dead = case_when(
      id_period %in% dead_animals & end==max(end) ~ 1,
      TRUE ~ 0
    )) %>%
    ungroup()
  
  return(stretched_df)
}
