source('functions/getData.R')

# gets all data (testing)
aldh1a1_all <- get_all_data('ALDH1A1')

# get data for tcga only
aldh1a1 <- get_tcga('ALDH1A1')

# each data set includes
# X -- expression matrix for gene
# Y -- processed phenotype with (possibly) risk, time, death, and gender
# p -- complete phenotye data
# NAME -- the name of the dataset
# NOTES -- if included, notes about survival information
names(aldh1a1)

# boxplot across risk groups
plot_favorable(x = as.double(aldh1a1$X), 
               y = aldh1a1$Y$risk, 
               title = paste0(aldh1a1$NAME, ', ALDH1A1'))

# km curve
plot.shiny.km(time = aldh1a1$Y$time,
              death = aldh1a1$Y$death,
              x = as.double(aldh1a1$X),
              title = paste0(aldh1a1$NAME, ', ALDH1A1'),
              continuous = FALSE,
              col = c('darkblue', 'darkred'))
