## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PupillometryR)

## -----------------------------------------------------------------------------
data("pupil_data")

#Check that IDs are not numeric
pupil_data$ID <- as.character(pupil_data$ID)
#remove participant number 8, who had problematic data
pupil_data <- subset(pupil_data, ID != 8)
#blinks were registered as -1, so replace with NAs
pupil_data$LPupil[pupil_data$LPupil == -1] <- NA
pupil_data$RPupil[pupil_data$RPupil == -1] <- NA


## -----------------------------------------------------------------------------
library(ggplot2)
theme_set(theme_classic(base_size = 12))

## -----------------------------------------------------------------------------
Sdata <- make_pupillometryr_data(data = pupil_data,
                                 subject = ID,
                                 trial = Trial,
                                 time = Time,
                                 condition = Type)

## -----------------------------------------------------------------------------
new_data <- replace_missing_data(data = Sdata)

head(new_data)

## -----------------------------------------------------------------------------
plot(new_data, pupil = LPupil, group = 'condition')

plot(new_data, pupil = LPupil, group = 'subject') 

## -----------------------------------------------------------------------------
regressed_data <- regress_data(data = new_data,
                               pupil1 = RPupil,
                               pupil2 = LPupil)


## -----------------------------------------------------------------------------
mean_data <- calculate_mean_pupil_size(data = regressed_data, 
                                       pupil1 = RPupil, 
                                       pupil2 = LPupil)

plot(mean_data, pupil = mean_pupil, group = 'subject')

## -----------------------------------------------------------------------------
mean_data <- downsample_time_data(data = mean_data,
                              pupil = mean_pupil,
                              timebin_size = 50,
                              option = 'median')

## -----------------------------------------------------------------------------
missing <- calculate_missing_data(mean_data, 
                                  mean_pupil)
head(missing)

## ----message = T--------------------------------------------------------------
mean_data2 <- clean_missing_data(mean_data,
                                 pupil = mean_pupil,
                                 trial_threshold = .75,
                                 subject_trial_threshold = .75)

## -----------------------------------------------------------------------------
filtered_data <- filter_data(data = mean_data2,
                             pupil = mean_pupil,
                             filter = 'median',
                             degree = 11)

plot(filtered_data, pupil = mean_pupil, group = 'subject')

## -----------------------------------------------------------------------------
int_data <- interpolate_data(data = filtered_data,
                             pupil = mean_pupil,
                             type = 'linear')

plot(int_data, pupil = mean_pupil, group = 'subject')

## -----------------------------------------------------------------------------
base_data <- baseline_data(data = int_data,
                           pupil = mean_pupil,
                           start = 0,
                           stop = 100)

plot(base_data, pupil = mean_pupil, group = 'subject')

## -----------------------------------------------------------------------------
window <- create_window_data(data = base_data,
                             pupil = mean_pupil)

plot(window, pupil = mean_pupil, windows = F, geom = 'boxplot')

head(window)

## -----------------------------------------------------------------------------
t.test(mean_pupil ~ Type, paired = T, data = window)

## -----------------------------------------------------------------------------
timeslots <- create_time_windows(data = base_data,
                                 pupil = mean_pupil,
                                 breaks = c(0, 2000, 4000, 6000, 8000, 10000))

plot(timeslots, pupil = mean_pupil, windows = T, geom = 'raincloud')

head(timeslots)

## -----------------------------------------------------------------------------
lm(mean_pupil ~ Window * Type, data = timeslots)

## -----------------------------------------------------------------------------
library(mgcv)

base_data$IDn <- as.numeric(base_data$ID)
base_data$Typen <- ifelse(base_data$Type == 'Easy', .5, -.5)
base_data$Trialn <- as.numeric(substr(base_data$Trial, 5, 5))
base_data$Trialn <- ifelse(base_data$Typen == .5, base_data$Trialn, base_data$Trialn + 3)
base_data$ID <- as.factor(base_data$ID)
base_data$Trial <- as.factor(base_data$Trial)

## -----------------------------------------------------------------------------
m1 <- bam(mean_pupil ~ s(Time) +
            s(Time,  by = Typen),
          data = base_data,
          family = gaussian)

summary(m1)

## -----------------------------------------------------------------------------
plot(base_data, pupil = mean_pupil, group = 'condition', model = m1)

## -----------------------------------------------------------------------------
qqnorm(resid(m1))

itsadug::acf_resid(m1)

## -----------------------------------------------------------------------------
base_data$Event <- interaction(base_data$ID, base_data$Trial, drop = T)

model_data <- base_data
model_data <- itsadug::start_event(model_data,
                          column = 'Time', event = 'Event')

model_data <- droplevels(model_data[order(model_data$ID,
                                          model_data$Trial,
                                          model_data$Time),])

## -----------------------------------------------------------------------------
m2 <- bam(mean_pupil ~ Typen +
            s(Time,  by = Typen) +
            s(Time, Event, bs = 'fs', m = 1),
          data = base_data,
          family = gaussian,
          discrete = T,
          AR.start = model_data$start.event, rho = .6)

# m2 <- bam(mean_pupil ~ 
          #   s(Time,  by = Typen) +
          #   s(Time, Event, bs = 'fs', m = 1),
          # data = base_data,
          # family = scat,
          # discrete = T,
          # AR.start = model_data$start.event, rho = .6)

summary(m2)
qqnorm(resid(m2))
itsadug::acf_resid(m2)
plot(base_data, pupil = mean_pupil, group = 'condition', model = m2)

## -----------------------------------------------------------------------------
differences <- create_difference_data(data = base_data,
                                      pupil = mean_pupil)

plot(differences, pupil = mean_pupil, geom = 'line')

## -----------------------------------------------------------------------------
spline_data <- create_functional_data(data = differences,
                                      pupil = mean_pupil,
                                      basis = 10,
                                      order = 4)


plot(spline_data, pupil = mean_pupil, geom = 'line', colour = 'blue')

## -----------------------------------------------------------------------------
ft_data <- run_functional_t_test(data = spline_data,
                                 pupil = mean_pupil,
                                 alpha = 0.05)


plot(ft_data, show_divergence = T, colour = 'red', fill = 'orange')

