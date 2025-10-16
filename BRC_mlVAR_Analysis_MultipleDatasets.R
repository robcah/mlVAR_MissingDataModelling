# packages (first install packages through "install.packages()")

library("tidyverse")
library("psych")
library("tseries")
library("mlVAR")
library("glue")
library("RColorBrewer")
library("nortest")
library("qgraph")

### Setting up the output folder
out_folder = '.'

# Assumption checking and mlVAR analysis for original data.
BRC_data <- read.csv('BRC_data_wide_BlanchardSchumacherProtocol.csv'
                     , header=T
                     , sep=","
) 

BRC_data$encounter_date <- as.Date(BRC_data$encounter_date,
)
BRC_data %>% head()

### Factor Analysis
factor_data <- (BRC_data 
                %>% dplyr::select(where(is.numeric) 
                                  & !participant_id
                                  & !day_n
                )
)

fa.parallel(factor_data, 
            fm="ml",
)
efaResPos <- (factor_data 
              %>% fa(nfactors=4,
                     rotate="oblimin",
                     fm="ml",
              )
)

fa.diagram(efaResPos)


# Stationarity tests
## Visual inspection

valid_participants=(BRC_data 
                    %>% group_by(participant_id)
                    %>% count()
                    %>% arrange(desc(n))
)$participant_id
constructs=colnames(factor_data)
colours = colorRampPalette(brewer.pal(8
                                      , 'Dark2'
)
)(length(constructs))

for (i in seq_along(constructs))
{
  construct = constructs[i]
  y_label=sprintf("%s rating", construct)
  plot=(BRC_data 
        %>% subset(participant_id
                   %in% valid_participants
        )
        %>% ggplot(aes(x=day_n,
                       y=.data[[construct]]
        )
        )
        + geom_line(color=colours[i])
        + scale_x_continuous("Date time")
        + scale_y_continuous(construct)
        + facet_wrap(~participant_id
                     , ncol=6
                     , shrink=FALSE
        )
        + ylab(y_label)
        + theme(panel.grid.minor=element_blank() 
                , panel.grid.major=element_blank()
                , strip.text=element_text(size=5
                                          , margin=margin()
                )
                , axis.text.x=element_text(size=6)
                , axis.text.y=element_text(size=6)
        )
  )
  print(plot)
}

## Kwiatkowski-Phillips-Schmidt-Shin test is used for testing a null hypothesis

### that an observable time series is stationary around a deterministic trend 
### (i.e. trend-stationary) against the alternative of a unit root
### Indicated by the parameter null="Trend" of tseries::kpss.test.
### H0: trend-stationary, H1: non-stationary

constructs = colnames(factor_data)
participants = unique(BRC_data$participant_id)
m = length(participants)

for (construct in constructs)
{
  kpps_p_construct <- c()
  for(p in participants)
  {
    BRC_construct_participant = (BRC_data
                                 [[construct]]
                                 [BRC_data$participant_id==p]
    )
    randomness <- tseries::kpss.test(na.exclude(BRC_construct_participant)
                                     , lshort=TRUE
                                     , null="Trend"
    )
    kpps_p_construct <- c(kpps_p_construct,
                          randomness$p.value
    )
  }
  kpps_p_construct  <- cbind(participants
                             , kpps_p_construct
  )
  cat(paste('KPSS test Bonferroni(alpha<0.05,'
            , glue('m={m}) for {construct},')
            , 'is it trend-stationary?:'
            , '\n'
            , glue('{all(kpps_p_construct[,2] < 0.05/m)}')
            , '\n'
  )
  )
}

# Normality

## Visual inspection
### Histogram vs Normal distribution

constructs=colnames(factor_data)
colours=colorRampPalette(brewer.pal(8
                                    , 'Dark2'
                                    )
)(length(constructs))
png("HistvNorm_OriginalData.png"
    , width=800
    , height=400
)
par(mfrow=c(2, 4))
for (i in seq_along(constructs))
{
  construct <- constructs[[i]]
  g <- density(BRC_data[[construct]])
  
  plot(g
       , ylab='density'
       , xlab="values"
       , main=construct
  )
  polygon(g, col=colours[[i]])
  xfit <- seq(min(g$x)
              , max(g$x)
              , length=50
  ) 
  yfit <- dnorm(xfit
                , mean=mean(g$x)
                , sd=sd(g$x)
  )
  lines(xfit
        , yfit
        , col='black'
        , lwd=2
  )
}
dev.off()

### QQ Plot
png("QQPlot_OriginalData.png"
    , width=800
    , height=400
)
par(mfrow=c(2, 4))
for (i in seq_along(constructs))
{
  construct <- constructs[[i]]
  qqnorm(BRC_data[[construct]]
         , main=construct
         , col=colours[[i]]
         )
  qqline(BRC_data[[construct]])
}
dev.off()



## Kolmogorov-Smirnov: Test to assess whether the data is sampled from a normal 
### distribution, namely we assume normality.
### H0: normality, H1: non-normality
### H0 is true if p-value > 0.05

for (construct in constructs)
{
  ks_construct = ks.test(BRC_data
                         [[construct]]
                         , rnorm(length(BRC_data
                                        [[construct]]
                         )
                         )
  )$p.value
  cat(paste(glue('KS test for {construct},')
            , 'is it normally distributed?:'
            , '\n'
            , glue('{all(ks_construct > 0.05)}')
            , '\n'
  )
  )
}

# ## Shapiro-Wilk test: Test to assess whether the data is sampled from a normal
# ### distribution, namely we assume normality. This test is more robust than
# ### Kolmogorov-Smirnov but it is limited to samples <= 5,000.To overcome this
# ### last limitation 5,000 samples (without replacement) are performed.
# ### H0: normality, H1: non-normality
# ### H0 is true if p-value > 0.05
# 
# for (construct in constructs)
# {
#   sw_construct = shapiro.test(BRC_data[[construct]]
#   )
#   cat(paste(glue('SW test for {construct}, is it normally distributed?:')
#             , '\n'
#             , glue('{all(sw_construct$p.value > 0.05)}')
#             , '\n'
#   )
#   )
# }

## Anderson-Darling normality: This test to assess whether the data is sampled 
### from a normal distribution, namely we assume normality. This test is not as 
### robust as Shapiro-Wilk but it is not limited to samples <= 5,000
### H0: normality, H1: non-normality
### H0 is true if p-value > 0.05

for (construct in constructs)
{
  # print(construct)
  ad_construct = ad.test(BRC_data[[construct]])
  cat(paste(glue('AD test for {construct}, '
                 , 'is it normally distributed?:'
  )
  , '\n'
  , glue('{all(ad_construct$p.value > 0.05)}')
  , '\n'
  )
  )
}

# Main mlVAR analysis
## Detrending: Despite stationarity results seeming ok (at least for 0.05 alpha) we were concerned as participants used the app for such a long time period, wrote a loop which detrends each participant individually and takes residuals.

consec_detrend=(BRC_data
                %>% arrange(dataset,
                            participant_id,
                            day_n
                )
                %>% dplyr::select(dataset,
                                  participant_id,
                                  day_n,
                                  everything()
                )
)

# Sorted by volume
participants=(BRC_data 
              %>% group_by(participant_id)
              %>% count()
              %>% arrange(desc(n)) # _n_ is the counts column
)$participant_id

constructs = colnames(BRC_data)[-1:-4]

# Based on S Allen et al
for (construct in constructs)
{
  for(participant in valid_participants)
  {
    detrended_data <- (lm(na.exclude(consec_detrend[[construct]]
                                     ~ consec_detrend$day_n,
                                     data=consec_detrend
                                     [consec_detrend$participant_id==participant]
    )
    )
    %>% resid()
    )
  }
  consec_detrend[[construct]][!is.na(consec_detrend[[construct]])] <- detrended_data
}

## mlVAR networks analysis

constructs = colnames(factor_data)

network_detrend <- mlVAR(consec_detrend
                         , vars=constructs
                         , idvar="participant_id"
                         , estimator="default"
                         , contemporaneous="correlated"
                         , temporal="correlated"
                         , beepvar="day_n"
                         , lags=1
)

## Plotting network results

for (i in mapply(list,
                 c(0.00625, 0.003125, 0.00625),
                 c('temporal'
                   , 'between'
                   , 'contemporaneous'
                 ),
                 SIMPLIFY=F)
)
{
  a = i[[1]]
  type = i[[2]]
  xi_detrend <- plot(network_detrend
                     , type
                     , layout="circle"
                     , rule="and"
                     , legend=FALSE
                     , edge.labels=T
                     , alpha=a
                     , cut=0.1
                     , theme="TeamFortress"
  )
  qgraph(xi_detrend
         , filetype='pdf'
         , filename=file.path(out_folder
                              , glue('mlvar_{type}'
                              )
         )
  )
}

## Saving network results as CSV

s = network_detrend%>%summary()

for (type in c('contemporaneous'
               , 'temporal'
               , 'between'
)
)
{
  write.csv(s[[type]],
            file.path(out_folder,
                      glue('mlvar_{type}.csv')
                      
            ),
            row.names=F,
  )
}

# Assumption checking and mlVAR analysis for data with synthetic constructs
# modeling data gaps plus last-value imputation.
BRC_data <- read.csv('Originals+Synthetics.csv'
                     , header=T
                     , sep=","
                     ) 

BRC_data$encounter_date <- as.Date(BRC_data$encounter_date,
)
BRC_data %>% head()

### Factor Analysis
factor_data <- (BRC_data 
                %>% dplyr::select(where(is.numeric) 
                                  & !participant_id
                                  & !day_n
                )
)

fa.parallel(factor_data, 
            fm="ml",
)
efaResPos <- (factor_data 
              %>% fa(nfactors=4,
                     rotate="oblimin",
                     fm="ml",
              )
)

fa.diagram(efaResPos)


# Stationarity tests
## Visual inspection

valid_participants=(BRC_data 
                    %>% group_by(participant_id)
                    %>% count()
                    %>% arrange(desc(n))
)$participant_id
constructs=colnames(factor_data)
colours = colorRampPalette(brewer.pal(8
                                      , 'Dark2'
)
)(length(constructs))

for (i in seq_along(constructs))
{
  construct = constructs[i]
  y_label=sprintf("%s rating", construct)
  plot=(BRC_data 
        %>% subset(participant_id
                   %in% valid_participants
        )
        %>% ggplot(aes(x=day_n,
                       y=.data[[construct]]
        )
        )
        + geom_line(color=colours[i])
        + scale_x_continuous("Date time")
        + scale_y_continuous(construct)
        + facet_wrap(~participant_id
                     , ncol=6
                     , shrink=FALSE
        )
        + ylab(y_label)
        + theme(panel.grid.minor=element_blank() 
                , panel.grid.major=element_blank()
                , strip.text=element_text(size=5
                                          , margin=margin()
                )
                , axis.text.x=element_text(size=6)
                , axis.text.y=element_text(size=6)
        )
  )
  print(plot)
}

## Kwiatkowski-Phillips-Schmidt-Shin test is used for testing a null hypothesis that an observable time series is stationary around a deterministic trend (i.e. trend-stationary) against the alternative of a unit root
# Indicated by the parameter null="Trend" of tseries::kpss.test.
# H0: trend-stationary, H1: non-stationary

constructs = colnames(factor_data)
participants = unique(BRC_data$participant_id)
m = length(participants)

for (construct in constructs)
{
  kpps_p_construct <- c()
  for(p in participants)
  {
    BRC_construct_participant = (BRC_data
                                 [[construct]]
                                 [BRC_data$participant_id==p]
    )
    randomness <- tseries::kpss.test(na.exclude(BRC_construct_participant)
                                     , lshort=TRUE
                                     , null="Trend"
    )
    kpps_p_construct <- c(kpps_p_construct,
                          randomness$p.value
    )
  }
  kpps_p_construct  <- cbind(participants
                             , kpps_p_construct
  )
  cat(paste('KPSS test Bonferroni(alpha<0.05,'
            , glue('m={m}) for {construct},')
            , 'is it trend-stationary?:'
            , '\n'
            , glue('{all(kpps_p_construct[,2] < 0.05/m)}')
            , '\n'
  )
  )
}

# Normality
## Visual inspection
### Histogram vs Normal distribution
constructs=colnames(factor_data)
colours = colorRampPalette(brewer.pal(8
                                      , 'Dark2'
)
)(length(constructs))
png("HistvNorm_ImputedData.png"
    , width=800
    , height=400
)
par(mfrow=c(2, 4))
for (i in seq_along(constructs))
{
  construct <- constructs[[i]]
  g <- density(BRC_data[[construct]])
  
  plot(g
       , ylab='density'
       , xlab="values"
       , main=construct
  )
  polygon(g, col=colours[[i]])
  xfit <- seq(min(g$x)
              , max(g$x)
              , length=50
  ) 
  yfit <- dnorm(xfit
                , mean=mean(g$x)
                , sd=sd(g$x)
  )
  lines(xfit
        , yfit
        , col='black'
        , lwd=2
  )
}
dev.off()

### QQ Plot
png("QQPlot_ImputedData.png"
    , width=800
    , height=400
)
par(mfrow=c(2, 4))
for (i in seq_along(constructs))
{
  construct <- constructs[[i]]
  qqnorm(BRC_data[[construct]]
         , main=construct
         , col=colours[[i]]
  )
  qqline(BRC_data[[construct]])
}
dev.off()


## Kolmogorov-Smirnov: Test to assess whether the data is sampled from a normal distribution, namely we assume normality.
### H0: normality, H1: non-normality
### H0 is true if p-value > 0.05

for (construct in constructs)
{
  ks_construct = ks.test(BRC_data
                         [[construct]]
                         , rnorm(length(BRC_data
                                        [[construct]]
                         )
                         )
  )$p.value
  cat(paste(glue('KS test for {construct},')
            , 'is it normally distributed?:'
            , '\n'
            , glue('{all(ks_construct > 0.05)}')
            , '\n'
  )
  )
}

# ## Shapiro-Wilk test: Test to assess whether the data is sampled from a normal distribution, namely we assume normality. This test is more robust than Kolmogorov-Smirnov but it is limited to samples <= 5,000.To overcome this last limitation 5,000 samples (without replacement) are performed.
# ### H0: normality, H1: non-normality
# ### H0 is true if p-value > 0.05
# 
# for (construct in constructs)
# {
#   sw_construct = shapiro.test(BRC_data[[construct]]
#                               %>% sample(5000)
#                               )
#   cat(paste(glue('SW test for {construct}, is it normally distributed?:')
#             , '\n'
#             , glue('{all(sw_construct$p.value > 0.05)}')
#             , '\n'
#   )
#   )
# }

## Anderson-Darling normality: This test to assess whether the data is sampled from a normal distribution, namely we assume normality. This test is not as robust as Shapiro-Wilk but it is not limited to samples <= 5,000
### H0: normality, H1: non-normality
### H0 is true if p-value > 0.05

for (construct in constructs)
{
  # print(construct)
  ad_construct = ad.test(BRC_data[[construct]])
  cat(paste(glue('AD test for {construct}, '
                 , 'is it normally distributed?:'
  )
  , '\n'
  , glue('{all(ad_construct$p.value > 0.05)}')
  , '\n'
  )
  )
}

# Main mlVAR analysis
## Detrending: Despite stationarity results seeming ok (at least for 0.05 alpha) we were concerned as participants used the app for such a long time period, wrote a loop which detrends each participant individually and takes residuals.

consec_detrend=(BRC_data
                %>% arrange(dataset,
                            participant_id,
                            day_n
                )
                %>% dplyr::select(dataset,
                                  participant_id,
                                  day_n,
                                  everything()
                )
)

# Sorted by volume
participants=(BRC_data 
              %>% group_by(participant_id)
              %>% count()
              %>% arrange(desc(n)) # _n_ is the counts column
)$participant_id

constructs = colnames(BRC_data)[-1:-4]

# Based on S Allen et al
for (construct in constructs)
{
  for(participant in valid_participants)
  {
    detrended_data <- (lm(na.exclude(consec_detrend[[construct]]
                                     ~ consec_detrend$day_n,
                                     data=consec_detrend
                                     [consec_detrend$participant_id==participant]
    )
    )
    %>% resid()
    )
  }
  consec_detrend[[construct]][!is.na(consec_detrend[[construct]])] <- detrended_data
}

## mlVAR networks analysis

constructs = colnames(factor_data)

network_detrend <- mlVAR(consec_detrend
                         , vars=constructs
                         , idvar="participant_id"
                         , estimator="default"
                         , contemporaneous="correlated"
                         , temporal="correlated"
                         , beepvar="day_n"
                         , lags=1
                         )

## Plotting network results

for (i in mapply(list,
                 c(0.00625, 0.003125, 0.00625),
                 c('temporal'
                   , 'between'
                   , 'contemporaneous'
                 ),
                 SIMPLIFY=F)
     )
{
  a = i[[1]]
  type = i[[2]]
  xi_detrend <- plot(network_detrend
                     , type
                     , layout="circle"
                     , rule="and"
                     , legend=FALSE
                     , edge.labels=T
                     , alpha=a
                     , cut=0.1
                     , theme="TeamFortress"
                     )
  qgraph(xi_detrend
         , filetype='pdf'
         , filename=file.path(out_folder
                              , glue('mlvar_{type}'
                                     )
                              )
         )
}

## Saving network results as CSV

s = network_detrend%>%summary()

for (type in c('contemporaneous'
               , 'temporal'
               , 'between'
               )
     )
  {
  write.csv(s[[type]],
            file.path(out_folder,
                      glue('mlvar_{type}.csv')
                      
                                  ),
            row.names=F,
            )
  }