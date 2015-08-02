\name{dailyDCOffsetMetric}

\alias{dailyDCOffsetMetric}


\title{DC Offset Detection}

\description{
The dailyDCOffsetMetric() function identifies days with a jump in the signal mean.
}

\usage{
dailyDCOffsetMetric(df, offsetDays=5,
                    outlierWindow=7,
                    outlierThreshold=6.0,
                    outlierSelectivity=0.1)
}

\arguments{
  \item{df}{a dataframe containing \code{sample_mean} values obtained with getSingleValueMetrics()}
  \item{offsetDays}{number of days used in calculating weighting factors}
  \item{outlierWindow}{window size passed to findOutliers() function in the \pkg{seismic} package}
  \item{outlierThreshold}{detection threshold passed to findOutliers() function in the \pkg{seismic} package}
  \item{outlierSelectivity}{selectivity passed to findOutliers() function in the \pkg{seismic} package}
}

\details{
This algorithm calculates lagged differences of the daily mean timeseries over a window of \code{offsetDays} days.
Shifts in the mean that are persistent and  larger than the typical standard deviation of daily means will
generate higher metric values.

%% max preformatted length ---------------------------------------------------
Details of the algorithm are as follows
\preformatted{
# data0 = download 12+ months of daily means (in the 'df' dataframe)
# data1 = remove outliers using MAD outlier detection with the 'outlier' arguments specified
# data2 = replace outliers with rolling median values using a 7 day window
# weights = calculate absolute lagged differences with 1-N day lags (big jumps have large values)
# metric0 = multiply the lagged differences together and take the N'th root
# stddev0 = calculate the rolling standard deviation of data2 with a N-day window
# METRIC = divide metric0 by the median value of stddev0
}
}

\value{
A list is returned with a \code{MultipleTimeValueMetric} object for each day in the incoming dataframe..
}

%% \references{ }

\author{
Jonathan Callahan \email{jonathan.s.callahan@gmail.com}
}

\note{
After initial testing on stations in the IU network, a metric value > 10 appears to be indicative of a DC offset shift.
}

\seealso{
\code{\link{getSingleValueMetrics}}
}

%% \examples{}

\keyword{metrics}
