\name{PSDMetric}

\alias{PSDMetric}

\title{Power Spectral Density of a signal}

\description{
The PSDMetric() function performs spectral analysis on a seismic signal
and returns 'PSD' metrics with discretized spectral components as well as other
metrics based on PSD's
}

\usage{
PSDMetric(st,
          expLoPeriod=4/(st@traces[[1]]@stats@sampling_rate),
          expHiPeriod=100,
          linLoPeriod=16/(st@traces[[1]]@stats@sampling_rate),
          linHiPeriod=50)
}

\arguments{
  \item{st}{a \code{Stream} object containing a seismic signal}
  \item{expLoPeriod}{the low end of the period band use for calculating the exponential dead channel detector}
  \item{expHiPeriod}{the high end of the period band use for calculating the exponential dead channel detector}
  \item{linLoPeriod}{the low end of the period band use for calculating the linear dead channel detector}
  \item{linHiPeriod}{the high end of the period band use for calculating the linear dead channel detector}
}

\details{
This metric calculates average power spectra for a seismic signal as described in the McNamara paper.
See the \code{McNamaraPSD} method of \code{Stream} objects in the \pkg{seismic} package for details.

\emph{Uncorrected} spectral density values are returned in \code{spectrumMetricList} in units of dB.

Other metrics calculated from the PSDs are returned in \code{svMetricList}.  These metrics are:

\describe{

  \item{\bold{pct_above_nhnm} -- "percent above New High Noise Model"}{Percentage of PSD values that
  are above the New High Noise Model for their frequency. This value is calculated over the entire time period.}
  

  \item{\bold{pct_below_nlnm} -- "percent below New Low Noise Model"}{Percentage of PSD values that are below the
  New Low Noise Model for their frequency. This value is calculated over the entire time period.}

  \item{\bold{dead_channel_exp} -- "dead channel metric - exponential fit"}{A "dead channel" metric is calculated
  from the mean of all the PSDs generated. (Typically 47 for a 24 hour period.) Values of the PSD mean line over
  the band (expLoPeriod:expHiPeriod) are fit to an exponential.  The \code{dead_channel_exp} metric is the standard
  deviation of the fit residuals.  Lower numbers indicate a better fit and a higher likelihood that the mean PSD is
  exponential -- an indication of a "dead channel".}

  \item{\bold{dead_channel_lin} -- "dead channel metric - linear fit"}{ A "dead channel" metric is calculated from
  the mean of all the PSDs generated. (Typically 47 for a 24 hour period.) Values of the PSD mean line over the band
  (linLoPeriod:linHiPeriod) are fit to a line.  The \code{dead_channel_lin} metric is the standard deviation of the
  fit residuals.  Lower numbers indicate a better fit and a higher likelihood that the mean PSD is linear -- an
  indication of a "dead channel".}

}

}

\value{
  A list of lists is returned containing:
  \itemize{ 
    \item{\code{spectrumMetricList} = list of \code{SpectrumMetric} objects}
    \item{\code{svMetricList} = list of \code{SingleValueMetric} objects:}
    \itemize{
      \item{\code{pct_above_nhnm}}
      \item{\code{pct_below_nlnm}}
      \item{\code{dead_channel_exp}}
      \item{\code{dead_channel_lin}}
    }
  }
}

\references{
\href{http://pubs.usgs.gov/of/2005/1438/pdf/OFR-1438.pdf}{Seismic Noise Analysis System Using Power Spectral Density Probability Density Functions} (McNamara and Boaz 2005)
}

\author{
Jonathan Callahan \email{jonathan.s.callahan@gmail.com}
}

%% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpectrumMetric}}
\code{\link{SingleValueMetric}}
}

\examples{
  \dontrun{
# Open a connection to IRIS DMC webservices
iris <- new("IrisClient")

# NOTE:  The following trace has 1.728 million points.
# NOTE:  Downloading and calculating PSD may take a few seconds.
starttime <- as.POSIXct("2010-02-27",tz="GMT")
endtime <- as.POSIXct("2010-02-28",tz="GMT")
  
# Get the waveform
st <- getDataselect(iris,"IU","ANMO","00","BHZ",starttime,endtime)

# Calculate the PSD metric and show the SingleValueMetric results
listOfLists <- PSDMetric(st)
svMetricList <- listOfLists[['svMetricList']]

dummy <- lapply(svMetricList, show)
  }
}

\keyword{metrics}
