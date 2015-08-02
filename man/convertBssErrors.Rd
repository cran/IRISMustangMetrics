\name{convertBssErrors}

\alias{convertBssErrors}


\title{Generate Human Readable MUSTANG Errors}

\description{
The MUSTANG database is in charge of storing the results of metrics calculations
and is accessed through a webservice API.
The\code{convertBssErrors} function extracts pertinent error information from the HTML returned
by MUSTANG on error conditions.
}
\usage{
convertBssErrors(err_msg)
}

\arguments{
  \item{err_msg}{error text received from the MUSTANG}
}

%% \details{ }

\value{
Aext string with the \code{root cause} extracted from 
the MUSTANG HTML Java error dump.
}

%% \references{ }

\author{
Jonathan Callahan \email{jonathan.s.callahan@gmail.com}
}

%% \note{ }

%% ~Make other sections like Warning with \section{Warning}{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{SingleValueMetric-class}},
\code{\link{metricList2Xml}},
\code{\link{getMetricsXml}},
\code{\link{getBssMetricList}},
}

%% \examples{ }

