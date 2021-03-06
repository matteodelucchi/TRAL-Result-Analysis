% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_AAratio_comparison.R
\name{plot_AAratio_comparison}
\alias{plot_AAratio_comparison}
\title{Compare the AA frequency in two groups.}
\usage{
plot_AAratio_comparison(df = AAfreq, plot_title = NULL, ylim = 1)
}
\arguments{
\item{df}{data.frame. AA frequency count obtained through \code{\link{AAfreq_in_TR}} with an additional column "aa_ratio_overall".}

\item{plot_title}{str. Plot title.}

\item{ylim}{int. Scale the y-axis maximum value.}
}
\value{
ggplot object. The values from 'aa_ratio' are displayed in green, the values from 'aa_ratio_overall' in red.
}
\description{
Compare the AA frequency in two groups.
}
