\name{shrinkIt}
\alias{shrinkIt}
\title{Perform shrinkage of subject-level values towards group mean}
\usage{
shrinkIt(observation, obs1 = NULL, obs2 = NULL, noiseVar = "common")
}
\arguments{
  \item{observation}{The observation matrix (see
  description) on which shrinkage will be performed}

  \item{obs1}{Observation matrix (see description) of the
  first observation for each subject}

  \item{obs1}{Observation matrix (see description) of the
  second observation for each subject}

  \item{noiseVar}{Either a character string equal to
  "common", "individual" or "scaled", in which case the
  noise variance will be computed, or a numeric,
  non-negative scalar to be used as the global noise
  variance.}
}
\value{
A list containing the shrunken observation matrix, shrinage
parameter (lambda) values, noise variance values, and
signal variance values.  The shrinkage parameter values
range between 0 and 1, with 1 representing complete
shrinkage towards the group mean and 0 representing no
shrinkage towards the group mean.
}
\description{
This function performs shrinkage of subject values towards
the group mean. A common (group) noise variance, individual
noise variance, or individually-scaled noise variance may
be computed.  The noise variance is computed for each
variable 1 to m. Alternatively, a pre-specified global
noise variance may be used.
}
\details{
Each observation matrix is m-by-n, where m is the number of
variables observed for each subject, and n is the number of
subjects.  The variables should match across subjects.  For
example, in an fMRI context, m might be the number of
elements in the upper triangle of each subject's V-by-V
correlation matrix.  If noise variance is to be computed,
obs1 and obs2 are used to compute the noise variance.  If a
global noise variance value is provided, obs1 and obs2 will
not be used.
}
\examples{
\dontrun{

}
}

