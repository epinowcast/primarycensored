# primarycensored: Primary Event Censored Distributions

Provides functions for working with primary event censored distributions
and 'Stan' implementations for use in Bayesian modeling. Primary event
censored distributions are useful for modeling delayed reporting
scenarios in epidemiology and other fields (Charniga et al. (2024)
[doi:10.48550/arXiv.2405.08841](https://doi.org/10.48550/arXiv.2405.08841)
). It also provides support for arbitrary delay distributions, a range
of common primary distributions, and allows for truncation and secondary
event censoring to be accounted for (Park et al. (2024)
[doi:10.1101/2024.01.12.24301247](https://doi.org/10.1101/2024.01.12.24301247)
). A subset of common distributions also have analytical solutions
implemented, allowing for faster computation. In addition, it provides
multiple methods for fitting primary event censored distributions to
data via optional dependencies.

## See also

Useful links:

- <https://primarycensored.epinowcast.org>

- <https://github.com/epinowcast/primarycensored>

- Report bugs at <https://github.com/epinowcast/primarycensored/issues>

## Author

**Maintainer**: Sam Abbott <contact@samabbott.co.uk>
([ORCID](https://orcid.org/0000-0001-8057-8037)) \[copyright holder\]

Authors:

- Sam Brand <usi1@cdc.gov>
  ([ORCID](https://orcid.org/0000-0003-0645-5367))

- James Mba Azam <james.azam@lshtm.ac.uk>
  ([ORCID](https://orcid.org/0000-0001-5782-7330))

- Carl Pearson <carl.ab.pearson@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-0701-7860))

- Sebastian Funk <sebastian.funk@lshtm.ac.uk>
  ([ORCID](https://orcid.org/0000-0002-2842-3406))

- Kelly Charniga <kelly.charniga@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-7648-7041))

Other contributors:

- Adam Howes <adamthowes@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-2386-4031)) \[contributor\]
