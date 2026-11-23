# AI assistance statement

This statement is provided under the Bioconductor policy on AI-generated and
third-party code
(<https://contributions.bioconductor.org/ai-policy-third-party.html>).

## What was assisted

The code, documentation, unit tests and vignette of `BreedingExperiment` were
written with assistance from a large language model during the package's
initial development, up to and including its submission to Bioconductor.

* **Tool:** Claude Code, the command-line coding assistant from Anthropic.
* **Model:** Anthropic's Claude, principally Claude Opus 5.

## What was not assisted

No part of the package incorporates code copied from another software
project. Everything in it was written for this package and is redistributable
under the MIT licence declared in `DESCRIPTION`.

The statistical methods are not novel and were not devised by a model. Each
is an established, published method, implemented from its original
description and cited in the documentation of the function that implements
it:

* the numerator relationship matrix of Henderson (1976),
  <doi:10.2307/2529339>;
* the genomic relationship matrix of VanRaden (2008),
  <doi:10.3168/jds.2007-0980>;
* the dominance parameterisation of Vitezica, Varona and Legarra (2013),
  <doi:10.1534/genetics.113.155176>;
* the single-step matrix of Legarra, Aguilar and Misztal (2009),
  <doi:10.3168/jds.2009-2061>, with the adjustments of Christensen and Lund
  (2010), <doi:10.1186/1297-9686-42-2>, as applied by Aguilar and others
  (2010), <doi:10.3168/jds.2009-2730>.

## Verification

Assistance was directed and reviewed by the package author, who is
responsible for the correctness, reliability and continued maintenance of all
code in the package.

The implementations are checked against independent references rather than
taken on trust. The test suite in `tests/testthat` verifies, among other
things, that the pedigree relationship matrix reproduces textbook values
(parent-offspring 0.5, grandparent-grandchild 0.25, and an inbreeding
coefficient of 0.25 for the offspring of full sibs), that the dominance
coding reproduces the three genotype values published by Vitezica and others
exactly, that blending the single-step matrix fully towards the pedigree
returns the pedigree submatrix, and that the single-step matrix is positive
definite. The vignette additionally validates the package end to end by
predicting the breeding values of a withheld generation and comparing them
with the simulated truth.

## Maintainer

Muhammad Farooqi, Independent Researcher
(ORCID [0000-0003-4918-9791](https://orcid.org/0000-0003-4918-9791)).