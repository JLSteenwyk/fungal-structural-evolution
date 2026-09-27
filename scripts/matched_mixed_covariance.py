#!/usr/bin/env python3
"""Nested background/family covariance plus a low-rank species term.

Numerical likelihood evaluator only. Optimizer, model adequacy and uncertainty
calibration are separate requirements. Components are variance values, not SDs.
"""
import numpy as np
from scipy.linalg import cho_factor, cho_solve


class MatchedCovariance:
    def __init__(self, background, family, species_factor,
                 residual=1., background_variance=0., family_variance=0., species_variance=0.):
        variances = np.array([residual, background_variance, family_variance, species_variance], dtype=float)
        if not np.isfinite(variances).all() or residual <= 0 or np.any(variances[1:] < 0):
            raise ValueError('Residual variance must be positive; other variances must be nonnegative and finite')
        background, family = np.asarray(background), np.asarray(family)
        if background.ndim != 1 or family.shape != background.shape or len(family) == 0:
            raise ValueError('Expected nonempty equally sized background and family labels')
        self.n = len(family)
        _, self.bg = np.unique(background, return_inverse=True)
        _, self.fam = np.unique(family, return_inverse=True)
        self.ng = int(self.bg.max()) + 1
        self.nf = int(self.fam.max()) + 1
        membership = np.full(self.ng, -1, dtype=int)
        for bg, fam in zip(self.bg, self.fam):
            if membership[bg] not in [-1, fam]:
                raise ValueError('A background occurs in multiple family components; nested solver is inapplicable')
            membership[bg] = fam
        self.residual = float(residual)
        self.bvar = float(background_variance)
        counts = np.bincount(self.bg)
        denominator = residual + background_variance * counts
        self.background_coefficient = background_variance / (residual * denominator)
        self.family_vector = 1. / denominator[self.bg]
        family_information = np.bincount(self.fam, weights=self.family_vector)
        self.family_coefficient = family_variance / (1. + family_variance * family_information)
        self.logdet = float((self.n-self.ng)*np.log(residual) + np.log(denominator).sum()
                            + np.log1p(family_variance*family_information).sum())
        self.factor = np.asarray(species_factor, dtype=float)
        if self.factor.ndim != 2 or self.factor.shape[0] != self.n or not np.isfinite(self.factor).all():
            raise ValueError('Species factor must be a finite n-by-r matrix')
        self.pvar = float(species_variance)
        self.phylogenetic_cholesky = None
        if species_variance > 0 and self.factor.shape[1]:
            self.base_inverse_factor = self._base_solve(self.factor)
            core = np.eye(self.factor.shape[1]) + species_variance * (self.factor.T @ self.base_inverse_factor)
            core = (core + core.T) / 2
            self.phylogenetic_cholesky = cho_factor(core, lower=True, check_finite=True)
            self.logdet += float(2 * np.log(np.diag(self.phylogenetic_cholesky[0])).sum())

    @staticmethod
    def _group_sum(values, labels, groups):
        output = np.zeros((groups, values.shape[1]))
        np.add.at(output, labels, values)
        return output

    def _base_solve(self, values):
        grouped = self._group_sum(values, self.bg, self.ng)
        solved = values / self.residual - self.background_coefficient[self.bg, None] * grouped[self.bg]
        family_sum = self._group_sum(solved, self.fam, self.nf)
        return solved - (self.family_vector * self.family_coefficient[self.fam])[:, None] * family_sum[self.fam]

    def solve(self, values):
        values = np.asarray(values, dtype=float)
        vector = values.ndim == 1
        if vector:
            values = values[:, None]
        if values.ndim != 2 or values.shape[0] != self.n or not np.isfinite(values).all():
            raise ValueError('Right-hand side must be finite with n rows')
        result = self._base_solve(values)
        if self.phylogenetic_cholesky is not None:
            correction = cho_solve(self.phylogenetic_cholesky, self.factor.T @ result)
            result -= self.pvar * (self.base_inverse_factor @ correction)
        return result[:, 0] if vector else result


def profiled_reml(covariance, design, response):
    """Profile an overall variance multiplier at fixed covariance component ratios.

    Returns a likelihood criterion and conditional coefficients, not a completed
    fit or calibrated confidence interval. Design columns must be full rank.
    """
    x, y = np.asarray(design, dtype=float), np.asarray(response, dtype=float)
    if x.ndim != 2 or y.shape != (covariance.n,) or x.shape[0] != covariance.n:
        raise ValueError('Incompatible design or response dimensions')
    n,p = x.shape
    if p == 0 or n <= p or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Need a finite nonempty design with positive residual degrees of freedom')
    if np.linalg.matrix_rank(x) != p:
        raise ValueError('Fixed-effect design is rank deficient')
    inverse_x = covariance.solve(x)
    information = x.T @ inverse_x
    information = (information + information.T)/2
    chol = cho_factor(information, lower=True)
    beta = cho_solve(chol, x.T @ covariance.solve(y))
    residual = y - x @ beta
    quadratic = float(residual @ covariance.solve(residual))
    if not np.isfinite(quadratic) or quadratic <= 0:
        raise ValueError('Nonpositive profiled residual quadratic form')
    scale = quadratic/(n-p)
    logdet_information = float(2*np.log(np.diag(chol[0])).sum())
    objective = .5*(covariance.logdet + logdet_information + (n-p)*(1.+np.log(2*np.pi*scale)))
    return dict(negative_profiled_reml=float(objective), beta=beta,
                profiled_scale=float(scale), residual_quadratic=quadratic,
                conditional_beta_covariance=scale*cho_solve(chol,np.eye(p)),
                residual_degrees_of_freedom=n-p)
