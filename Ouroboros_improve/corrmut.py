#!/usr/bin/python
"""
Module containing functions for analysis of correlated mutations combined with
expectation-maximization.

@author: Miguel Correa Marrero
"""

import os

import numpy as np
from tqdm import tqdm

from copy import copy
import warnings

from sklearn.linear_model import SGDClassifier

from dummyestimator import DummyEstimator
from globalvars import ALPHA_RANGE
import output
from contacts import compute_couplings, get_interacting, normalize_contact_mtx
from helpers import round_labels


##################################
# Hidden variables / convergence #
##################################

def get_random_labels(no_seqs, int_frac, mode):
    """
    Function to create random labels for interaction / non-interaction

    Arguments
    ---------
    no_seqs: int, number of proteins in the alignments
    int_frac: float, assumed fraction of interacting proteins
    mode: str, whether we are performing hard or soft EM

    Returns
    ---------
    labels: list, values of the hidden variables
    """

    if mode == 'hard':
        labels = np.random.choice([0, 1], no_seqs, p=[1 - int_frac, int_frac])
    elif mode == 'soft':
        labels = np.random.choice([0.1, 0.9], no_seqs, p=[
                                  1 - int_frac, int_frac])
    return labels


def update_labels(alt_llhs, null_llhs, int_frac, mode):
    """
    Function to update the hidden variables after each iteration of the EM
    loop.

    Arguments
    ---------
    alt_llhs:   array-like, contains log-likelihoods of the alternative model
    null_llhs:  array-like, contains log-likelihoods of the null model
    int_frac:   float, prior probability of fraction of interacting proteins
    mode:       string, type of expectation-maximization

    Returns
    --------
    labels:     array-like, values of the hidden variables
    """
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        # Try-except block to catch floating point underflow
        try:
            numerator = np.exp(alt_llhs) * int_frac
            denominator = (np.exp(alt_llhs) * int_frac) + \
                (np.exp(null_llhs) * (1 - int_frac))
            labels = numerator / denominator
        except Warning:  # Increase floating precision to 80 bits
            try:
                numerator = np.exp(alt_llhs, dtype='longdouble') * int_frac
                denominator = (np.exp(alt_llhs, dtype='longdouble') * int_frac) + \
                    (np.exp(null_llhs, dtype='longdouble') * (1 - int_frac))
                labels = numerator / denominator
            except ZeroDivisionError:
                # This should not happen unless:
                #    - int_frac >= 1 or int_frac <= 0 (will raise exception at the start)
                #    - floating point underflow after trying with 80 bits precision
                #      (possibly due to extremely long proteins)
                raise Exception(
                    "Unexpected error: denominator in update equation is 0. Contact developer.")

    if mode == 'hard':
        labels = round_labels(labels)
    return labels


def has_converged(labels, pre_labels, mode, tol=0.005):
    """
    Function to determine whether the expectation-maximization loop has reached
    convergence

    Arguments
    ---------
    labels:     array-like, contains the vector of hidden variables at the
                current iteration
    pre_labels: array-like, contains the vector of hidden variables at the
                previous iteration
    mode:       str, mode in which the algorithm is being run
    tol:        float, tolerance, the difference threshold to be considered
                when using the algorithm in soft mode

    Returns
    -------
    converged: bool, whether the algorithm has converged or not
    """

    converged = False
    sum_abs_diffs = np.sum(np.absolute(np.subtract(labels, pre_labels)))
    if mode == 'hard':
        # Convergence reached when there is no difference with the previous
        # step
        if sum_abs_diffs == 0:
            converged = True

    elif mode == 'soft':
        # Convergence reached when the difference between the two vectors is
        # smaller than the tolerance
        if sum_abs_diffs / (len(labels)) < tol:
            converged = True

    return converged


#############################
# Alternative model fitting #
#############################


def calc_alt_llhs(num_mtx_a, bin_mtx_b, models_a, num_mtx_b, bin_mtx_a,
                  models_b, out_dir, iters):
    """
    Given two multiple sequence alignments in numeric matrix and binary matrix
    form, calculate sequence pair alternative log-likelihoods.

    Arguments
    ---------
    num_mtx_a, num_mtx_b: array-like. Contains a multiple sequence alignment as
                          a numeric matrix with all proteins under
                          consideration, interacting and non-interacting.
                          Each matrix is of size
                          (number of sequences x number of MSA columns)
    bin_mtx_a, bin_mtx_b: array-like. Contains a multiple sequence alignment as
                          a binary matrix with all proteins under
                          consideration, interacting and non-interacting.
                          Each matrix is of size
                          (number of sequences x
                          (number of columns x number of allowed amino acids))
    models_a, models_b:   list, contains one fitted SGDClassifer object per MSA
                          column

    Returns
    -------
    alt_llhs:             list, contains alternative sequence pair
                          log-likelihoods
    """
    alt_mtx1 = get_alt_model(num_mtx_a, bin_mtx_b, models_a)
    alt_mtx2 = get_alt_model(num_mtx_b, bin_mtx_a, models_b)
    concat_mtx = np.concatenate((alt_mtx1, alt_mtx2), axis=1)
    np.savetxt(os.path.join(out_dir, "".join(
        ["alt_llhs_mtx_", str(iters), ".csv"])), concat_mtx, delimiter=',')
    alt_llhs = np.sum(concat_mtx, axis=1)
    return alt_llhs


def get_alt_model(num_mtx, bin_mtx, models, pc=np.log(1 / 210)):
    """
    Given a multiple sequence alignment in numeric matrix form, the other
    alignment in binary matrix form, and the fitted models, this function
    will find an appropiate value of lambda for each model and return the
    fitted probabilities according to the logistic models.

    Arguments
    ---------
    num_mtx: array-like. Contains a multiple sequence alignment as a numeric
             matrix with all proteins under consideration,
             interacting and non-interacting.
    bin_mtx: array-like. Contains a multiple sequence alignment as a binary
             matrix with all proteins under consideration,
             interacting and non-interacting.
    models:  list, contains one fitted model object per MSA column
    pc:      float, pseudocount for when a residue was not present in the
             training data, as its probability cannot be estimated properly

    Returns
    ---------
    alt_mtx: array-like. Contains the values of the log-probability of the data
             according to the logistic models, element-wise.
    """
    alt_mtx = np.zeros_like(num_mtx,dtype='float64').T

    # Iterate over columns and their corresponding models
    for i, col in enumerate(num_mtx.T):
        cur_model = models[i]

        # Get model predictions
        log_probs = cur_model.predict_log_proba(bin_mtx)

        for j, res in enumerate(col):
            if res in cur_model.classes_:
                # Get index of residue type in self.classes_
                res_idx = np.asscalar(np.where(cur_model.classes_ == res)[0])
                log_prob = log_probs[:, res_idx][j]
                # Check whether it has 0 probability
                if log_prob != float('-inf'):
                    alt_mtx[i][j] = log_prob  # Use predicted probability
                else:  # If it is 0, use a pseudocount
                    alt_mtx[i][j] = pc
            # If the residue was not in the training data, use pseudocount
            # (possible in hard EM)
            else:
                alt_mtx[i][j] = pc

    # Return alternative model matrix
    return alt_mtx.T


def select_interacting(num_mtx, bin_mtx, labels):
    """
    Auxiliary function for fit_msa_mdels.
    Used for fitting the models in hard EM; selects observations with a hidden
    variable value of 1.
    """
    if labels is None:
        # This is the case when initializing the models
        return num_mtx, bin_mtx, labels
    else:
        # This is the case inside the EM loop
        labels = np.asarray(labels)
        idxs = np.where(np.asarray(labels) == 1)[0]

        int_num = np.take(num_mtx, idxs, axis=0)
        int_bin = np.take(bin_mtx, idxs, axis=0)
        weights = np.take(labels, idxs)

        return int_num, int_bin, weights


def fit_msa_models(num_mtx, bin_mtx, mode, seqs_weight, fixed_alphas=None, n_jobs=2,
                   sample_weights=None, l1_ratio=0.99, dfmax=100,
                   random_state=42, sgd_tol=1e-3):
    """
    Given two MSAs, one in numeric matrix format and another in binary matrix
    format, fit logistic regressions for each column in the numeric matrix
    using the binary matrix as predictors.

    Arguments
    ---------
    num_mtx:            array-like, MSA in numeric matrix form
    bin_mtx:            arrray-like, MSA in binary matrix form
    mode:               string, whether we are performing 'soft' or 'hard' EM
    fixed_alphas:       list, values of alpha to use in model fitting
    n_jobs:             int, number of CPUs to use in model fitting
    seqs_weight:        list, weight of each observation when they are in one cluster
    sample_weights:     list, weight for each observation
    l1_ratio:           float, elastic net mixing parameter
    dfmax:              int, maximum number of degrees of freedom allowed in
                        the models
    random_state:       int, random state for stochastic gradient descent
    sgd_tol:            float, tolerance for stochastic gradient descent

    Returns
    -------
    models:             list of fitted SGDClassifier objects
    alpha_per_col:      list, selected values of alpha; only returned if no
                        value was passed to fixed_alphas
    """
    models = []
    n_obs = num_mtx.shape[0]
    alpha_per_col = []

    # Select cases with a hidden variable value of 1 in hard EM
    if mode == 'hard':
        num_mtx, bin_mtx,\
            sample_weights = select_interacting(num_mtx, bin_mtx,
                                                sample_weights)

    # Fit models for each column of the MSA
    for idx, col in enumerate(tqdm(num_mtx.T)):
        if len(np.unique(col)) > 1:  # Column contains more than one class
            col_models = []
            col_dfs = []
            col_bics = []
            # Initialization: if no predefined values of the regularization
            # strenght are given, train models on a range of them and select one
            if fixed_alphas is None:
                for alpha in ALPHA_RANGE:
                    clf = SGDClassifier(loss='log', penalty='elasticnet',
                                        alpha=alpha, l1_ratio=l1_ratio,
                                        n_jobs=n_jobs, max_iter=100,
                                        random_state=random_state, tol=sgd_tol)
                    # now the sample weights is none
                    clf.fit(bin_mtx, col, sample_weight=sample_weights)
                    
                    # Discard models with a number of degrees of freedom above
                    # a certain threshold
                    dfs = calc_degrees_freedom(clf)
                    if dfs <= dfmax:
                        col_models.append(clf)
                        col_dfs.append(dfs)
                    else: # Stop once models are overtly complex
                        break

                # Select regularization strenght by choosing the model with the
                # minimum Bayesian Information Criterion
                for j, model in enumerate(col_models):
                    posterior_logprobs = get_posterior_logprobs(col,
                                                                bin_mtx, model)
                    bic = calc_bic(posterior_logprobs, col_dfs[j], n_obs)
                    col_bics.append(bic)
                best_idx = col_bics.index(min(col_bics))
                models.append(col_models[best_idx])
                alpha_per_col.append(ALPHA_RANGE[best_idx])

            else:
                # EM iterations after initialization: if predefined values of
                # the regularization strength are given, use those to fit models
                clf = SGDClassifier(loss='log', penalty='elasticnet',
                                    alpha=fixed_alphas[idx], l1_ratio=l1_ratio,
                                    n_jobs=n_jobs, max_iter=1000,
                                    random_state=random_state, tol=sgd_tol)
                clf.fit(bin_mtx, col, sample_weight=np.multiply(sample_weights, seqs_weight))
                models.append(clf)
        else:  # Column contains only one class; use a dummy model
            # Can happen in hard EM
            clf = DummyEstimator(prob=0.99 - (1 / 210))
            clf.fit(bin_mtx, col)
            models.append(clf)
            if fixed_alphas is None:
                # Commonly selected value, strong regularization
                alpha_per_col.append(0.01)

    if fixed_alphas:
        return models, None
    else:
        return models, alpha_per_col




def get_posterior_logprobs(col, bin_mtx, model, pc=np.log(1 / 210)):
    """
    Given a model, calculate the log probability of the observations.
    """
    log_probs = model.predict_log_proba(bin_mtx)
    posterior_logprobs = []
    for idx, res in enumerate(col):
        # Check if the residue was in the training data
        if res in model.classes_:
            # Get index of residue in self.classes_
            res_idx = np.asscalar(np.where(model.classes_ == res)[0])
            # Retrieve the log-probability
            log_prob = log_probs[:, res_idx][idx]
            # Check whether the probability equals 0
            if log_prob != float('-inf'):
                posterior_logprobs.append(log_prob)
            else:  # If it, use a pseudocount
                posterior_logprobs.append(pc)
        else:  # If the residue was not in the training data, use a pseudocount
            posterior_logprobs.append(pc)

    return posterior_logprobs


def calc_degrees_freedom(model):
    """
    Compute the number of non-zero coefficients.
    """
    all_nonzero = []
    n_classes = model.coef_.shape[0]
    # Iterate over the coefficient 3D matrix, going through the coefficients
    # for predicting each class
    for i in range(n_classes):
        # Get the positions of the non-zero coefficients and add them to a list
        nonzero = list(np.where(model.coef_[i, :] != 0)[0])
        all_nonzero.extend(nonzero)
    # The number of degrees of freedom is the number of non-zero coefficients
    # accross all of the 3D matrix, not counting those that are repeated
    dfs = len(np.unique(all_nonzero))

    return dfs


def calc_bic(posterior_logprobs, dfs, n_obs):
    """
    Calculate the Bayesian Information Criterion for a given model.

    Arguments
    ---------
    lhs:    array-like, log-probability of the observations according to the
            models
    dfs:    int, degrees of freedom (non-zero coefficients) of the model
    n_obs:  int, number of observations

    Returns
    ---------
    bic: float, value of the Bayesian Information Criterion
    """
    lhs = np.sum(posterior_logprobs)
    bic = dfs * np.log(n_obs) - 2 * lhs
    return bic

##############
# Null model #
##############


def calc_null_llhs(a1, a2, mode, weights, out_path, iters, pc_null=1 / 2100):
    """
    TODO: unit test
    Calculate the log-likelihood of each sequence pair in the MSAs
    under the assumption of independent evolution: this is, the probability
    of each amino acid depends only on its frequency in a particular column.

    Writes a matrix of probabilities according to the null model to the disk.

    Arguments
    ---------
    a1:        array-like, multiple sequence alignment as numeric matrix
    a2:        array-like, multiple sequence alignment as numeric matrix
    mode:      string. whether the algorithm is using 'hard' or 'soft' EM
    weights:   array-like, values of the hidden variable; used to multiply
               residue probability by the weight of the sequence
    out_path:  string, path to output directory
    iters:     int, number of iterations
    pc_null:   float, pseudocount for residues that did not appear in examples
               used to build the null model

    Returns
    ---------
    null_llhs: array-like, (n_samples)
               Each element is the likelihood of a sequence pair under the null
               model
    """

    # In the case of hard EM, select non-interacting cases
    if mode == 'hard':
        na1, na2, nw = select_noninteracting(a1, a2, weights)
        inv_weights = 1 - np.asarray(nw)
    # In soft EM, just make new references to arrays to simplify things
    elif mode == 'soft':
        na1, na2 = a1, a2
        inv_weights = 1 - np.asarray(weights)

    null_1 = get_null_model(na1, inv_weights, pc_null)
    null_2 = get_null_model(na2, inv_weights, pc_null)
    null_mtx_1 = score_null(a1, null_1, pc_null)
    null_mtx_2 = score_null(a2, null_2, pc_null)

    concat_mtx = np.concatenate((null_mtx_1, null_mtx_2), axis=1)
    np.savetxt(os.path.join(out_path,
                            ''.join(['null_llhs_mtx_', str(iters), '.csv'])),
               concat_mtx, delimiter=',')
    null_llhs = np.sum(concat_mtx, axis=1)

    return null_llhs


def get_null_model(array, weights, pc_null):
    """
    Given a multiple sequence alignment in numeric matrix format, calculates
    the log-probability (weighted according to the values of the hidden
    variables) of each residue in each column. It returns this information in a
    position weight matrix-like format.
    Only call from within calc_null_llhs().

    Arguments
    ---------
    array:   array-like, multiple sequence alignment as numeric matrix
    weights: array-like, values of the hidden variable; used to multiply
             residue probability by the weight of the sequence
    pc_null: float, pseudocount for residues that did not appear in examples
             used to build the null model

    Returns
    ---------
    null_model: list of dicts. Each dictionary contains the null model for each
                column. Each key is an amino acid, as it appears in the numeric
                matrix, and the value is it's weighted probability.
                Can be thought of as a position weight matrix.
    """

    null_model = []
    # Number of observations with z=1 when doing hard-EM,
    # just the sum of the weights when doing soft EM!
    total = np.sum(weights, dtype="longdouble")

    for i, val in enumerate(array.T):
        # Obtain unique column elements
        items = np.unique(val)
        col_dict = {}
        for item in items:
            # Get indexes where this element appears
            idxs = np.where(val == item)[0]
            # Calculate (weighted) element probability
            prob = np.sum([weights[idx]
                           for idx in idxs], dtype="longdouble") / total
            # Check if probability is 0; use a pseudocount in that case
            # This can happen in soft EM in the edge case that all cases that
            # contain the residue have an associated hidden variable equal to 1
            if prob == 0:
                col_dict[str(item)] = np.log(pc_null)
            else:
                col_dict[str(item)] = np.log(prob)
        null_model.append(col_dict)

    return null_model


def score_null(array, null_model, pc_null):
    """
    Calculate the probability of each position of each sequence in a multiple
    sequence alignment according to the null model.

    Arguments
    ---------
    array:      array-like, multiple sequence alignment as numeric matrix
    null_model: list of dicts. Each dictionary contains the null model for each
                column. Each key is an amino acid, as it appears in the numeric
                matrix, and the value is it's weighted probability.
                Can be thought of as a position weight matrix.
    pc_null:    float, pseudocount for residues that did not appear in examples
                used to build the null model
    Returns
    --------
    null_mtx:   array-like, contains residue log-probabilities according to the
                null model for each position in each sequence
    """
    null_mtx = np.zeros_like(array, dtype='float64').T

    # Iterate over array columns and fill log-probability matrix
    for idx, val in enumerate(array.T):
        for j, item in enumerate(val):
            try:
                null_mtx[idx][j] = null_model[idx][str(item)]
            except KeyError:
                # Residue not present in examples used to build the null model;
                # use a pseudocount
                null_mtx[idx][j] = pc_null

    return null_mtx.T


def select_noninteracting(num_mtx_a, num_mtx_b, labels):
    """
    Auxiliary function to select pairs of putatively non-interacting proteins
    before updating the null models in hard EM.

    Arguments
    ---------
    num_mtx_a, num_mtx_b: array-like. Contains a multiple sequence alignment as
                          a numeric matrix with all proteins under
                          consideration, interacting and non-interacting.
                          Each matrix is of size
                          (number of sequences x number of MSA columns)

    labels:               array-like, values of the hidden variables

    Returns
    -------
    nonint_num_a, nonint_num_b: array-like. Selected rows of num_mtx_a and
                                num_mtx_b
    nonint_labels:              array-like. Selected items of labels
    """
    if not all(np.unique(labels) == np.unique(np.array([0, 1]))):
        raise Exception("""Hidden variable values other than 0 or 1 in call to
            corrmut.select_noninteracting(). This function is only called in
            hard EM, where values should be rounded to 0 or 1.""")
    idxs = np.where(np.asarray(labels) == 0)
    nonint_num_a = np.squeeze(np.take(num_mtx_a, idxs, axis=0))
    nonint_num_b = np.squeeze(np.take(num_mtx_b, idxs, axis=0))
    nonint_labels = np.squeeze(np.take(labels, idxs, axis=0))
    return nonint_num_a, nonint_num_b, nonint_labels


######################
# Contact prediction #
######################

def contact_prediction(num_mtx_a, bin_mtx_b, num_mtx_b, bin_mtx_a,
                       labels, seqs_weight, mode, n_jobs, dfmax):
    """
    Function for a final round of contact prediction
    """
    # In hard EM, select observations with z=1; in soft EM, select observations
    # with z>0.5

    int_num_a, int_bin_b, int_num_b, int_bin_a, weights = get_interacting(
        num_mtx_a, bin_mtx_b, num_mtx_b, bin_mtx_a, labels, mode)
    # Remove constant columns that might have appeared; keep the indexes of the
    # constant columns

    # Fit new models of the MSAs, allowing to optimize the regularization
    # strengths again
    models_a, alphas_a = fit_msa_models(
        int_num_a, int_bin_b, mode, seqs_weight, sample_weights=weights, n_jobs=n_jobs,
        dfmax=dfmax)
    models_b, alphas_b = fit_msa_models(
        int_num_b, int_bin_a, mode, seqs_weight, sample_weights=weights, n_jobs=n_jobs,
        dfmax=dfmax)

    couplings, contact_mtx = compute_couplings(models_a, models_b)

    return couplings, contact_mtx


##############################
# Initialization and EM loop #
##############################

def init_model(num_mtx_a, bin_mtx_b, num_mtx_b, bin_mtx_a, seqs_weight, mode,
               init, int_frac, out_dir, n_jobs, dfmax):
    """
    Calculate initial values for the hidden variables before starting the
    EM loop, either randomly or by warm initialization.

    Writes files to disk.

    Arguments
    ---------
    num_mtx_a:  array-like, MSA in numeric matrix form
    bin_mtx_b:  array-like, MSA in binary matrix form
    num_mtx_b:  array-like, MSA in numeric matrix form
    bin_mtx_a:  array-like, MSA in binary matrix form
    mode:       str, whether we are performing hard or soft EM
    init:       str, method to initialize the hidden variables
    int_frac:   float, assumed fraction of interacting proteins
    out_dir:    str, output path
    n_jobs:     int, number of CPUs to use to fit the models
    dfmax:      int, maximum number of degrees of freedom allowed

    Returns
    ---------
    init_labels: array-like, initial values for the hidden variables
    """

    if init == 'random':
        init_labels = get_random_labels(num_mtx_a.shape[0], int_frac, mode)
        alt_llhs = None
        null_llhs = None
        norm_contact_mtx = None
        alphas_a = None
        alphas_b = None
    elif init == 'warm':
        # Fit initial models
        # Assume all sequence pairs are interacting: do not exclude anything,
        # do not pass sample weights
        print('Fitting models for MSA A...')
        models_a, alphas_a = fit_msa_models(num_mtx_a, bin_mtx_b, mode, seqs_weight, n_jobs=n_jobs,
                                            dfmax=dfmax)
        print('Fitting models for MSA B...')
        models_b, alphas_b = fit_msa_models(num_mtx_b, bin_mtx_a, mode, seqs_weight, n_jobs=n_jobs,
                                            dfmax=dfmax)

        couplings, contact_mtx = compute_couplings(models_a, models_b)
        np.savetxt(os.path.join(out_dir, ''.join(
            ['contact_mtx_', 'init', '.csv'])), contact_mtx, delimiter=',')
        norm_contact_mtx = normalize_contact_mtx(contact_mtx)
        np.savetxt(os.path.join(out_dir, ''.join(
            ['norm_contact_mtx_', 'init', '.csv'])), contact_mtx, delimiter=',')

        np.savetxt(os.path.join(out_dir, ''.join(
            ['fixed_alphas_a_iter_', str('init'), '.csv'])), alphas_a)
        np.savetxt(os.path.join(out_dir, ''.join(
            ['fixed_alphas_b_iter_', str('init'), '.csv'])), alphas_b)

        alt_llhs = calc_alt_llhs(num_mtx_a, bin_mtx_b, models_a, num_mtx_b,
                                 bin_mtx_a, models_b, out_dir, iters='init')

        # Observation weights passed to this call is a list of 0s;
        # every sequence pair gets the maximum weight
        null_llhs = calc_null_llhs(num_mtx_a, num_mtx_b, mode,
                                   num_mtx_a.shape[0] * [0],
                                   out_dir, 'init')

        init_labels = update_labels(
            alt_llhs, null_llhs, int_frac, mode=mode)

    return init_labels, alt_llhs, null_llhs, norm_contact_mtx, alphas_a, alphas_b


def em_loop(num_mtx_a, num_mtx_b, bin_mtx_a, bin_mtx_b, labels, seqs_weight,
            int_frac, mode, out_dir, n_jobs,
            max_iters=20, tol=0.005,
            true_labels=None, dfmax=100, fixed_alphas_a=None, fixed_alphas_b=None):
    """
    Main function for carrying out expectation-maximization.

    Note that aside from the returned values, it also writes several files to
    disk when calling certain functions (calc_alt_llhs(),
    calc_null_llhs())

    Arguments
    ---------
    num_mtx_a, num_mtx_b: array-like. Contains a multiple sequence alignment as
                          a numeric matrix with all proteins under
                          consideration, interacting and non-interacting.
    bin_mtx_a, bin_mtx_b: array-like. Contains a multiple sequence alignment as
                          a binary matrix with all proteins under
                          consideration, interacting and non-interacting.
    labels:               array-like, initial values of the hidden variables
    int_frac:             float, prior fraction of interacting proteins
    mode:                 str, whether to perform soft or hard
                          expectation-maximization
    out_dir:              str, output path
    n_jobs:               int, number of CPUs to use to fit the models
    max_iters:            int, maximum number of EM iterations
    tol:                  float, difference threshold for convergence check
    true_labels:          list, contains ground truth labels

    Returns
    ---------
    labels_per_iter:        list, contains the values of the hidden variables
                            for each iteration
    alt_llhs_per_iter:      list, contains alternative model log-likelihoods
                            of all sequence pairs for each iteration
    null_llhs_per_iter:     list, contains null model log-likelihoods of all
                            sequence pairs for each iteration

    """

    labels_per_iter = []
    alt_llhs_per_iter = []
    null_llhs_per_iter = []
    contacts_per_iter = []

    iters = 0
    converged = False

    while (iters < max_iters) and (converged is not True):

        print(f'Starting EM iteration number {iters+1}')

        # =====================================================================
        # Maximization step: update co-evolutionary and null models
        # =====================================================================
        # First, select putatively interacting proteins using the labels
        # derived from the previous iteration. In hard EM, these are sequence
        # pairs with a hidden variable of 1. In soft EM, they are all the
        # sequences passed, no matter the hidden variable value.

        # Fit logistic models for each column
        # Optimize values of alpha only for random initialization
        if iters == 0 and (not fixed_alphas_a and not fixed_alphas_b):
            print('Maximization step: fitting models for MSA A...')

            models_a, fixed_alphas_a = fit_msa_models(num_mtx_a, bin_mtx_b,
                                                      mode,
                                                      seqs_weight,
                                                      sample_weights=labels,
                                                      n_jobs=n_jobs,
                                                      dfmax=dfmax)
            print('Maximization step: fitting models for MSA B...')
            models_b, fixed_alphas_b = fit_msa_models(num_mtx_b, bin_mtx_a,
                                                      mode,
                                                      seqs_weight,
                                                      sample_weights=labels,
                                                      n_jobs=n_jobs,
                                                      dfmax=dfmax)

            # Dump values of alpha
            np.savetxt(os.path.join(out_dir, ''.join(
                ['fixed_alphas_a_iter_', str(iters), '.csv'])), fixed_alphas_a)
            np.savetxt(os.path.join(out_dir, ''.join(
                ['fixed_alphas_b_iter_', str(iters), '.csv'])), fixed_alphas_b)

        else:
            print('Maximization step: fitting models for MSA A...')
            models_a, _ = fit_msa_models(num_mtx_a, bin_mtx_b, mode,
                                        seqs_weight,
                                         fixed_alphas=fixed_alphas_a,
                                         sample_weights=labels,
                                         n_jobs=n_jobs, dfmax=dfmax)
            print('Maximization step: fitting models for MSA B...')
            models_b, _ = fit_msa_models(num_mtx_b, bin_mtx_a, mode,
                                        seqs_weight,
                                         fixed_alphas=fixed_alphas_b,
                                         sample_weights=labels,
                                         n_jobs=n_jobs, dfmax=dfmax)

        # =====================================================================
        # Expectation step: update labels based on the new co-evolutionary and
        # null models
        # =====================================================================
        print('Expectation step: updating hidden variables...')
        # Calculate likelihoods for the alternative model using all available
        # pairs of proteins

        # Use these to update the alternative and null model
        alt_llhs = calc_alt_llhs(num_mtx_a, bin_mtx_b, models_a, num_mtx_b,
                                 bin_mtx_a, models_b, out_dir, iters)
        null_llhs = calc_null_llhs(num_mtx_a, num_mtx_b, mode, labels, out_dir,
                                   iters)

        # Save previous labels for convergence calculations; update labels
        pre_labels = labels
        labels = update_labels(alt_llhs, null_llhs,
                               int_frac, mode=mode)

        # Predict contacts and dump contact matrix
        couplings, contact_mtx = compute_couplings(models_a, models_b)
        np.savetxt(os.path.join(out_dir, ''.join(
            ['contact_mtx_', str(iters), '.csv'])), contact_mtx, delimiter=',')
        norm_contact_mtx = normalize_contact_mtx(contact_mtx)
        np.savetxt(os.path.join(out_dir, ''.join(
            ['norm_contact_mtx_', str(iters), '.csv'])),
            norm_contact_mtx, delimiter=',')

        # Add new information to function output
        labels_per_iter.append(labels)
        alt_llhs_per_iter.append(alt_llhs)
        null_llhs_per_iter.append(null_llhs)
        contacts_per_iter.append(norm_contact_mtx)

        # Check whether the EM has converged
        converged = has_converged(labels, pre_labels, mode, tol)
        iters += 1

    return labels_per_iter, alt_llhs_per_iter, null_llhs_per_iter, contacts_per_iter


def em_wrapper(num_mtx_a, num_mtx_b, bin_mtx_a, bin_mtx_b, n_starts,
               int_frac, mode, seqs_weight, results_dir, n_jobs, dfmax, test,
               em_args, true_labels=None):
    """
    Function for repeated calling of the expectation-maximization loop.
    This is used to carry out multiple random starts.

    Writes files to disk.

    Arguments
    ---------
    num_mtx_a, num_mtx_b: array-like. Contains a multiple sequence alignment as
                          a numeric matrix with all proteins under
                          consideration, interacting and non-interacting.
    bin_mtx_a, bin_mtx_b: array-like. Contains a multiple sequence alignment as
                          a binary matrix with all proteins under
                          consideration, interacting and non-interacting.
    n_starts:             int, number of random starts
    int_frac:             float, prior fraction of interacting proteins
    mode:                 str, whether to perform soft or hard
                          expectation-maximization
    results_dir:          str, path for input/output
    test:                 str, whether to activate test options or not
    n_jobs:               int, number of CPUs to use to fit the models
    em_args:              dict, contains keyword arguments for em_loop()
    true_labels:          list, contains ground truth labels


    Returns
    -------
    None
    """
    for i in range(n_starts):
        print('Now in random start: {}'.format(str(i)))
        start_path = os.path.join(results_dir, ''.join(['n_start', str(i)]))
        checks_path = os.path.join(start_path, 'checks')
        os.mkdir(start_path)
        os.mkdir(checks_path)

        print('Initializing model...')
        init_labels, *_ = init_model(num_mtx_a, bin_mtx_b, num_mtx_b,
                                     bin_mtx_a, seqs_weight, mode, 'random',
                                     int_frac, checks_path, n_jobs, dfmax)

        print('Starting EM loop...')
        labels_per_iter, alt_llhs_per_iter, \
            null_llhs_per_iter, contacts_per_iter = em_loop(num_mtx_a, num_mtx_b,
                                                            bin_mtx_a, bin_mtx_b,
                                                            init_labels,
                                                            seqs_weight,
                                                            int_frac, mode,
                                                            checks_path, n_jobs,
                                                            **em_args)
        alt_int_per_iter, \
            null_nonint_per_iter = compute_llhs(labels_per_iter,
                                                alt_llhs_per_iter,
                                                null_llhs_per_iter)
        # Create output for the current iteration
        if test:
            alt_true_per_iter, \
                null_true_per_iter = compute_llhs(
                    [em_args['true_labels']] * len(labels_per_iter),
                    alt_llhs_per_iter, null_llhs_per_iter)
            labels_per_iter.insert(0, init_labels)
            output.create_output(labels_per_iter, alt_llhs_per_iter,
                                 null_llhs_per_iter, alt_int_per_iter,
                                 null_nonint_per_iter, mode, start_path, test,
                                 true_labels, alt_true_per_iter,
                                 null_true_per_iter)
        else:
            labels_per_iter.insert(0, init_labels)
            output.create_output(labels_per_iter, alt_llhs_per_iter,
                                 null_llhs_per_iter, alt_int_per_iter,
                                 null_nonint_per_iter, mode, start_path, test)


def compute_llhs(labels_per_iter, alt_llhs_per_iter, null_llhs_per_iter):
    """
    Function for calculating log-likelihoods of the data according to the models.

    In the case of hard EM, these are the alternative model log-likelihoods of
    putatively interacting sequence pairs and the null model log-likelihood
    of putatively non-interacting sequence pairs.

    In the case of soft EM, these are the alternative model log-likelihoods of
    all sequence pairs, weighted according to their associated hidden variable,
    and the null model model log-likelihoods of all sequence pairs, weighted
    according to the their associated hidden variable.

    Both of these amount to a matrix multiplication of the hidden variable
    values and the log-likelihoods.

    Arguments
    ---------
    labels_per_iter:        list, contains the values of the hidden variables
                            for each iteration
    alt_llhs_per_iter:      list, contains alternative model log-likelihoods
                            of all sequence pairs for each iteration
    null_llhs_per_iter:     list, contains null model log-likelihoods of all
                            sequence pairs for each iteration

    Returns
    -------
    alt_int_per_iter:       list, contains alternative model log-likelihoods of
                            putatively interacting sequence pair for each
                            iteration
    null_nonint_per_iter:   list, contains null model log-likelihoods of
                            putatively non-interacting sequence pairs for each
                            iterations
    """

    zs = np.array(labels_per_iter)
    alts = np.array(alt_llhs_per_iter)
    nulls = np.array(null_llhs_per_iter)

    alt_int_llh = np.multiply(alts, zs)
    null_nonint_llh = np.multiply(nulls, (1-zs))

    return alt_int_llh, null_nonint_llh
