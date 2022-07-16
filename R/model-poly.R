#' Model list of sequences
#'
#' Analyses a concurrent sequences using a PPM model.
#'
#' @param model
#' A PPM model object as produced by (for example)
#' \code{\link{new_ppm_simple}} or \code{\link{new_ppm_decay}}.
#'
#' @param seqs
#' A list of integer vectors defining the input
#' sequences (equivalently a numeric vector containing solely integers,
#' or a factor vector, both of which which will be coerced to integer vectors).
#'
#' @param times
#' (NULL or list of numeric vectors for each sequence)
#' Time points corresponding to each element of the sequence.
#' Used to align events in consecutive sequences and for decay-based models.
#'
#' @param zero_indexed
#' (Logical scalar)
#' Whether or not sequences in \code{seqs} argument are 0-indexed
#' (i.e. drawn from an alphabet with a minimum value of 0).
#' If FALSE, it is assumed that the sequence is 1-indexed
#' (i.e. drawn from an alphabet with a minimum value of 1).
#'
#' @param train
#' (Logical scalar or vector)
#' Whether or not the model should learn from the incoming sequences,
#' or individual sequences.
#'
#' @param predict
#' (Logical scalar or vector)
#' Whether or not to generate predictions for each element of
#' the incoming sequences, or individual sequences.
#'
#' @param return_distribution
#' (Logical scalar)
#' Whether or not to return the conditional distribution over each
#' potential continuation as part of the model output
#' (ignored if \code{predict = FALSE}).
#'
#' @param return_entropy
#' (Logical scalar)
#' Whether or not to return the entropy of each event prediction
#' (ignored if \code{predict = FALSE}).
#'
#' @md
#' @rdname model_poly
#' @export
model_poly <- function(model,
                       seqs,
                       times = NULL,
                       zero_indexed = FALSE,
                       train = TRUE,
                       predict = TRUE,
                       return_distribution = TRUE,
                       return_entropy = TRUE) {
  if (!is.null(times)
      && !identical(lapply(seqs, length), lapply(times, length))) {
    stop("if times are not NULL, equal numbers of sequences and times ",
         "must be given, with an equal number of items in respective ",
         "sequences/times.")
  }

  if (any(lapply(seqs, is.character))) {
    stop("sequences in 'seqs' cannot be a character vector; ",
         "please provide a factor representation instead.")
  }

  seqs <- lapply(seqs, as.integer)

  stopifnot(is_ppm(model))
  lapply(seqs, (function(x) checkmate::qassert(x, "X")))
  checkmate::qassert(zero_indexed, "B1")
  checkmate::qassert(train, "B")
  checkmate::qassert(predict, "B")
  checkmate::qassert(return_distribution, "B1")
  checkmate::qassert(return_entropy, "B1")
  stopifnot(is.null(times) || all(lapply(times, is.numeric)))

  for (seq in seqs) {
    if (is.factor(seq)
        && length(model$alphabet_levels > 0)
        && !identical(levels(seq), model$alphabet_levels)) {
      warning("sequence's factor levels seemed inconsistent with ",
              "model$alphabet_levels")
    }
  }

  if (any(unlist(lapply(times, (function(x) diff(x))) < 0))) {
    stop("decreasing values of time are not permitted")
  }

  if (is.null(times)) {
    times <- lapply(seqs, (function(x) seq_len(length(x))))
  }

  if (length(train) == 1) train <- rep(train, length(seqs))

  if (length(predict) == 1) predict <- rep(predict, length(seqs))

  if (!zero_indexed) seqs <- lapply(seqs, function(x) x - 1L)

  res <- model$model_poly(
    x = seqs,
    times = times,
    train = train,
    predict = predict,
    return_distribution = return_distribution,
    return_entropy = return_entropy
  )

  res
}
