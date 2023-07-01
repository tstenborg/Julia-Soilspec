# Title:       cssfunction
#
# Description: Assesses the adequate calibration set size for infrared
#              spectroscopy.
#
# Parameters:
#   S            A matrix of the scores of the principal components.
#   k            A vector containing the sample set sizes to be evaluated.
#   method       The sampling algorithm. Options are:
#                - "kss" (Kennard-Stone Sampling),
#                - "kms" (K-Means Sampling), the default,
#                - "clhs" (Conditioned Latin Hypercube Sampling).
#   repetitions  The number of times that the sampling must be carried out for
#                for each sample size to be evaluated. The results of the
#                final msd is the average of the ones obtained at each
#                iteration. Note that since the "kss" method is deterministic
#                and always returns the same results, there is no need for
#                repetitions.
#   n            The number of equally spaced points at which the probability
#                densities are to be estimated (see the density function of the
#                R package, stats).
#   from, to     A vector of the left and right-most points of the grid at which
#                the densities are to be estimated. Default is the minimums and
#                maximums of the variables in S.
#   bw           A vector containing the the smoothing bandwidth to be use for
#                the probability densities (see density function of the package
#                stats).
#   ...          Arguments to be passed to the calibration sampling algorithms,
#                i.e. additional aruments to be used for the clhs, kenStone or
#                naes functions which run inside this function.
#
# Info:        A function for assessing the adequate calibration set size for
#              - Kennard-Stone Sampling,
#              - K-Means Sampling,
#              - Conditioned Latin Hypercube Sampling.
#
#              This function works by comparing the probability density function
#              (pdf) of the population and the pdf of the sample set in order to
#              assess the representativeness of the sample set. See
#              Ramirez-Lopez et al. (2014) for more details.
#
# Returns:     A table with the following columns
#              - css: the sample set size (k),
#              - msd,
#              - msd_sd: the standard deviation of the msd for all the
#                        repetitions (doesn't apply to "kss" since it always
#                        returns the same results).
#
# Reference:
#  Ramirez-Lopez, L., Schmidt, K., Behrens, T., van Wesemael, B.,
#    Demattê, J. A., Scholten, T. (2014), Sampling optimal calibration sets in
#    soil infrared spectroscopy. Geoderma, 226, 140-150.

function css(S, k, method = "kms", repetitions = 10, n = 512, from, to, bws, ...)

  requireNamespace("clhs")
  requireNamespace("matrixStats")

  if missing(from)
    min_sc = matrixStats::colMins(S)
  else
    min_sc = from
  end
  if missing(from)
    max_sc = matrixStats::colMaxs(S)
  else
    max_sc = to
  end

  if length(min_sc) != ncol(S)
    stop("Argument 'from' needs to be a vector with the same number of variables (columns) as in S ")
  end

  if length(max_sc) != ncol(S)
    stop("Argument 'to' needs to be a vector with the same number of variables (columns) as in S ")
  end

  if missing(bw)
    d_bandwidths = apply(S, 2, bw.nrd0)
  else
    d_bandwidths = bw
  end

  if length(d_bandwidths) != ncol(S)
    stop("Argument 'bw' needs to be a vector with the same number of variables (columns) as in S ")
  end


  # Matrix where the density values will be stored.
  sc_dens = NULL
  for i in 1:length(min_sc)
    i_sc_dens = data.frame(x = seq(min_sc[i], max_sc[i], length = n),
                           densc = rep(NA, n), pc = paste("PC-", i, sep = ""))
    sc_dens = rbind(sc_dens, i_sc_dens)
  end

  # Estimate the density distribution of each variable.
  names(d_bandwidths) = colnames(S)
  for i in 1:length(min_sc)
    idsty = density(S[:,i],
                    bw = "nrd0",
                    n = n, from = min_sc[i], to = max_sc[i],
                    kernel = "gaussian")
    sc_dens[sc_dens$pc == paste("PC-", i, sep = ""), "densc"] = idsty$y
    #d_bandwidths[i] <- idsty$bw
  end

  # 3. Define the different sample set sizes.

  # 4. Sample with the specified algorithm.
  # For each sample size the sampling is repeated 10 times and the differences
  # between the density distribution of the whole set is compared against the
  # density distribution of the sample set.
  #
  results_ss = data.frame(css = k,
                          msd = rep(NA, length(k))
  )
  if method == "kss" & repetitions > 1
    warning("For Kennard-Stone Sampling repetitions are not necessary. Only one repetition will be executed.")
    repetitions = 1
  end

  for i in 1:repetitions
    results_ss[:, -1] = NA
    fn = paste(method, "_temp_results_rep", i,".txt", sep = "")

    for j in 1:length(k)
      set.seed(j)

      if method == "kms"
        i_calidx = prospectr::naes(X = S,
                         k = k[j],
                         method = 0,
                         .center = FALSE,
                         .scale = FALSE, ...)$model
      end

      if method == "kss"
        i_calidx = prospectr::kenStone(X = S,
                             k = k[j],
                             metric = "mahal",
                             .center = FALSE,
                             .scale = FALSE, ...)$model
      end

      if method == "clhs"
        i_calidx = clhs::clhs(as.data.frame(S),
                         size = k[j],
                         simple = TRUE,
                         progress = FALSE, ...)
      end


      m_sc_dens = sc_dens
      for m in 1:length(min_sc)

        # Use the same bandwidth (bw) as in the whole set of candidates.
        slc = sc_dens$pc == paste("PC-", m, sep = "")
        m_sc_dens[slc, "densc"] = density(S[i_calidx, m],
                                          bw = d_bandwidths[m],
                                          n = n, from = min_sc[m],
                                          to = max_sc[m],
                                          kernel = "gaussian")$y
      end
      results_ss$msd[j] = mean((m_sc_dens$densc - sc_dens$densc)^2, na.rm = T)

      # Write the results to a table.
      if method == "kss"
        write.table(results_ss,
                    file = paste(method, "_final_results_.txt", sep = ""),
                    row.names = FALSE, sep = "\t")
        final_ss = results_ss

      else
        write.table(results_ss,
                    file = fn,
                    row.names = FALSE, sep = "\t")
      end
    end
  end

  if method != "kss"
    # 5. Read the iteration results from the generated files and compute the
    #    mean of the iterations.
    nmsreps = paste(method, "_temp_results_rep", 1:repetitions, ".txt",
                    sep = "")
    final_ss = 0
    for i in nmsreps
      iter = which(i == nmsreps)
      results_ss = read.table(i, header = T, sep = "\t")
      final_ss = final_ss + results_ss
      if (iter == length(nmsreps))
        final_ss = final_ss/iter
      end
    end

    # 6. Read the iteration results from the generated files and compute the
    #    standard deviation of the iterations.
    final_ss_sd = 0
    for i in nmsreps
      iter = which(i == nmsreps)
      results_ss = read.table(i, header = T, sep = "\t")
      final_ss_sd = (results_ss - final_ss_sd)^2
      if (iter == length(nmsreps))
        final_ss_sd = (final_ss_sd / iter)^0.5
      end
    end
    if repetitions > 2
      final_ss = data.frame(final_ss, msd_sd = final_ss_sd[:, 2])
    end
    write.table(final_ss, file = paste(method, "_final_results_.txt", sep = ""),
                sep = "\t", row.names = FALSE)
  end
  return(final_ss)
end
