# Title:       chBLCext
#
# Description: Fits a convex hull for regions of interest in a spectrum.
#
# Parameters:
#   spectra    [numeric] An array of the spectrum's intensities.
#   type       [character] Either "R" or "A" to indicate the type of the input
#                pectrum. Use "R" for reflectance, "A" for absorption.
#   wav        [numeric] An array of the spectrum's wavelengths. The number of
#                values should match those in the "spectra" parameter.
#
# Returns:     A fitted convex hull.
#
# Info:        R version inspired by the "continuumRemoval" function of the
#              prospectr package.
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function chBLCext(spectra, type = c("R", "A"), wav, ...)

  # If absorbance, take the inverse to find the convex hull points.
  if type == "A"
    spectra = 1 / spectra
  end

  function cHullFun(x, wav)
    cHull = sort(chull(c(wav[1] - 1, wav, wav[length(wav)] + 1), c(0, x, 0)))
    cHull = cHull[-c(1, length(cHull))] - 1
    return(approx(x = wav[cHull], y = x[cHull], xout = wav,
      method = "linear")$y)
  end
  cont = cHullFun(spectra, wav)

  if type == "A"
    # Subtraction: absorbance (Fig. 5 Clark & Roush (1984)).
    hullSpectra = 1 + spectra - cont
  else
    # Division: reflectance (Fig. 5 Clark & Roush (1984)).
    hullSpectra = spectra / cont
  end

  if type == "A"
    # Back transform the CR spectra.
    hullSpectra = 1 / hullSpectra
    # Back transform to show the CH line.
    cont = 1 / cont
    # Back transform the spectra.
    spectra = 1 / spectra
  end

  # Prepare xy for the polygon.
  pol = cbind(wav, as.numeric(hullSpectra[1,]))

  retval = list(wave = wav,
                cHull = hullSpectra,
                rawSpec = spectra,
                continuum = cont,
                polygon = pol)
  return(retval)
end
