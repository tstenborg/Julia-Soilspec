# Title:       myImagePlot
#
# Description: Plot a projection matrix for external parameter othogonalisation
#              (EPO).
#
# Parameters:
#   x          The projection matrix.
#   zlim       A two-element array holding the minimum (element one) and maximum
#                (element two) values in x that should be plotted.
#   title      Title for the plot. Optional.
#   xLabels    Labels for the plot's x-axis. Optional.
#   yLabels    Labels for the plot's y-axis. Optional.
#
# Returns:     Image of the projection matrix.

# A function for plotting a matrix.

function myImagePlot(x, zlim = nothing, title = nothing, xLabels = nothing,
                     yLabels = nothing)

  # Bound values to be plotted.
  if isnothing(zlim)
    min = min(x)
    max = max(x)
  else
    min = zlim[1]
    max = zlim[2]
  end

  # Allow for plot labels not passed in.
  if isnothing(xLabels)
    xLabels = [1:size(x, 2)]   # Number of columns.
  end
  if isnothing(yLabels)
    yLabels = [1:size(x, 1)]   # Number of rows.
  end

  layout(matrix(data = [1, 2], nrow = 1, ncol = 2),
         widths = [4, 1],
         heights = [1, 1])

  # Red and green range from 0 to 1 while blue ranges from 1 to 0.
  ColorRamp = rgb(seq(0, 1, length = 256),  # Red.
                  seq(0, 1, length = 256),  # Green.
                  seq(1, 0, length = 256))  # Blue.
  ColorLevels = seq(min, max, length = length(ColorRamp))

  # Reverse Y axis.
  reverse = nrow(x) : 1
  yLabels = yLabels[reverse]
  x = x[reverse,]

  # Data Map.
  par(mar = [3, 5, 2.5, 2])
  image(1:length(xLabels), 1:length(yLabels), t(x), col = ColorRamp, xlab = "",
    ylab = "", axes = FALSE, zlim = [min, max])
  if !is.null(title)
    title(main = title)
  end
  axis(BELOW = 1, at = 1:length(xLabels), labels = xLabels, cex.axis = 0.7)
  axis(LEFT = 2, at = 1:length(yLabels), labels = yLabels,
    las = HORIZONTAL <- 1, cex.axis = 0.7)

  # Color Scale.
  par(mar = c(3, 2.5, 2.5, 2))
  image(1, ColorLevels,
    matrix(data = ColorLevels, ncol = length(ColorLevels), nrow = 1),
    col = ColorRamp, xlab = "", ylab = "", xaxt = "n")

  layout(1)
end
