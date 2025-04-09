#' PCM Friction Model
#'
#' Implements the PCM friction model for runout simulation.
#'
#' @param mu Numeric. Sliding friction coefficient (default: 0.1).
#' @param md Numeric. Mass-to-drag ratio (default: 40).
#' @param v_p Numeric. Initial velocity (m/s) (default: 1).
#' @param theta_p Numeric. Slope angle of previous grid cell (degrees) (default: 30).
#' @param theta_i Numeric. Slope angle of current grid cell (degrees) (default: 20).
#' @param l Numeric. Distance between grid cells (m) (default: 12.5).
#'
#' @return Numeric. Computed velocity in the local grid cell. Returns `NaN` if the velocity is not physically possible (e.g., stopping condition).
#'
#' @details
#' The PCM model calculates velocity propagation across a terrain grid based on friction, slope, and mass-to-drag ratio. 
#' It includes velocity correction for concave slope transitions as per Wichmann (2017).
#'
#' @references
#' Wichmann, V.: The Gravitational Process Path (GPP) model (v1.0) – a GIS-based simulation framework for gravitational processes, Geosci. Model Dev., 10, 3309–3327, https://doi.org/10.5194/gmd-10-3309-2017, 2017.
#' Perla, R., Cheng, T. T., and McClung, D. M.: A Two–Parameter Model of Snow–Avalanche Motion, J. Glaciol., 26, 197–207, https://doi.org/10.3189/S002214300001073X, 1980.
#'
#' @examples
#' pcm(mu = 0.1, md = 40, v_p = 2, theta_p = 35, theta_i = 25, l = 12.5)
#'
#' @export
pcm <- function(mu = 0.1, md = 40, v_p = 1, theta_p = 30, theta_i = 20, l = 12.5) {
  g <- 9.80665  # Gravitational acceleration (m/s^2)
  
  alpha <- g * (sin(theta_i * pi / 180) - mu * cos(theta_i * pi / 180))  # Acceleration
  beta <- -2 * l / md  # Adjustment factor
  
  # Velocity correction for concave transitions
  delta_theta <- if (theta_p > theta_i) theta_p - theta_i else 0
  
  # Compute velocity
  velocity_sq <- alpha * md * (1 - exp(beta)) + (v_p^2 * exp(beta) * cos(delta_theta * pi / 180))
  
  v_i <- if (velocity_sq < 0) NaN else sqrt(velocity_sq)
  
  return(v_i)
}
