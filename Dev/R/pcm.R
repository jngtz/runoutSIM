
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
