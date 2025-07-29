
rockfall <- function(mu = 0.1, v_p = 1, h_i = 300, h_p = 400, 
                     K=0.25, method = 'rolling', l = 12.5){
  
  # K: reduction percentage 0.25 = 75%
  # mu: either 'rolling' or 'sliding' friction coefficient depending on method
  # method: 'sliding' or 'rolling'
  # Scheidegger 1975
  # see Wichmann 2017 in GMD for more details
  
  g <- 9.80665  # Gravitational acceleration (m/s^2)
  h_f <- h_p - h_i # height difference
  
  if(is.null(v_p)){
    # Initial velocity from start cell
    
    v_i <- sqrt(2*g*h_f*K)
  
  } else {

    if(method == "rolling"){
      
      vel <- v_p^2 + 2 * g * (h_f - mu * l)
    
      v_i <- if (vel < 0 | is.nan(vel)) NaN else sqrt(vel)
      
    } else if(method == "sliding"){
      
      vel <- v_p^2 + 10/7 * g * (h_f - mu * l)
      
      v_i <- if (vel < 0 | is.nan(vel)) NaN else sqrt(vel)
      
    }
  }
  
  return(v_i)
  
}

library(ggplot2)

# Distance along the profile (0 to 100 m)
x <- 0:100

# Parameters
cliff_height <- 30
steep_part <- -cliff_height * (1 - exp(-0.15 * x))  # Exponential drop

# Flattening towards the bottom
flattened <- -cliff_height / (1 + exp(-0.3 * (x - 30)))  # Logistic curve

# Combine and blend (weighted sum)
z <- 0.7 * steep_part + 0.3 * flattened

# Add talus slope (gentle), starting around x = 40
z[x > 40] <- z[x > 40] + (x[x > 40] - 40) * (-0.15)

# Make slope flat after x = 60
z[x > 60] <- z[60]

z = z + abs(min(z))

# Plot
profile_df <- data.frame(Distance = x, Elevation = z)

ggplot(profile_df, aes(x = Distance, y = Elevation)) +
  geom_line(color = "black", size = 0.6) +
  labs(title = "Synthetic Hillslope Profile (Cliff to Talus Slope)",
       x = "Distance (m)", y = "Elevation (m)") +
  theme_minimal() +
  coord_fixed(ratio = 1)

Ks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
k_plots <- list()

for (k in 1:length(Ks)){
  
  print(k)
  vel_store <- rep(NA, times = 100)
  
  v_p = NULL
  h_p = max(z) + 5
  h_i = z[1]
  
  for(i in 1:length(x)){
    h_i = z[i]
    vel <- rockfall(mu = 0.6, v_p = v_p, h_i = h_i, h_p = h_p, 
                    K = Ks[k], method = 'sliding', l = 1)
    
    v_p = vel
    h_p = h_i
    
    if(is.nan(v_p)){
      vel = 0
    }
     vel_store[[i]] <- vel
    
  }
  
  
  # Find where runout stops
  stop_index <- which(vel_store == 0)[1]
  stop_point <- data.frame(Distance = x[stop_index], Elevation = z[stop_index])
  
  k_plots[[k]] <- ggplot(profile_df, aes(x = Distance, y = Elevation)) +
    geom_line(color = "black", size = 0.6) +
    labs(title = paste("K = ", Ks[k]),
         x = "Distance (m)", y = "Elevation (m)") +
    geom_point(data = stop_point) +
    theme_minimal() +
    coord_fixed(ratio = 1)
  
}

k_plots
