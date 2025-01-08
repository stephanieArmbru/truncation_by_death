

## AR(1) Stationarity Investigation

## NON-STATIONARY
# Set parameters
phi <- 2              # AR(1) coefficient
sigma_epsilon <- 1       # Standard deviation of white noise
n <- 100                 # Number of time steps

# Simulate AR(1) using arima.sim
y <- arima.sim(model = list(ar = phi), n = n, sd = sigma_epsilon)

# Plot the result
plot(y, type = "l",
     main = "Simulated AR(1) Process",
     xlab = "Time", ylab = "Y_t")

## STATIONARITY
# Set parameters
phi <- 0.9             # AR(1) coefficient
sigma_epsilon <- 1       # Standard deviation of white noise
n <- 100                 # Number of time steps

# Simulate AR(1) using arima.sim
y <- arima.sim(model = list(ar = phi), n = n, sd = sigma_epsilon)

# Plot the result
plot(y, type = "l",
     main = "Simulated AR(1) Process",
     xlab = "Time", ylab = "Y_t")
