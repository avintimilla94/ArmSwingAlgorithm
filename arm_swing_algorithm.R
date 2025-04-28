# arm_swing_algorithm.R

#The following script calculates several arm swing parameters from IMU(s) on the wrist during walking. Original code was written by Elke Warmerdam, Kiel University on MatLab and converted to R code by Antonio Vintimilla, PT, DPT, PhD. 
#Original source code can be found at <https://github.com/EWarmerdam/ArmSwingAlgorithm/blob/master/arm_swing_algorithm.m> 
#The following code and definitions have been adapted directly from the aformentioned source.

# arm_swing_algorithm Calculates several arm swing parameters from IMU(s) on the wrist during walking
#
#   Definitions:
#   swing = either a backward or a forward swing
#   cycle = 2 swings (backward and forward)

#   input:
#          ang_vel     = raw angular velocity (rad/s) from IMU (x forward (direction of thumb), y left and z
#                        vertical in n pose (Nx3 when using one IMU and Nx6 when using two IMUs
#                        (first three channels left, last three right arm))
#          fs          = sample frequency of the IMU
#          TH_min_ampl = Amplitude threshold. Every swing below this
#                        threshold will be discarded
#          ...         = Optional input parameter for the number of swings
#                        that should be removed at the start and end of the walk when
#                        only analysing steady state walking ([start, end])
#   output:
#           amplitude  = amplitude per swing (range of motion) [deg]
#           pk_ang_vel = peak angular velocity per swing [deg/s]
#           start_idx  = sample number at which a swing starts
#           end_idx    = sample number at which a swing ends
#     regularity_angle = regularity of the angular signal (similarity of
#                        neighbouring swings; 1 = similar to neighbouring swings) (0-1)
#   regularity_ang_vel = regularity of the angular velocity signal
#       pk_vel_forward = average peak angular velocity of all the forward swings
#      pk_vel_backward = average peak angular velocity of all the backward swings
#      perc_time_swing = time during walking that there were swings detected [%]
#            frequency = frequency of the arm cycle [Hz]
#      perc_both_swing = percentage time during walking bout that there is
#                        arm swing detected in both arms [%]
#  amplitude_asymmetry = asymmetry of the amplitude between left and right swings (0% means no asymmetry) [%]
# peak_velocity_asymmetry= asymmetry of the peak angular velocity between left and right swings [%]
#     coordination_max = coordination of the timing between left and right swings (1 when the arms
#                        move exactly out of phase with each other) (0-1)
#

# Required libraries
library(signal)
library(pracma)
library(tibble)

# ----- Tukey Window (needed for autocorrelation) -----
tukeywin <- function(n, r = 0.5) {
  if (r <= 0) return(rep(1, n))
  if (r >= 1) return(hanning(n))
  w <- numeric(n)
  per <- floor(r * (n - 1) / 2)
  for (i in 1:n) {
    if (i <= per) {
      w[i] <- 0.5 * (1 + cos(pi * ((2 * i) / (r * (n - 1)) - 1)))
    } else if (i <= (n - per)) {
      w[i] <- 1
    } else {
      w[i] <- 0.5 * (1 + cos(pi * ((2 * i) / (r * (n - 1)) - (2 / r) + 1)))
    }
  }
  return(w)
}

# ----- Autocorrelation-based Regularity -----
auto_cor_wrist <- function(signal, fs) {
  r <- round(20 / 1000 * fs)
  N_win <- round(4500 / 1000 * fs)
  N_acf <- round(2500 / 1000 * fs)
  acf_start_idx <- round(300 / 1000 * fs)
  
  sig_ext <- c(rep(0, round(N_win / 2)), signal, rep(0, round(N_win / 2)))
  acf_vals <- c()
  
  for (i in seq(1, length(sig_ext) - N_win, by = r)) {
    win <- sig_ext[i:(i + N_win - 1)]
    win <- win * tukeywin(N_win, 0.3)
    acf <- acf(win, lag.max = N_acf, plot = FALSE, type = "correlation")$acf
    acf_vals <- c(acf_vals, max(acf[(acf_start_idx + 1):length(acf)], na.rm = TRUE))
  }
  
  mean(acf_vals, na.rm = TRUE)
}

# ----- Core Processing for Each Arm -----
process_arm <- function(data_arm, fs, TH_min_ampl) {
  pca_input <- scale(data_arm[, 1:2], center = TRUE, scale = FALSE)
  pca_result <- prcomp(pca_input)
  ang_vel_pca <- pca_result$x[, 1]
  if (pca_result$rotation[2, 1] < 0) ang_vel_pca <- -ang_vel_pca
  
  angle <- cumtrapz(ang_vel_pca) / fs
  wts <- c(1 / (2 * fs), rep(1 / fs, fs - 1), 1 / (2 * fs))
  mov_avg <- stats::filter(angle, wts, sides = 2)
  idx <- (fs %/% 2 + 1):(length(angle) - fs %/% 2)
  angle_pca <- angle[idx] - mov_avg[idx]
  ang_vel_pca <- ang_vel_pca[idx]
  
  peaks <- findpeaks(angle_pca, minpeakheight = TH_min_ampl)
  if (is.null(peaks) || nrow(peaks) < 2) {
    return(list(
      amplitude = NA, pk_ang_vel = NA, start_idx = NA, end_idx = NA,
      regularity_angle = NA, regularity_ang_vel = NA, pk_vel_forward = NA,
      pk_vel_backward = NA, perc_time_swing = NA, frequency = NA
    ))
  }
  
  pk_idx <- peaks[, 2]
  ampl <- c()
  pk_ang_vel <- c()
  start_idx <- c()
  end_idx <- c()
  
  for (i in 1:(length(pk_idx) - 1)) {
    idx1 <- pk_idx[i]
    idx2 <- pk_idx[i + 1]
    ampl <- c(ampl, abs(angle_pca[idx2] - angle_pca[idx1]))
    pk_ang_vel <- c(pk_ang_vel, max(abs(ang_vel_pca[idx1:idx2])))
    start_idx <- c(start_idx, idx1)
    end_idx <- c(end_idx, idx2)
  }
  
  reg_angle <- auto_cor_wrist(angle_pca, fs)
  reg_vel <- auto_cor_wrist(ang_vel_pca, fs)
  pk_vel_forward <- mean(pk_ang_vel[seq(1, length(pk_ang_vel), 2)])
  pk_vel_backward <- mean(pk_ang_vel[seq(2, length(pk_ang_vel), 2)])
  perc_time_swing <- sum(end_idx - start_idx) / (end_idx[length(end_idx)] - start_idx[1]) * 100
  freq <- length(ampl) / ((end_idx[length(end_idx)] - start_idx[1]) / fs)
  
  list(
    amplitude = ampl,
    pk_ang_vel = pk_ang_vel,
    start_idx = start_idx,
    end_idx = end_idx,
    regularity_angle = reg_angle,
    regularity_ang_vel = reg_vel,
    pk_vel_forward = pk_vel_forward,
    pk_vel_backward = pk_vel_backward,
    perc_time_swing = perc_time_swing,
    frequency = freq
  )
}

# ----- Main Algorithm Wrapper -----
arm_swing_algorithm <- function(ang_vel, fs, TH_min_ampl = 10) {
  if (!is.matrix(ang_vel)) ang_vel <- as.matrix(ang_vel)
  
  butter_lp <- butter(2, 3 / (0.5 * fs), type = "low")
  ang_vel_filt <- apply(ang_vel, 2, function(col) filtfilt(butter_lp, col))
  ang_vel_deg <- ang_vel_filt * 180 / pi
  
  if (ncol(ang_vel_deg) == 3) {
    nr_imu <- 1
  } else if (ncol(ang_vel_deg) == 6) {
    nr_imu <- 2
  } else stop("Invalid number of channels")
  
  arm_l <- process_arm(ang_vel_deg[, 1:3], fs, TH_min_ampl)
  if (nr_imu == 2) {
    arm_r <- process_arm(ang_vel_deg[, 4:6], fs, TH_min_ampl)
  } else {
    arm_r <- NULL
  }
  
  if (!is.null(arm_r) && !all(is.na(arm_r$amplitude))) {
    perc_both <- min(arm_l$perc_time_swing, arm_r$perc_time_swing)
    ampl_asym <- (mean(arm_l$amplitude) - mean(arm_r$amplitude)) /
      max(mean(arm_l$amplitude), mean(arm_r$amplitude)) * 100
    pkvel_asym <- (mean(arm_l$pk_ang_vel) - mean(arm_r$pk_ang_vel)) /
      max(mean(arm_l$pk_ang_vel), mean(arm_r$pk_ang_vel)) * 100
    coord <- NA
  } else {
    perc_both <- NA
    ampl_asym <- NA
    pkvel_asym <- NA
    coord <- NA
  }
  
  list(
    arm_swing_l = arm_l,
    arm_swing_r = arm_r,
    perc_both_swing = perc_both,
    amplitude_asymmetry = ampl_asym,
    peak_velocity_asymmetry = pkvel_asym,
    coordination_max = coord
  )
}

# ----- Return Single Arm Metrics as Tibble -----
get_single_arm_metrics <- function(arm) {
  if (is.null(arm) || is.na(arm$amplitude[1])) {
    return(tibble(
      regularity_angle = NA,
      regularity_ang_vel = NA,
      pk_vel_forward = NA,
      pk_vel_backward = NA,
      perc_time_swing = NA,
      frequency = NA
    ))
  }
  
  tibble(
    regularity_angle = arm$regularity_angle,
    regularity_ang_vel = arm$regularity_ang_vel,
    pk_vel_forward = arm$pk_vel_forward,
    pk_vel_backward = arm$pk_vel_backward,
    perc_time_swing = arm$perc_time_swing,
    frequency = arm$frequency
  )
}

# ----- Return Bilateral Metrics as Tibble -----
get_bilateral_metrics <- function(result) {
  if (is.na(result$perc_both_swing)) {
    return(tibble(
      perc_both_swing = NA,
      amplitude_asymmetry = NA,
      peak_velocity_asymmetry = NA,
      coordination_max = NA
    ))
  }
  
  tibble(
    perc_both_swing = result$perc_both_swing,
    amplitude_asymmetry = result$amplitude_asymmetry,
    peak_velocity_asymmetry = result$peak_velocity_asymmetry,
    coordination_max = result$coordination_max
  )
}
