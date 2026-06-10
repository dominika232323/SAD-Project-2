set.seed(42) # Zapewnienie powtarzalności wyników

# ---  Parametry symulacji ---
N <- 1000          # Długość sekwencji (liczba próbek)
alpha <- 0.05      # Poziom istotności (prawdopodobieństwo fałszywego alarmu)
sigma <- 2.0       # Odchylenie standardowe szumu (wariancja = 4.0)
A <- 0.3           # Amplituda sygnału (przypadek krytyczny: A << sigma^2)
M <- 5000          # Liczba powtórzeń Monte Carlo dla każdej hipotezy

# ---  Wyznaczenie teoretycznego progu detekcji ---
# Pod H0 statystyka G ~ Chi-square(df=1)
gamma_thresh <- qchisq(1 - alpha, df = 1)
cat(paste("Teoretyczny próg detekcji (gamma) dla alpha =", alpha, "wynosi:", round(gamma_thresh, 4), "\n"))

# ---  Generowanie kodu binarnego s(i) ---
# Kod s(i) ze zbioru {-1, 1}
s <- sample(c(-1, 1), N, replace = TRUE)

# ---  Funkcja obliczająca statystykę GLRT G(x) ---
compute_glrt <- function(x, s_code) {
  N_samples <- length(x)
  sum_x2 <- sum(x^2)
  sum_xs <- sum(x * s_code)
  
  # Mianownik ułamka pod logarytmem (S_0 - A_hat^2 * N)
  denominator <- sum_x2 - (1 / N_samples) * (sum_xs)^2
  
  # Statystyka G(x)
  G <- N_samples * log(sum_x2 / denominator)
  return(G)
}

# ---  Pętla Monte Carlo ---
G_H0 <- numeric(M)
G_H1 <- numeric(M)

for (m in 1:M) {
  # Generowanie czystego szumu gaussowskiego
  w <- rnorm(N, mean = 0, sd = sigma)
  
  # Hipoteza H0: Tylko szum
  x_H0 <- w
  G_H0[m] <- compute_glrt(x_H0, s)
  
  # Hipoteza H1: Sygnał + Szum
  x_H1 <- A * s + w
  G_H1[m] <- compute_glrt(x_H1, s)
}

# ---  Analiza wyników i empiryczna weryfikacja ---
empiryczny_alpha <- mean(G_H0 > gamma_thresh)
empiryczna_moc   <- mean(G_H1 > gamma_thresh)

# Teoretyczna moc testu przy użyciu niecentralnego chi-kwadrat
lambda_param <- N * (A^2) / (sigma^2)
teoretyczna_moc <- 1 - pchisq(gamma_thresh, df = 1, ncp = lambda_param)

# ---  Wyświetlenie raportu ---
cat("\n=== WYNIKI SYMULACJI ===\n")
cat(sprintf("Zadany poziom istotności (alpha):     %.4f\n", alpha))
cat(sprintf("Empiryczny fałszywy alarm (H0):       %.4f\n", empiryczny_alpha))
cat("----------------------------------------\n")
cat(sprintf("Parametr niecentralności (lambda):    %.4f\n", lambda_param))
cat(sprintf("Teoretyczna moc testu (P_D):          %.4f\n", teoretyczna_moc))
cat(sprintf("Empiryczna moc testu (z symulacji):   %.4f\n", empiryczna_moc))

# ---  Wizualizacja gęstości rozkładów ---
plot(density(G_H0), col = "blue", lwd = 2, xlim = c(0, max(G_H1)),
     main = "Rozkład statystyki GLRT pod hipotezą H0 i H1",
     xlab = "Wartość statystyki G(x)", ylab = "Gęstość prawdopodobieństwa")
lines(density(G_H1), col = "red", lwd = 2)
abline(v = gamma_thresh, col = "darkgreen", lty = 2, lwd = 2)
legend("topright", legend = c("H0 (Tylko szum)", "H1 (Sygnał+Szum)", "Próg detekcji (gamma)"),
       col = c("blue", "red", "darkgreen"), lty = c(1, 1, 2), lwd = 2)


# ----------------------------------------------------------------------
# Analiza teoretycznej mocy detektora GLRT w funkcji A oraz N
# ----------------------------------------------------------------------

# Czyszczenie przestrzeni i ustawienie parametrów stałych
graphics.off()
alpha <- 0.05
sigma <- 2.0
gamma_thresh <- qchisq(1 - alpha, df = 1)

# Ustawienie podziału okna graficznego na 2 wykresy obok siebie
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# ======================================================================
# WYKRES 1: Moc testu vs Amplituda (A) dla różnych N
# ======================================================================
A_seq <- seq(0, 0.6, length.out = 200) # Zakres amplitud do analizy
N_values <- c(250, 500, 1000, 2000)    # Różne długości sekwencji
kolory_N <- c("orange", "green3", "blue", "purple")

# Inicjalizacja pustego wykresu
plot(NULL, xlim = c(0, max(A_seq)), ylim = c(0, 1),
     xlab = "Amplituda sygnału (A)", 
     ylab = "Prawdopodobieństwo wykrycia (P_D)",
     main = "Moc testu w funkcji amplitudy",
     cex.lab = 1.2, cex.main = 1.1)
abline(h = alpha, col = "gray", lty = 3) # Poziom fałszywego alarmu
abline(h = 1.0, col = "gray", lty = 3)

# Rysowanie krzywych dla każdego N
for (i in 1:length(N_values)) {
  N_curr <- N_values[i]
  # Obliczenie parametru niecentralności lambda dla wektora amplitud
  lambda_seq <- N_curr * (A_seq^2) / (sigma^2)
  # Teoretyczna moc testu
  Pow_seq <- 1 - pchisq(gamma_thresh, df = 1, ncp = lambda_seq)
  
  lines(A_seq, Pow_seq, col = kolory_N[i], lwd = 2.5)
}

legend("bottomright", legend = paste("N =", N_values),
       col = kolory_N, lty = 1, lwd = 2.5, bty = "n", title = "Długość kodu")

# ======================================================================
# WYKRES 2: Moc testu vs Liczba próbek (N) dla różnych A
# ======================================================================
N_seq <- seq(10, 3000, length.out = 300) # Zakres liczby próbek do analizy
A_values <- c(0.1, 0.2, 0.3, 0.4)        # Różne wartości amplitudy
kolory_A <- c("red", "chocolate", "darkblue", "magenta4")

# Inicjalizacja pustego wykresu
plot(NULL, xlim = c(0, max(N_seq)), ylim = c(0, 1),
     xlab = "Liczba próbek (N)", 
     ylab = "Prawdopodobieństwo wykrycia (P_D)",
     main = "Moc testu w funkcji liczby próbek",
     cex.lab = 1.2, cex.main = 1.1)
abline(h = alpha, col = "gray", lty = 3)
abline(h = 1.0, col = "gray", lty = 3)

# Rysowanie krzywych dla każdej amplitudy A
for (i in 1:length(A_values)) {
  A_curr <- A_values[i]
  # Obliczenie parametru niecentralności lambda dla wektora N
  lambda_seq <- N_seq * (A_curr^2) / (sigma^2)
  # Teoretyczna moc testu
  Pow_seq <- 1 - pchisq(gamma_thresh, df = 1, ncp = lambda_seq)
  
  lines(N_seq, Pow_seq, col = kolory_A[i], lwd = 2.5)
}

legend("bottomright", legend = paste("A =", A_values),
       col = kolory_A, lty = 1, lwd = 2.5, bty = "n", title = "Amplituda")

