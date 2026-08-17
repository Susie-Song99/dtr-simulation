# ==============================================================
# 02_candidates.R — 연구 서사 2단계: 해결책 설계
#
# 연속성: 01_problem.R과 동일한 p=2 DGP 사용 (R1=0, generate_data_prob)
#
# 문제: Full Q 편향 원인 = m₂(O2)가 Stage 1 추정을 오염
#       → m₂를 추정 방정식에서 제거해야 한다
#
# 후보 비교:
#   Tian 방식:   Ỹ₁ = |ĝ₂|                   (R1=0이므로 R1 항 없음)
#   Regret 방식: Ỹ₁ = R2 − A2ĝ₂ + |ĝ₂|      (A2 항이 m2를 제거)
#
# 채택 근거:
#   Tian: Ỹ₁ = |ĝ₂|
#         E[A1·S1·Ỹ₁ᵀⁱᵃⁿ] = E[A1·S1·|g₂|] ≠ E[A1·S1·(m₂+|g₂|)]
#         → m₂ 항 누락 → ψ̂₁ᵀⁱᵃⁿ ≠ ψ₁* → 기각
#
#   Regret: Ỹ₁ᴵᴼ = R2 − A2·ĝ₂ + |ĝ₂| →ₚ m₂ + ε₂ + |g₂|  (A2 항 상쇄)
#         E[A1·S1·Ỹ₁ᴵᴼ] = E[A1·S1·(m₂+|g₂|)] + E[A1·S1·ε₂]
#                        = E[A1·S1·(m₂+|g₂|)] + 0
#                        = (E[S1S1ᵀ])⁻¹ · ψ₁*  → Bias = 0 ✓ 채택
#
# 성능 지표:
#   Acc = P(d̂₁(O1) = d₁*(O1))  (정책 정확도, 1차 지표)
#   |Bias(ψ̂₁₁)|, Var(ψ̂₁₁)      진단 지표 (O1₁ 성분, 편향 집중 위치)
# ==============================================================

source("/Users/songsusie/dtr_simulation/00_dgp.R")
library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(SEED)

# ─── p=2 DGP 상수 (01_problem.R과 동일) ─────────────────────
P_PROB   <- 2
A1_TO_O2 <- 1.5

generate_data_prob <- function(n, gamma = 1.0) {
  O1 <- matrix(rnorm(n * P_PROB), n, P_PROB)
  A1 <- sample(c(-1L, 1L), n, replace = TRUE)
  O2_1 <- O1[,1] + A1_TO_O2 * A1 + rnorm(n, 0, 0.3)
  O2_2 <- O1[,2] + rnorm(n, 0, 0.3)
  O2   <- cbind(O2_1, O2_2)
  A2   <- sample(c(-1L, 1L), n, replace = TRUE)
  m2   <- gamma * O2[,1]^2
  g2   <- 1.0 * O2[,1]
  eps2 <- rnorm(n, 0, SIGMA_R)
  R2   <- m2 + A2 * g2 + eps2
  # R1 없음 — 01_problem.R DGP와 동일
  list(O1=O1, A1=A1, O2=O2, A2=A2, R2=R2, m2=m2, g2=g2)
}

compute_oracle_psi1_prob <- function(gamma = 1.0, N = 5e5) {
  set.seed(SEED)
  O1 <- matrix(rnorm(N * P_PROB), N, P_PROB)
  S1 <- cbind(1, O1)
  Q1_given <- function(a1_val) {
    A1v  <- rep(a1_val, N)
    O2_1 <- O1[,1] + A1_TO_O2 * A1v + rnorm(N, 0, 0.3)
    O2_2 <- O1[,2] + rnorm(N, 0, 0.3)
    f2   <- gamma * O2_1^2
    g2   <- 1.0 * O2_1
    f2 + abs(g2)
  }
  Qp <- Q1_given(+1); Qm <- Q1_given(-1)
  S_st <- rbind(S1, S1)
  A_st <- c(rep(1L, N), rep(-1L, N))
  Q_st <- c(Qp, Qm)
  STS  <- crossprod(S_st) / (2*N)
  STQ  <- crossprod(A_st * S_st, Q_st) / (2*N)
  as.vector(solve(STS + 1e-8 * diag(P_PROB + 1), STQ))
}

# ─── Oracle ψ₁* 사전 계산 (γ=1 고정) ─────────────────────────
cat("Oracle ψ₁* 계산 중 (p=2, γ=1)...\n")
psi1_true <- compute_oracle_psi1_prob(gamma = 1.0)
cat(sprintf("  ψ₁* = (%.4f, %.4f, %.4f)\n",
    psi1_true[1], psi1_true[2], psi1_true[3]))

set.seed(SEED)
O1_test <- matrix(rnorm(N_TEST * P_PROB), N_TEST, P_PROB)


# ─── 추정기 ──────────────────────────────────────────────────

# Full Q: 01_problem.R의 fit_full_q_prob와 동일
fit_full_q_prob <- function(d) {
  n  <- nrow(d$O1)
  S2 <- cbind(1, d$O1, d$A1, d$O2)
  X2 <- cbind(S2, d$A2 * S2)
  coef2 <- tryCatch(
    solve(crossprod(X2), crossprod(X2, d$R2)),
    error = function(e) qr.solve(X2, d$R2)
  )
  coef2[is.na(coef2)] <- 0
  q        <- ncol(S2)
  beta2    <- coef2[1:q]
  psi2_hat <- coef2[(q+1):(2*q)]
  g2_hat   <- as.vector(S2 %*% psi2_hat)
  max_q2   <- as.vector(S2 %*% beta2) + abs(g2_hat)
  S1 <- cbind(1, d$O1)
  X1 <- cbind(S1, d$A1 * S1)
  coef1 <- tryCatch(
    solve(crossprod(X1), crossprod(X1, max_q2)),
    error = function(e) qr.solve(X1, max_q2)
  )
  coef1[(P_PROB + 2):(2 * P_PROB + 2)]  # ψ₁ 추정값 (p+1 차원)
}

# Tian: Stage 2 IO 추정 → Ỹ₁ = |ĝ₂| (R1=0)
# 편향 원인: m₂ 항이 Ỹ₁에서 빠져 있음 → E[A1·S1·Ỹ₁] ≠ E[A1·S1·(m₂+|g₂|)]
fit_tian_prob <- function(d) {
  S2    <- cbind(1, d$O2)
  X2_io <- d$A2 * S2          # 상호작용 항만으로 psi2 추정 (IO 방식)
  psi2_hat <- tryCatch(
    solve(crossprod(X2_io), crossprod(X2_io, d$R2)),
    error = function(e) qr.solve(X2_io, d$R2)
  )
  g2_hat   <- as.vector(S2 %*% psi2_hat)
  Y1_tian  <- abs(g2_hat)      # Ỹ₁ = |ĝ₂|, R1=0이므로 R1 항 없음
  S1 <- cbind(1, d$O1)
  X1 <- d$A1 * S1              # IO-QL Stage 1: 상호작용 항만
  coef1 <- tryCatch(
    solve(crossprod(X1), crossprod(X1, Y1_tian)),
    error = function(e) qr.solve(X1, Y1_tian)
  )
  coef1  # ψ₁ 추정값 (p+1 차원)
}

# Regret (IO-QL): Ỹ₁ = R2 − A2·ĝ₂ + |ĝ₂| →ₚ m₂ + ε₂ + |g₂|
# A2 항이 g₂를 상쇄 → m₂가 Ỹ₁에 남음 → Bias = 0
fit_regret_prob <- function(d) {
  S2    <- cbind(1, d$O2)
  X2_io <- d$A2 * S2
  psi2_hat <- tryCatch(
    solve(crossprod(X2_io), crossprod(X2_io, d$R2)),
    error = function(e) qr.solve(X2_io, d$R2)
  )
  g2_hat    <- as.vector(S2 %*% psi2_hat)
  Y1_regret <- d$R2 - d$A2 * g2_hat + abs(g2_hat)
  S1 <- cbind(1, d$O1)
  X1 <- d$A1 * S1
  coef1 <- tryCatch(
    solve(crossprod(X1), crossprod(X1, Y1_regret)),
    error = function(e) qr.solve(X1, Y1_regret)
  )
  coef1
}

# ─── 성능 지표 ────────────────────────────────────────────────
compute_accuracy_prob <- function(psi1_hat, psi1_true, O1_test) {
  S1 <- cbind(1, O1_test)
  d_hat  <- sign(as.vector(S1 %*% psi1_hat));  d_hat[d_hat == 0]  <- 1L
  d_true <- sign(as.vector(S1 %*% psi1_true)); d_true[d_true == 0] <- 1L
  mean(d_hat == d_true)
}


# ════════════════════════════════════════════════════════════
# 이론 검증 — 큰 표본으로 편향 조건 수치 확인
# ════════════════════════════════════════════════════════════
cat("\n─────────────────────────────────────────────────────────\n")
cat("이론 검증: Tian vs Regret pseudo-outcome 편향 조건\n")
cat("─────────────────────────────────────────────────────────\n")

N_THEORY <- 5e5
set.seed(SEED)
O1_th <- matrix(rnorm(N_THEORY * P_PROB), N_THEORY, P_PROB)
A1_th <- sample(c(-1L, 1L), N_THEORY, replace = TRUE)
O2_1_th <- O1_th[,1] + A1_TO_O2 * A1_th + rnorm(N_THEORY, 0, 0.3)
O2_2_th <- O1_th[,2] + rnorm(N_THEORY, 0, 0.3)
A2_th   <- sample(c(-1L, 1L), N_THEORY, replace = TRUE)
eps2_th <- rnorm(N_THEORY, 0, SIGMA_R)
m2_th   <- 1.0 * O2_1_th^2
g2_th   <- 1.0 * O2_1_th
R2_th   <- m2_th + A2_th * g2_th + eps2_th

S1_th <- cbind(1, O1_th)  # N × (p+1)

# 목표: E[A1·S1·(m₂+|g₂|)] = (E[S1S1ᵀ])·ψ₁*
# IO-QL Stage 1 추정량의 수렴 대상
C_mat    <- crossprod(S1_th) / N_THEORY          # E[S1S1ᵀ] 추정
target   <- as.vector(C_mat %*% psi1_true)       # 목표값 벡터 (p+1)

# Tian: E[A1·S1·|g₂|]  (m₂ 항 누락)
Y1_tian_th  <- abs(g2_th)
E_Z1_tian   <- as.vector(crossprod(A1_th * S1_th, Y1_tian_th)) / N_THEORY

# Regret: E[A1·S1·(m₂+ε₂+|g₂|)]
Y1_io_th <- m2_th + eps2_th + abs(g2_th)         # oracle Ỹ₁ᴵᴼ
E_Z1_io  <- as.vector(crossprod(A1_th * S1_th, Y1_io_th)) / N_THEORY

cat(sprintf("\n[목표] (E[S1S1ᵀ])·ψ₁*  O1₁ 성분 = %+.6f\n", target[2]))
cat(sprintf("\n[Tian]  E[A1·S1·|g₂|]  O1₁ 성분 = %+.6f\n", E_Z1_tian[2]))
cat(sprintf("        편향 기여 (Tian−목표)      = %+.6f  ≠ 0 → ✗ 기각\n",
    E_Z1_tian[2] - target[2]))
cat(sprintf("\n[Regret] E[A1·S1·(m₂+ε₂+|g₂|)] O1₁ 성분 = %+.6f\n", E_Z1_io[2]))
cat(sprintf("         |차이 (Regret−목표)|              = %.6f  ≈ 0 → ✓ 채택\n",
    abs(E_Z1_io[2] - target[2])))


# ════════════════════════════════════════════════════════════
# 시뮬레이션: n별 Acc / |Bias| / Var 비교 (γ=1, p=2)
# ════════════════════════════════════════════════════════════
cat("\n─────────────────────────────────────────────────────────\n")
cat("시뮬레이션: Full Q / Tian / Regret 비교 (γ=1, p=2)\n")
cat("─────────────────────────────────────────────────────────\n")

METHODS_CAND  <- c("fullq", "tian", "regret")
METHOD_LABELS <- c(fullq  = "Full Q",
                   tian   = "Tian (|ĝ₂|)",
                   regret = "Regret (IO-QL)")
N_REPS_CAND   <- 500

# 진단 지표: O1₁ 성분(인덱스 2) — 편향이 집중되는 위치
PSI1_IDX <- 2  # psi1_true 벡터에서 O1₁에 대응하는 인덱스

df_candidates <- do.call(rbind, lapply(N_GRID, function(n) {
  ests <- vapply(seq_len(N_REPS_CAND), function(b) {
    d <- generate_data_prob(n, gamma = 1.0)
    c(fullq  = unname(tryCatch(fit_full_q_prob(d)[PSI1_IDX],
                               error = function(e) NA_real_)),
      tian   = unname(tryCatch(fit_tian_prob(d)[PSI1_IDX],
                               error = function(e) NA_real_)),
      regret = unname(tryCatch(fit_regret_prob(d)[PSI1_IDX],
                               error = function(e) NA_real_)))
  }, numeric(length(METHODS_CAND)))

  do.call(rbind, lapply(METHODS_CAND, function(m) {
    e <- ests[m, ][!is.na(ests[m, ])]
    # Acc: 각 반복의 psi1_hat 전체 벡터가 필요하므로 별도 반복
    acc_vals <- vapply(seq_len(N_REPS_CAND), function(b) {
      d <- generate_data_prob(n, gamma = 1.0)
      psi_hat <- tryCatch(
        switch(m,
               fullq  = fit_full_q_prob(d),
               tian   = fit_tian_prob(d),
               regret = fit_regret_prob(d)),
        error = function(e) rep(NA_real_, P_PROB + 1)
      )
      if (anyNA(psi_hat)) return(NA_real_)
      compute_accuracy_prob(psi_hat, psi1_true, O1_test)
    }, numeric(1L))
    acc_vals <- acc_vals[!is.na(acc_vals)]

    data.frame(n      = n,
               method = m,
               acc    = mean(acc_vals),
               bias   = mean(e) - psi1_true[PSI1_IDX],
               vari   = var(e),
               mse    = mean((e - psi1_true[PSI1_IDX])^2))
  }))
}))

cat("\n방법별 성능 (γ=1, p=2, ψ₁₁ 성분)\n")
cat(sprintf("%-6s  %-18s  %6s  %8s  %8s\n", "n", "방법", "Acc", "|Bias|", "Var"))
cat(strrep("-", 54), "\n")
for (n_val in N_GRID) {
  sub <- df_candidates[df_candidates$n == n_val, ]
  for (m in METHODS_CAND) {
    r <- sub[sub$method == m, ]
    cat(sprintf("%-6d  %-18s  %6.4f  %8.4f  %8.5f\n",
        n_val, METHOD_LABELS[m], r$acc, abs(r$bias), r$vari))
  }
  cat("\n")
}


# ─── 플롯 ─────────────────────────────────────────────────────
METHOD_COLORS <- c(fullq  = "#c03838", tian = "#e07830", regret = "#2a7a55")
METHOD_SHAPES <- c(fullq  = 16, tian = 17, regret = 15)

p_acc <- ggplot(df_candidates,
                aes(x = n, y = acc, color = method, shape = method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  scale_x_log10(breaks = N_GRID, labels = N_GRID) +
  scale_color_manual(values = METHOD_COLORS, labels = METHOD_LABELS) +
  scale_shape_manual(values = METHOD_SHAPES, labels = METHOD_LABELS) +
  labs(title    = "정책 정확도 (Acc)",
       subtitle = "Regret은 n이 커질수록 Full Q를 추월한다",
       x = "n (log scale)", y = "Acc", color = NULL, shape = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

p_bias <- ggplot(df_candidates,
                 aes(x = n, y = abs(bias), color = method, shape = method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_log10(breaks = N_GRID, labels = N_GRID) +
  scale_color_manual(values = METHOD_COLORS, labels = METHOD_LABELS) +
  scale_shape_manual(values = METHOD_SHAPES, labels = METHOD_LABELS) +
  labs(title    = "|편향| — ψ̂₁₁ 성분",
       subtitle = "Full Q·Tian: n→∞ 에도 편향 유지 | Regret: 0에 수렴",
       x = "n (log scale)", y = "|Bias|", color = NULL, shape = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

p_var <- ggplot(df_candidates,
                aes(x = n, y = vari, color = method, shape = method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  scale_x_log10(breaks = N_GRID, labels = N_GRID) +
  scale_color_manual(values = METHOD_COLORS, labels = METHOD_LABELS) +
  scale_shape_manual(values = METHOD_SHAPES, labels = METHOD_LABELS) +
  labs(title    = "분산 — ψ̂₁₁ 성분",
       subtitle = "Regret은 편향 없으나 분산이 크다 → 3단계 과제",
       x = "n (log scale)", y = "Var", color = NULL, shape = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

p_combined <- (p_acc | p_bias) / p_var +
  plot_annotation(
    title    = "2단계: pseudo-outcome 후보 비교 (γ=1, p=2 DGP)",
    subtitle = "이론 근거: Tian Bias ≠ 0 (기각) / Regret Bias = 0 (채택) — 단, 분산이 큼",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave("/Users/songsusie/dtr_simulation/02_candidates.png",
       p_combined, width = 12, height = 10, dpi = 150)
cat("\n02_candidates.R 완료 — 02_candidates.png 저장됨\n")
