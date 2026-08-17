# ==============================================================
# 03_variance.R — 연구 서사 3단계: 새로운 문제 발견
#
# 연속성: 01~02와 동일한 p=2 DGP (R1=0, generate_data_prob)
#
# 관찰: Regret(IO-QL)은 Bias=0이지만 분산이 Full Q보다 크다
#
# 분석:
#   3A. ψ̂₁ 추정 분포: 편향 vs 분산 trade-off 시각화
#   3B. Acc vs n: 두 방법 수렴 행태 비교
#   3C. 분산 원천 분해:
#       R1=0이므로 ε₁·m̃₁ 없음 → 구조가 더 단순하다
#
#       IO-QL η₁ = M̃ + ε₂ + G̃
#         M̃ = m₂ − E[m₂]             (Stage 2 주효과 전체 변동)
#         ε₂                          (Stage 2 잡음)
#         G̃ = |g₂| − E[|g₂|]         (블립 변동, 비조건부)
#
#       Full Q η₁ ≈ L + G
#         L  = m₂_linear − E[m₂_linear|O1,A1]  (m₂ 선형근사 변동)
#         G  = |g₂| − E[|g₂||O1,A1]            (블립 변동, 조건부)
#
#       차이: IO-QL이 M̃ 전체를, Full Q는 선형근사 L만 노출
#             IO-QL이 G̃(비조건부), Full Q는 G(조건부) → G̃ ≥ G
#             IO-QL이 ε₂ 추가 노출
# ==============================================================

source("/Users/songsusie/dtr_simulation/00_dgp.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

set.seed(SEED)

# ─── p=2 DGP 상수 (01~02와 동일) ─────────────────────────────
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
    gamma * O2_1^2 + abs(O2_1)
  }
  Qp <- Q1_given(+1); Qm <- Q1_given(-1)
  S_st <- rbind(S1, S1)
  A_st <- c(rep(1L, N), rep(-1L, N))
  STS  <- crossprod(S_st) / (2*N)
  STQ  <- crossprod(A_st * S_st, c(Qp, Qm)) / (2*N)
  as.vector(solve(STS + 1e-8 * diag(P_PROB + 1), STQ))
}

# ─── Oracle ψ₁* (γ=1) ──────────────────────────────────────
cat("Oracle ψ₁* 계산 중 (p=2, γ=1)...\n")
psi1_true <- compute_oracle_psi1_prob(gamma = 1.0)
cat(sprintf("  ψ₁* = (%.4f, %.4f, %.4f)\n",
    psi1_true[1], psi1_true[2], psi1_true[3]))

set.seed(SEED)
O1_test <- matrix(rnorm(N_TEST * P_PROB), N_TEST, P_PROB)

PSI1_IDX <- 2  # O1₁ 성분 — 편향·분산이 집중되는 위치

# ─── 추정기 (02와 동일) ───────────────────────────────────────
fit_full_q_prob <- function(d) {
  S2 <- cbind(1, d$O1, d$A1, d$O2)
  X2 <- cbind(S2, d$A2 * S2)
  coef2 <- tryCatch(
    solve(crossprod(X2), crossprod(X2, d$R2)),
    error = function(e) qr.solve(X2, d$R2)
  )
  coef2[is.na(coef2)] <- 0
  q        <- ncol(S2)
  beta2    <- coef2[1:q]; psi2_hat <- coef2[(q+1):(2*q)]
  g2_hat   <- as.vector(S2 %*% psi2_hat)
  max_q2   <- as.vector(S2 %*% beta2) + abs(g2_hat)
  S1 <- cbind(1, d$O1)
  X1 <- cbind(S1, d$A1 * S1)
  coef1 <- tryCatch(
    solve(crossprod(X1), crossprod(X1, max_q2)),
    error = function(e) qr.solve(X1, max_q2)
  )
  coef1[(P_PROB + 2):(2 * P_PROB + 2)]
}

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

compute_accuracy_prob <- function(psi1_hat, psi1_true, O1_test) {
  S1 <- cbind(1, O1_test)
  d_hat  <- sign(as.vector(S1 %*% psi1_hat));  d_hat[d_hat == 0]  <- 1L
  d_true <- sign(as.vector(S1 %*% psi1_true)); d_true[d_true == 0] <- 1L
  mean(d_hat == d_true)
}


# ════════════════════════════════════════════════════════════
# Phase 3A — ψ̂₁ 추정 분포 시각화 (n=500, γ=1)
# Full Q: 편향으로 분포 이동 / IO-QL: 중심은 ψ₁*이나 분산 큼
# ════════════════════════════════════════════════════════════
cat("\n─────────────────────────────────────────\n")
cat("Phase 3A: ψ̂₁ 추정 분포 시각화\n")
cat("─────────────────────────────────────────\n")

N_REPS_DIST <- 1000
N_DIST      <- 500

psi1_fullq_dist  <- vapply(seq_len(N_REPS_DIST), function(b) {
  tryCatch(fit_full_q_prob(generate_data_prob(N_DIST, 1.0))[PSI1_IDX],
           error = function(e) NA_real_)
}, numeric(1L))
psi1_regret_dist <- vapply(seq_len(N_REPS_DIST), function(b) {
  tryCatch(fit_regret_prob(generate_data_prob(N_DIST, 1.0))[PSI1_IDX],
           error = function(e) NA_real_)
}, numeric(1L))
psi1_fullq_dist  <- psi1_fullq_dist[!is.na(psi1_fullq_dist)]
psi1_regret_dist <- psi1_regret_dist[!is.na(psi1_regret_dist)]

cat(sprintf("\nFull Q:  평균=%.4f, SD=%.4f, |Bias|=%.4f\n",
    mean(psi1_fullq_dist), sd(psi1_fullq_dist),
    abs(mean(psi1_fullq_dist) - psi1_true[PSI1_IDX])))
cat(sprintf("IO-QL:   평균=%.4f, SD=%.4f, |Bias|=%.4f\n",
    mean(psi1_regret_dist), sd(psi1_regret_dist),
    abs(mean(psi1_regret_dist) - psi1_true[PSI1_IDX])))
cat(sprintf("ψ₁* (O1₁ 성분) = %.4f\n", psi1_true[PSI1_IDX]))

df_dist <- data.frame(
  psi1   = c(psi1_fullq_dist, psi1_regret_dist),
  method = rep(c("Full Q", "IO-QL Regret"),
               c(length(psi1_fullq_dist), length(psi1_regret_dist)))
)

p_dist <- ggplot(df_dist, aes(x = psi1, fill = method, color = method)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_vline(xintercept = psi1_true[PSI1_IDX],
             linetype = "solid", color = "black", linewidth = 1.0) +
  annotate("text", x = psi1_true[PSI1_IDX] + 0.1, y = Inf,
           label = "ψ₁*", vjust = 1.5, size = 4, fontface = "bold") +
  geom_vline(xintercept = mean(psi1_fullq_dist),
             linetype = "dashed", color = "#c03838", linewidth = 0.7) +
  geom_vline(xintercept = mean(psi1_regret_dist),
             linetype = "dashed", color = "#2a7a55", linewidth = 0.7) +
  scale_fill_manual(values  = c("Full Q" = "#c03838", "IO-QL Regret" = "#2a7a55")) +
  scale_color_manual(values = c("Full Q" = "#c03838", "IO-QL Regret" = "#2a7a55")) +
  labs(title    = "Phase 3A: ψ̂₁₁ 추정 분포 (n=500, γ=1, B=1000)",
       subtitle = "Full Q: 편향으로 분포 이동 | IO-QL: 중심은 ψ₁*이나 분산이 큼",
       x = "ψ̂₁₁ (O1₁ 성분)", y = "밀도", fill = NULL, color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))


# ════════════════════════════════════════════════════════════
# Phase 3B — Acc vs n: 두 방법 수렴 행태 비교
# ════════════════════════════════════════════════════════════
cat("\n─────────────────────────────────────────\n")
cat("Phase 3B: Acc vs n — 수렴 행태 비교\n")
cat("─────────────────────────────────────────\n")

N_REPS_ACC <- 500

df_acc_n <- do.call(rbind, lapply(N_GRID, function(n) {
  acc_vals <- vapply(seq_len(N_REPS_ACC), function(b) {
    d <- generate_data_prob(n, gamma = 1.0)
    psi_fq <- tryCatch(fit_full_q_prob(d),  error = function(e) rep(NA, P_PROB+1))
    psi_io <- tryCatch(fit_regret_prob(d),  error = function(e) rep(NA, P_PROB+1))
    c(fullq  = if (anyNA(psi_fq)) NA_real_ else
                 compute_accuracy_prob(psi_fq, psi1_true, O1_test),
      regret = if (anyNA(psi_io)) NA_real_ else
                 compute_accuracy_prob(psi_io, psi1_true, O1_test))
  }, numeric(2L))

  psi_ests <- vapply(seq_len(N_REPS_ACC), function(b) {
    d <- generate_data_prob(n, gamma = 1.0)
    c(fullq  = tryCatch(fit_full_q_prob(d)[PSI1_IDX],  error = function(e) NA_real_),
      regret = tryCatch(fit_regret_prob(d)[PSI1_IDX],  error = function(e) NA_real_))
  }, numeric(2L))

  do.call(rbind, lapply(c("fullq","regret"), function(m) {
    e   <- psi_ests[m, ][!is.na(psi_ests[m, ])]
    acc <- mean(acc_vals[m, ][!is.na(acc_vals[m, ])])
    data.frame(n = n, method = m, acc = acc,
               bias = mean(e) - psi1_true[PSI1_IDX], vari = var(e))
  }))
}))

cat("\n[Acc vs n, γ=1, p=2]\n")
cat(sprintf("%-6s  %-14s  %6s  %8s  %8s\n", "n", "방법", "Acc", "|Bias|", "Var"))
cat(strrep("-", 50), "\n")
for (n_val in N_GRID) {
  sub <- df_acc_n[df_acc_n$n == n_val, ]
  for (m in c("fullq","regret")) {
    r <- sub[sub$method == m, ]
    label <- ifelse(m=="fullq", "Full Q", "IO-QL Regret")
    cat(sprintf("%-6d  %-14s  %6.4f  %8.4f  %8.5f\n",
        n_val, label, r$acc, abs(r$bias), r$vari))
  }
  cat("\n")
}

METHOD_COLORS2 <- c(fullq = "#c03838", regret = "#2a7a55")
METHOD_LABELS2 <- c(fullq = "Full Q", regret = "IO-QL Regret")

p_acc_n <- ggplot(df_acc_n, aes(x = n, y = acc, color = method, shape = method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  scale_x_log10(breaks = N_GRID, labels = N_GRID) +
  scale_color_manual(values = METHOD_COLORS2, labels = METHOD_LABELS2) +
  scale_shape_manual(values = c(fullq=16, regret=15), labels = METHOD_LABELS2) +
  labs(title    = "Phase 3B: Acc vs n (γ=1, p=2)",
       subtitle = "Full Q: 편향으로 수렴 안 됨 | IO-QL: 수렴하지만 n 작을 때 불안정",
       x = "n (log scale)", y = "Acc", color = NULL, shape = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))


# ════════════════════════════════════════════════════════════
# Phase 3C — 분산 원천 분해 (p=2, R1=0)
#
# IO-QL η₁ = M̃ + ε₂ + G̃  (R1=0이라 ε₁·m̃₁ 없음)
# Full Q η₁ ≈ L + G
# ════════════════════════════════════════════════════════════
cat("\n─────────────────────────────────────────\n")
cat("Phase 3C: 분산 원천 분해 (p=2, R1=0)\n")
cat("─────────────────────────────────────────\n")

N_DECOMP <- 1e5

decompose_variance_prob <- function(gamma) {
  set.seed(SEED)
  O1   <- matrix(rnorm(N_DECOMP * P_PROB), N_DECOMP, P_PROB)
  A1   <- sample(c(-1L, 1L), N_DECOMP, replace = TRUE)
  O2_1 <- O1[,1] + A1_TO_O2 * A1 + rnorm(N_DECOMP, 0, 0.3)
  O2_2 <- O1[,2] + rnorm(N_DECOMP, 0, 0.3)
  A2   <- sample(c(-1L, 1L), N_DECOMP, replace = TRUE)
  eps2 <- rnorm(N_DECOMP, 0, SIGMA_R)
  m2   <- gamma * O2_1^2
  g2   <- 1.0 * O2_1
  g2_abs <- abs(g2)

  # ── IO-QL Regret η₁ 원천 ──
  # Ỹ₁ᴵᴼ = m2 + ε2 + |g2|  (oracle)
  M_tilde <- m2 - mean(m2)            # Stage 2 주효과 변동 (전체)
  G_tilde <- g2_abs - mean(g2_abs)    # 블립 변동 (비조건부)

  # ── Full Q η₁ 원천 ──
  # m2 선형근사: E[m2 | (1, O2₁)] 의 예측값 변동
  X_m2   <- cbind(1, O2_1)
  b_m2   <- solve(crossprod(X_m2), crossprod(X_m2, m2))
  m2_lin <- as.vector(X_m2 %*% b_m2)
  L_term <- m2_lin - mean(m2_lin)     # 선형근사 변동

  # |g2| 조건부 변동: E[|g2| | O1, A1] 회귀 잔차
  X_cond  <- cbind(1, O1, A1, O1[,1]*A1, O1[,2]*A1)
  b_cond  <- solve(crossprod(X_cond), crossprod(X_cond, g2_abs))
  G_cond  <- g2_abs - as.vector(X_cond %*% b_cond)  # 조건부 잔차

  c(
    var_eta_io   = var(M_tilde + eps2 + G_tilde),
    var_eta_full = var(L_term + G_cond),
    var_eps2     = var(eps2),
    var_M_tilde  = var(M_tilde),
    var_G_tilde  = var(G_tilde),
    var_L        = var(L_term),
    var_G_cond   = var(G_cond)
  )
}

cat("\n분산 원천 분해 (n=1e5, oracle):\n\n")
cat(sprintf("%-24s  %8s  %8s  %8s  %8s  %8s\n",
    "원천", "γ=0", "γ=0.5", "γ=1.0", "γ=1.5", "γ=2.0"))
cat(strrep("-", 72), "\n")

decomp_list <- lapply(GAMMA_GRID, decompose_variance_prob)
decomp_mat  <- do.call(cbind, decomp_list)
colnames(decomp_mat) <- as.character(GAMMA_GRID)

rows_print <- c(
  "Var(η₁, IO-QL)"    = "var_eta_io",
  "Var(η₁, Full Q)"   = "var_eta_full",
  "  ε₂ (IO만)"       = "var_eps2",
  "  M̃=m₂-E[m₂] (IO)" = "var_M_tilde",
  "  G̃=|g₂|-E[|g₂|] (IO)" = "var_G_tilde",
  "  L=m₂_lin (FQ)"   = "var_L",
  "  G=|g₂|_cond (FQ)"= "var_G_cond"
)
for (label in names(rows_print)) {
  vals <- decomp_mat[rows_print[label], ]
  cat(sprintf("%-24s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n",
      label, vals[1], vals[2], vals[3], vals[4], vals[5]))
}
cat(strrep("-", 72), "\n")
diff_vals <- decomp_mat["var_eta_io",] - decomp_mat["var_eta_full",]
cat(sprintf("%-24s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n",
    "Var 차이 (IO−FQ)",
    diff_vals[1], diff_vals[2], diff_vals[3], diff_vals[4], diff_vals[5]))

# ─── 분산 원천 누적 바 차트 ──────────────────────────────────
SOURCE_LABELS <- c(
  eps2    = "ε₂ (Stage 2 잡음, IO만)",
  M_tilde = "M̃ (m₂ 전체 변동, IO)",
  G_tilde = "G̃ (블립 비조건부, IO)",
  L       = "L (m₂ 선형근사, Full Q)",
  G_cond  = "G (블립 조건부, Full Q)"
)
SOURCE_COLORS <- c(
  eps2    = "#f0b060",
  M_tilde = "#c07840",
  G_tilde = "#8060c0",
  L       = "#aaaaaa",
  G_cond  = "#606060"
)

df_decomp <- do.call(rbind, lapply(seq_along(GAMMA_GRID), function(gi) {
  g <- GAMMA_GRID[gi]
  d <- decomp_list[[gi]]
  rbind(
    data.frame(gamma = g, method = "IO-QL Regret",
               source = c("eps2","M_tilde","G_tilde"),
               value  = c(d["var_eps2"], d["var_M_tilde"], d["var_G_tilde"])),
    data.frame(gamma = g, method = "Full Q",
               source = c("L","G_cond"),
               value  = c(d["var_L"], d["var_G_cond"]))
  )
}))
df_decomp$source <- factor(df_decomp$source, levels = names(SOURCE_LABELS))
df_decomp$method <- factor(df_decomp$method, levels = c("Full Q","IO-QL Regret"))

p_decomp <- ggplot(df_decomp, aes(x = factor(gamma), y = value, fill = source)) +
  geom_col(alpha = 0.85) +
  facet_wrap(~method) +
  scale_fill_manual(values = SOURCE_COLORS, labels = SOURCE_LABELS, name = "분산 원천") +
  labs(title    = "Phase 3C: 분산 원천 분해 (p=2, R1=0)",
       subtitle = "IO-QL: ε₂ + M̃ + G̃ 세 원천 노출 | Full Q: L + G만 (단, 편향 있음)",
       x = "오설정 강도 γ", y = "Var(η₁) 기여량") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"),
        legend.text = element_text(size = 8))


# ════════════════════════════════════════════════════════════
# Phase 3D — IO-QL 오류 위치 분석
# Full Q(1C)와 동일 방식으로: 어디서 틀리는가?
# n=500 (분산 문제가 드러나는 중간 표본) 고정
# ════════════════════════════════════════════════════════════
cat("\n─────────────────────────────────────────\n")
cat("Phase 3D: IO-QL 오류 위치 분석 (n=500, γ=1)\n")
cat("─────────────────────────────────────────\n")

N_DIAG_IO <- 5000
set.seed(SEED)
d_diag_io <- generate_data_prob(N_DIAG_IO, gamma = 1.0)

# IO-QL 추정 (n=500 단일 표본으로)
set.seed(SEED + 1)
d_train_io    <- generate_data_prob(500, gamma = 1.0)
psi1_hat_io   <- tryCatch(
  fit_regret_prob(d_train_io),
  error = function(e) rep(NA_real_, P_PROB + 1)
)

S1_diag_io <- cbind(1, d_diag_io$O1)
d_true_io  <- sign(as.vector(S1_diag_io %*% psi1_true))
d_hat_io   <- sign(as.vector(S1_diag_io %*% psi1_hat_io))
d_true_io[d_true_io == 0] <- 1L
d_hat_io[d_hat_io == 0]   <- 1L
error_ind_io <- (d_hat_io != d_true_io)

# 결정 경계까지의 거리 (oracle ψ₁* 기준)
margin_io <- abs(as.vector(S1_diag_io %*% psi1_true)) /
             sqrt(sum(psi1_true[-1]^2))

cat(sprintf("\nIO-QL 전체 오류율: %.4f (n_diag=%d, train n=500)\n",
    mean(error_ind_io), N_DIAG_IO))
cat(sprintf("경계 근처 절반 오류율: %.4f  /  먼 절반: %.4f\n",
    mean(error_ind_io[margin_io <  median(margin_io)]),
    mean(error_ind_io[margin_io >= median(margin_io)])))

# 10분위 오류율
df_margin_io <- data.frame(
  margin = margin_io,
  error  = as.integer(error_ind_io)
) %>%
  mutate(bin = dplyr::ntile(margin, 10)) %>%
  group_by(bin) %>%
  summarize(error_rate = mean(error), margin_mid = mean(margin), .groups = "drop")

cat(sprintf("10분위 가장 가까운 빈 오류율: %.4f  /  가장 먼 빈: %.4f\n",
    df_margin_io$error_rate[which.min(df_margin_io$margin_mid)],
    df_margin_io$error_rate[which.max(df_margin_io$margin_mid)]))

# Full Q 1C와 비교 (같은 N_DIAG로)
set.seed(SEED + 1)
d_train_fq  <- generate_data_prob(500, gamma = 1.0)
psi1_hat_fq <- tryCatch(
  fit_full_q_prob(d_train_fq),
  error = function(e) rep(NA_real_, P_PROB + 1)
)
d_hat_fq <- sign(as.vector(S1_diag_io %*% psi1_hat_fq))
d_hat_fq[d_hat_fq == 0] <- 1L
error_ind_fq <- (d_hat_fq != d_true_io)

df_margin_fq <- data.frame(
  margin = margin_io,
  error  = as.integer(error_ind_fq)
) %>%
  mutate(bin = dplyr::ntile(margin, 10)) %>%
  group_by(bin) %>%
  summarize(error_rate = mean(error), margin_mid = mean(margin), .groups = "drop")

cat(sprintf("Full Q  전체 오류율: %.4f\n", mean(error_ind_fq)))
cat(sprintf("Full Q  10분위 가장 가까운: %.4f  /  가장 먼: %.4f\n",
    df_margin_fq$error_rate[which.min(df_margin_fq$margin_mid)],
    df_margin_fq$error_rate[which.max(df_margin_fq$margin_mid)]))

# scatter: IO-QL 오류 위치
df_diag_io <- data.frame(
  O1_1   = d_diag_io$O1[,1],
  O1_2   = d_diag_io$O1[,2],
  error  = factor(error_ind_io, labels = c("정확","오류")),
  margin = margin_io
)
boundary_io <- data.frame(x = seq(-3, 3, length.out = 100))
boundary_io$y <- -(psi1_true[1] + psi1_true[2] * boundary_io$x) / psi1_true[3]

p_scatter_io <- ggplot(df_diag_io %>% filter(abs(O1_1) < 3, abs(O1_2) < 3)) +
  geom_point(aes(x=O1_1, y=O1_2, color=error, alpha=error), size=0.8) +
  geom_line(data=boundary_io, aes(x=x, y=y),
            color="black", linewidth=1.2, linetype="solid") +
  annotate("text", x=2.2, y=1.5, label="결정 경계", size=3, color="gray30") +
  scale_color_manual(values = c("정확"="#a8c8e8", "오류"="#2a7a55")) +
  scale_alpha_manual(values = c("정확"=0.3, "오류"=0.9)) +
  coord_cartesian(xlim=c(-3,3), ylim=c(-3,3)) +
  labs(title    = "Phase 3D-1: IO-QL 오류 위치 (n=500)",
       subtitle = "오류(초록)가 어디에 분포하는가",
       x="O1₁", y="O1₂", color=NULL, alpha=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="none", plot.title=element_text(face="bold"))

# margin 라인: IO-QL vs Full Q 오류율 곡선 비교
df_margin_both <- rbind(
  df_margin_io %>% mutate(method = "IO-QL (n=500)"),
  df_margin_fq %>% mutate(method = "Full Q (n=500)")
)

p_margin_compare <- ggplot(df_margin_both,
    aes(x=margin_mid, y=error_rate, color=method, linetype=method)) +
  geom_line(linewidth=1.0) + geom_point(size=2.5) +
  scale_color_manual(values=c("IO-QL (n=500)"="#2a7a55", "Full Q (n=500)"="#c03838")) +
  labs(title    = "Phase 3D-2: 결정경계 거리 vs 오류율 비교",
       subtitle = "두 방법의 오류 위치 패턴이 같은가 다른가",
       x="결정 경계까지 거리", y="오류율",
       color=NULL, linetype=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", plot.title=element_text(face="bold"))

p_3d <- p_scatter_io | p_margin_compare


# ════════════════════════════════════════════════════════════
# Phase 3E — 분산 원천 ablation: M̃이 경계 불안정의 원인인가?
#
# 각 분산 원천을 oracle로 제거했을 때
# Var(ψ̂₁₁)과 경계 오류율이 얼마나 줄어드는지 직접 측정.
#
# 버전 정의:
#   [Full]     Ỹ₁ = R2 - A2·ĝ2 + |ĝ2|           (실제 추정)
#   [Oracle]   Ỹ₁ = m2 + ε2 + |g2|               (M̃ 구조 유지, 추정 오차만 제거)
#   [No-M̃]    Ỹ₁ = ε2 + |g2|                    (m2 제거 → M̃ 완전 소거)
#   [No-ε₂]   Ỹ₁ = m2 + |g2|                    (ε2 제거)
#   [No-G̃]    Ỹ₁ = m2 + ε2                      (|g2| 제거)
#
# ε2 = R2 - m2 - A2·g2  (저장된 변수에서 역산)
# ════════════════════════════════════════════════════════════
cat("\n─────────────────────────────────────────\n")
cat("Phase 3E: 분산 원천 ablation — M̃이 경계 불안정 원인인가?\n")
cat("─────────────────────────────────────────\n")

N_ABL      <- 500   # 훈련 표본 크기
N_REPS_ABL <- 1000  # 반복 수
N_DIAG_ABL <- 5000  # 오류율 측정용 진단 표본

ABLATION_VERSIONS <- c("Full", "Oracle", "No-M", "No-eps2", "No-G")

# 훈련 데이터로 ψ̂₁₁ 추정, 진단 데이터로 경계 오류율 측정
set.seed(SEED)
O1_abl_diag <- matrix(rnorm(N_DIAG_ABL * P_PROB), N_DIAG_ABL, P_PROB)
S1_abl_diag <- cbind(1, O1_abl_diag)
d_true_abl  <- sign(as.vector(S1_abl_diag %*% psi1_true))
d_true_abl[d_true_abl == 0] <- 1L
margin_abl  <- abs(as.vector(S1_abl_diag %*% psi1_true)) /
               sqrt(sum(psi1_true[-1]^2))
near_idx    <- margin_abl < median(margin_abl)  # 경계 근처 절반

abl_results <- vapply(seq_len(N_REPS_ABL), function(b) {
  d   <- generate_data_prob(N_ABL, gamma = 1.0)
  S2  <- cbind(1, d$O2)
  X2_io <- d$A2 * S2
  psi2_hat <- tryCatch(
    solve(crossprod(X2_io), crossprod(X2_io, d$R2)),
    error = function(e) qr.solve(X2_io, d$R2)
  )
  g2_hat <- as.vector(S2 %*% psi2_hat)

  # oracle 변수: ε2 = R2 - m2 - A2·g2 (저장된 변수 역산)
  eps2_true <- d$R2 - d$m2 - d$A2 * d$g2
  g2_true   <- d$g2

  Y1_list <- list(
    "Full"    = d$R2 - d$A2 * g2_hat + abs(g2_hat),     # 실제 IO-QL
    "Oracle"  = d$m2 + eps2_true + abs(g2_true),         # 추정 오차만 제거
    "No-M"    = eps2_true + abs(g2_true),                 # M̃ 완전 소거
    "No-eps2" = d$m2 + abs(g2_true),                     # ε2 제거
    "No-G"    = d$m2 + eps2_true                         # G̃ 제거
  )

  S1_tr <- cbind(1, d$O1)
  X1_tr <- d$A1 * S1_tr

  out <- vapply(ABLATION_VERSIONS, function(v) {
    Y1 <- Y1_list[[v]]
    coef1 <- tryCatch(
      solve(crossprod(X1_tr), crossprod(X1_tr, Y1)),
      error = function(e) rep(NA_real_, P_PROB + 1)
    )
    if (anyNA(coef1)) return(NA_real_)
    coef1[PSI1_IDX]
  }, numeric(1L))

  out
}, numeric(length(ABLATION_VERSIONS)))
# abl_results: [버전 × 반복] 행렬

# Var(ψ̂₁₁) 및 경계 오류율 집계
abl_summary <- do.call(rbind, lapply(ABLATION_VERSIONS, function(v) {
  est <- abl_results[v, ][!is.na(abl_results[v, ])]

  # 경계 오류율: 각 반복의 ψ̂₁ 전체 벡터 필요 → 별도 소규모 반복
  set.seed(SEED)
  near_err_rates <- vapply(seq_len(200), function(b) {
    d   <- generate_data_prob(N_ABL, gamma = 1.0)
    S2  <- cbind(1, d$O2)
    X2_io <- d$A2 * S2
    psi2_hat <- tryCatch(
      solve(crossprod(X2_io), crossprod(X2_io, d$R2)),
      error = function(e) qr.solve(X2_io, d$R2)
    )
    g2_hat    <- as.vector(S2 %*% psi2_hat)
    eps2_true <- d$R2 - d$m2 - d$A2 * d$g2
    g2_true   <- d$g2

    Y1 <- switch(v,
      "Full"    = d$R2 - d$A2 * g2_hat + abs(g2_hat),
      "Oracle"  = d$m2 + eps2_true + abs(g2_true),
      "No-M"    = eps2_true + abs(g2_true),
      "No-eps2" = d$m2 + abs(g2_true),
      "No-G"    = d$m2 + eps2_true
    )
    S1_tr <- cbind(1, d$O1)
    X1_tr <- d$A1 * S1_tr
    coef1 <- tryCatch(
      solve(crossprod(X1_tr), crossprod(X1_tr, Y1)),
      error = function(e) rep(NA_real_, P_PROB + 1)
    )
    if (anyNA(coef1)) return(NA_real_)
    d_hat_abl <- sign(as.vector(S1_abl_diag %*% coef1))
    d_hat_abl[d_hat_abl == 0] <- 1L
    mean((d_hat_abl != d_true_abl)[near_idx])  # 경계 근처 오류율
  }, numeric(1L))

  data.frame(
    version    = v,
    var_psi    = var(est),
    sd_psi     = sd(est),
    bias       = mean(est) - psi1_true[PSI1_IDX],
    near_err   = mean(near_err_rates, na.rm = TRUE)
  )
}))

cat("\n[ablation 결과 — γ=1, n=500, B=1000]\n")
cat(sprintf("%-12s  %8s  %8s  %8s  %10s\n",
    "버전", "Var(ψ̂₁₁)", "SD(ψ̂₁₁)", "|Bias|", "경계오류율"))
cat(strrep("-", 56), "\n")
for (i in seq_len(nrow(abl_summary))) {
  r <- abl_summary[i, ]
  cat(sprintf("%-12s  %8.4f  %8.4f  %8.4f  %10.4f\n",
      r$version, r$var_psi, r$sd_psi, abs(r$bias), r$near_err))
}

cat("\n[해석]\n")
cat(sprintf("Full → No-M:    Var %.4f → %.4f (%.0f%% 감소)\n",
    abl_summary$var_psi[abl_summary$version=="Full"],
    abl_summary$var_psi[abl_summary$version=="No-M"],
    100*(1 - abl_summary$var_psi[abl_summary$version=="No-M"] /
             abl_summary$var_psi[abl_summary$version=="Full"])))
cat(sprintf("Full → No-eps2: Var %.4f → %.4f (%.0f%% 감소)\n",
    abl_summary$var_psi[abl_summary$version=="Full"],
    abl_summary$var_psi[abl_summary$version=="No-eps2"],
    100*(1 - abl_summary$var_psi[abl_summary$version=="No-eps2"] /
             abl_summary$var_psi[abl_summary$version=="Full"])))
cat(sprintf("Full → No-G:    Var %.4f → %.4f (%.0f%% 감소)\n",
    abl_summary$var_psi[abl_summary$version=="Full"],
    abl_summary$var_psi[abl_summary$version=="No-G"],
    100*(1 - abl_summary$var_psi[abl_summary$version=="No-G"] /
             abl_summary$var_psi[abl_summary$version=="Full"])))

# ════════════════════════════════════════════════════════════
# Phase 3F — 단계별 전면 진단: Stage 2 + Stage 1 체계적 비교
#
# Full Q와 IO-QL 각각 Stage 2, Stage 1을 병렬로 진단해서
# "어느 단계가, 얼마나, 왜" 문제인지 수치로 확정한다.
#
# Stage 2 진단:
#   MSE(ĝ₂) = E[(ĝ₂ - g₂)²]          — blip 추정 오차
#   Bias(ψ̂₂₁) = E[ψ̂₂₁] - 1.0        — O2₁ 계수 편향 (진값=1.0)
#   Var(ψ̂₂₁)                          — O2₁ 계수 분산
#   MSE(m̂₂) (Full Q 한정)            — 주효과 추정 오차
#
# Stage 1 진단:
#   Bias(ψ̂₁₁) = E[ψ̂₁₁] - ψ₁₁*
#   Var(ψ̂₁₁)
#   Acc (n=500, γ=1)
# ════════════════════════════════════════════════════════════
cat("\n─────────────────────────────────────────\n")
cat("Phase 3F: 단계별 전면 진단 (Stage 2 + Stage 1)\n")
cat("─────────────────────────────────────────\n")

N_REPS_3F <- 1000
N_3F      <- 500

# Full Q Stage 2: S2 = (1, O1₁, O1₂, A1, O2₁, O2₂) → 6열
# psi2_hat의 O2₁ 계수 위치: S2 중 O2₁은 5번째 → PSI2_IDX_FQ = 5
PSI2_IDX_FQ <- 5   # Full Q S2 안에서 O2₁ 위치
PSI2_IDX_IO <- 2   # IO-QL S2 = (1, O2₁, O2₂) 안에서 O2₁ 위치
PSI2_TRUE   <- 1.0 # DGP에서 g2 = 1.0 * O2₁

set.seed(SEED)
diag_3f <- vapply(seq_len(N_REPS_3F), function(b) {
  d <- generate_data_prob(N_3F, gamma = 1.0)

  # ── Full Q Stage 2 ─────────────────────────────────────────
  S2_fq  <- cbind(1, d$O1, d$A1, d$O2)
  X2_fq  <- cbind(S2_fq, d$A2 * S2_fq)
  cf2_fq <- tryCatch(
    solve(crossprod(X2_fq), crossprod(X2_fq, d$R2)),
    error = function(e) { v <- qr.solve(X2_fq, d$R2); v[is.na(v)] <- 0; v })
  cf2_fq[is.na(cf2_fq)] <- 0
  q_fq        <- ncol(S2_fq)
  beta2_fq    <- cf2_fq[1:q_fq]
  psi2_fq     <- cf2_fq[(q_fq+1):(2*q_fq)]
  g2hat_fq    <- as.vector(S2_fq %*% psi2_fq)
  m2hat_fq    <- as.vector(S2_fq %*% beta2_fq)
  max_q2_fq   <- m2hat_fq + abs(g2hat_fq)

  # ── IO-QL Stage 2 ─────────────────────────────────────────
  S2_io  <- cbind(1, d$O2)
  X2_io  <- d$A2 * S2_io
  psi2_io <- tryCatch(
    solve(crossprod(X2_io), crossprod(X2_io, d$R2)),
    error = function(e) qr.solve(X2_io, d$R2))
  g2hat_io    <- as.vector(S2_io %*% psi2_io)
  Y1_regret   <- d$R2 - d$A2 * g2hat_io + abs(g2hat_io)

  # ── Full Q Stage 1 ─────────────────────────────────────────
  S1    <- cbind(1, d$O1)
  X1_fq <- cbind(S1, d$A1 * S1)
  cf1_fq <- tryCatch(
    solve(crossprod(X1_fq), crossprod(X1_fq, max_q2_fq)),
    error = function(e) rep(NA_real_, 2*(P_PROB+1)))
  psi1_fq <- if (anyNA(cf1_fq)) NA_real_ else cf1_fq[P_PROB+2]

  # ── IO-QL Stage 1 ─────────────────────────────────────────
  X1_io <- d$A1 * S1
  cf1_io <- tryCatch(
    solve(crossprod(X1_io), crossprod(X1_io, Y1_regret)),
    error = function(e) rep(NA_real_, P_PROB+1))
  psi1_io <- if (anyNA(cf1_io)) NA_real_ else cf1_io[PSI1_IDX]

  c(
    # Stage 2 blip 추정 (O2₁ 계수)
    psi2_fq_val  = psi2_fq[PSI2_IDX_FQ],
    psi2_io_val  = psi2_io[PSI2_IDX_IO],
    # Stage 2 MSE(ĝ₂)
    mse_g2_fq    = mean((g2hat_fq - d$g2)^2),
    mse_g2_io    = mean((g2hat_io - d$g2)^2),
    # Full Q Stage 2 MSE(m̂₂)  ← 오설정 잔차 핵심
    mse_m2_fq    = mean((m2hat_fq - d$m2)^2),
    # Stage 1 blip 추정
    psi1_fq_val  = psi1_fq,
    psi1_io_val  = psi1_io
  )
}, numeric(7L))
# diag_3f: [지표 × 반복] 행렬

psi2_fq_v <- diag_3f["psi2_fq_val", ][!is.na(diag_3f["psi2_fq_val", ])]
psi2_io_v <- diag_3f["psi2_io_val", ][!is.na(diag_3f["psi2_io_val", ])]
psi1_fq_v <- diag_3f["psi1_fq_val", ][!is.na(diag_3f["psi1_fq_val", ])]
psi1_io_v <- diag_3f["psi1_io_val", ][!is.na(diag_3f["psi1_io_val", ])]

cat("\n[Stage 2 진단 — γ=1, n=500, B=1000]\n")
cat(sprintf("%-20s  %8s  %8s  %10s  %10s\n",
    "항목", "Full Q", "IO-QL", "차이", ""))
cat(strrep("-", 65), "\n")
cat(sprintf("%-20s  %8.4f  %8.4f\n", "E[ψ̂₂₁] (진값=1.0)",
    mean(psi2_fq_v), mean(psi2_io_v)))
cat(sprintf("%-20s  %8.4f  %8.4f\n", "|Bias(ψ̂₂₁)|",
    abs(mean(psi2_fq_v) - PSI2_TRUE), abs(mean(psi2_io_v) - PSI2_TRUE)))
cat(sprintf("%-20s  %8.4f  %8.4f\n", "Var(ψ̂₂₁)",
    var(psi2_fq_v), var(psi2_io_v)))
cat(sprintf("%-20s  %8.4f  %8.4f\n", "MSE(ĝ₂)",
    mean(diag_3f["mse_g2_fq", ]), mean(diag_3f["mse_g2_io", ])))
cat(sprintf("%-20s  %8.4f  %8s\n", "MSE(m̂₂) [Full Q만]",
    mean(diag_3f["mse_m2_fq", ]), "(IO-QL은 추정 안 함)"))

cat("\n[Stage 1 진단 — γ=1, n=500, B=1000]\n")
cat(sprintf("%-20s  %8s  %8s\n", "항목", "Full Q", "IO-QL"))
cat(strrep("-", 42), "\n")
cat(sprintf("%-20s  %8.4f  %8.4f  ← 진값 = %.4f\n",
    "E[ψ̂₁₁]", mean(psi1_fq_v), mean(psi1_io_v), psi1_true[PSI1_IDX]))
cat(sprintf("%-20s  %8.4f  %8.4f\n", "|Bias(ψ̂₁₁)|",
    abs(mean(psi1_fq_v) - psi1_true[PSI1_IDX]),
    abs(mean(psi1_io_v) - psi1_true[PSI1_IDX])))
cat(sprintf("%-20s  %8.4f  %8.4f\n", "Var(ψ̂₁₁)",
    var(psi1_fq_v), var(psi1_io_v)))
cat(sprintf("%-20s  %8.4f  %8.4f\n", "SD(ψ̂₁₁)",
    sd(psi1_fq_v), sd(psi1_io_v)))

cat("\n[진단 해석]\n")
cat("Stage 2: Full Q vs IO-QL ψ̂₂₁ 편향 비교\n")
cat(sprintf("  Full Q  |Bias(ψ̂₂₁)| = %.4f → ĝ₂ 추정 자체는 편향 %s\n",
    abs(mean(psi2_fq_v) - PSI2_TRUE),
    ifelse(abs(mean(psi2_fq_v) - PSI2_TRUE) < 0.05, "작음", "있음")))
cat(sprintf("  IO-QL   |Bias(ψ̂₂₁)| = %.4f → ĝ₂ 추정 자체는 편향 %s\n",
    abs(mean(psi2_io_v) - PSI2_TRUE),
    ifelse(abs(mean(psi2_io_v) - PSI2_TRUE) < 0.05, "작음", "있음")))
cat(sprintf("  Full Q  MSE(m̂₂)     = %.4f → 주효과 오설정이 m̂₂ 오차를 만듦\n",
    mean(diag_3f["mse_m2_fq", ])))
cat("\nStage 1: Full Q 편향 원인이 Stage 2가 아니라 Stage 1 회귀 자체임을 확인\n")
cat(sprintf("  Full Q  Stage 1 |Bias| = %.4f (Stage 2 blip은 편향 없는데 Stage 1에서 편향)\n",
    abs(mean(psi1_fq_v) - psi1_true[PSI1_IDX])))
cat(sprintf("  IO-QL   Stage 1 |Bias| = %.4f, Var = %.4f\n",
    abs(mean(psi1_io_v) - psi1_true[PSI1_IDX]),
    var(psi1_io_v)))


# ─── 결론 출력 ─────────────────────────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat("3단계 결론\n")
cat(strrep("=", 60), "\n")
cat("IO-QL Regret: Bias=0 ✓, 단 Var(η₁ᴵᴼ) > Var(η₁ᶠᵘˡˡ)\n")
cat(sprintf("오류 위치: IO-QL 전체 오류율=%.4f, 경계 근처=%.4f vs 먼곳=%.4f\n",
    mean(error_ind_io),
    mean(error_ind_io[margin_io <  median(margin_io)]),
    mean(error_ind_io[margin_io >= median(margin_io)])))
cat("분산 원천 분해 결과 → 4단계 후보 고안 근거로 사용\n")
cat("→ 오류 패턴·분산 구조에서 무엇이 실제 해결 대상인지 확인 후 후보 도출\n")

p_combined <- (p_dist | p_acc_n) / p_decomp / p_3d +
  plot_annotation(
    title    = "3단계: IO-QL 분산 과잉 진단",
    subtitle = "분포 비교 → Acc 수렴 → 분산 분해 → 오류 위치 분석",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave("/Users/songsusie/dtr_simulation/03_variance.png",
       p_combined, width = 14, height = 18, dpi = 150)
cat("\n03_variance.R 완료 — 03_variance.png 저장됨\n")
