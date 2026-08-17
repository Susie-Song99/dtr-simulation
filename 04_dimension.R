# ==============================================================
# 04_dimension.R — 연구 서사 4단계: 원인 추적 — 고차원 추정 오류
#
# 질문: IO-QL에서 p(공변량 차원)가 커질 때 추정 오류가
#       어디서, 왜 발생하는가?
#
# 메커니즘:
#   Stage 2 OLS: ψ̂₂ ← regress R2 on A2·S2 (S2 dim = p+1)
#   진짜 신호: O2₁만 관련. 나머지 O2₂...O2_p는 잡음 변수.
#   p 증가 → 잡음 변수 p-1개 추가 → Var(ψ̂₂) 증가 → MSE(ĝ₂) 증가
#   → Ỹ₁ 추정 오차 증가 → Var(ψ̂₁) 증가
#
# 분석:
#   4A. Stage 2 추정 오차: MSE(ĝ₂) vs p
#   4B. Var(ψ̂₁₁) vs p: Full vs Oracle vs IO-QL Full
#   4C. Full-Oracle 갭 vs p (추정 오차 기여분 정량화)
#   4D. 경계 오류율 vs p
# ==============================================================

source("/Users/songsusie/dtr_simulation/00_dgp.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

set.seed(SEED)

# ─── 실험 상수 ─────────────────────────────────────────────────
GAMMA_DIM   <- 1.0         # γ 고정 (차원 효과만 분리)
N_FIXED     <- 500         # n 고정
N_REPS_DIM  <- 500         # 몬테카를로 반복
N_TEST_DIM  <- 5000        # 정책 정확도 평가 테스트셋
P_GRID      <- c(2, 5, 10, 20, 50)  # 차원 범위
A1_TO_O2    <- 1.5

# ─── 고차원 DGP ────────────────────────────────────────────────
# p 차원: O1 = (O1₁,...,O1_p), O2₁ = O1₁ + A1_TO_O2·A1 + noise
#          나머지 O2₂...O2_p는 순수 잡음 (신호 없음)
# 진짜 ψ₂*: O2₁ 성분만 0이 아님 → 고차원에서 OLS가 이를 놓침
generate_data_p <- function(n, p, gamma = GAMMA_DIM) {
  O1    <- matrix(rnorm(n * p), n, p)
  A1    <- sample(c(-1L, 1L), n, replace = TRUE)
  O2_1  <- O1[, 1] + A1_TO_O2 * A1 + rnorm(n, 0, 0.3)   # 신호 변수
  O2_noise <- matrix(rnorm(n * (p - 1), 0, 1), n, p - 1) # 잡음 변수
  O2    <- cbind(O2_1, O2_noise)
  A2    <- sample(c(-1L, 1L), n, replace = TRUE)
  m2    <- gamma * O2[, 1]^2    # 비선형 주효과: O2₁만 관련
  g2    <- 1.0 * O2[, 1]        # 블립: O2₁만 관련
  eps2  <- rnorm(n, 0, SIGMA_R)
  R2    <- m2 + A2 * g2 + eps2
  list(O1 = O1, A1 = A1, O2 = O2, A2 = A2,
       R2 = R2, m2 = m2, g2 = g2, eps2 = eps2)
}

# ─── Oracle ψ₁* (p별, MC 계산) ─────────────────────────────────
compute_oracle_psi1_p <- function(p, gamma = GAMMA_DIM, N = 2e5) {
  set.seed(SEED)
  O1 <- matrix(rnorm(N * p), N, p)
  S1 <- cbind(1, O1)
  Q1_given <- function(a1_val) {
    A1v   <- rep(a1_val, N)
    O2_1  <- O1[, 1] + A1_TO_O2 * A1v + rnorm(N, 0, 0.3)
    gamma * O2_1^2 + abs(O2_1)
  }
  Qp <- Q1_given(+1); Qm <- Q1_given(-1)
  S_st <- rbind(S1, S1)
  A_st <- c(rep(1L, N), rep(-1L, N))
  STS  <- crossprod(S_st) / (2 * N)
  STQ  <- crossprod(A_st * S_st, c(Qp, Qm)) / (2 * N)
  as.vector(solve(STS + 1e-8 * diag(p + 1), STQ))
}

# ─── IO-QL Stage 2 추정 (p차원 OLS) ───────────────────────────
fit_ioql_p <- function(d) {
  S2     <- cbind(1, d$O2)      # dim: n × (p+1)
  X2_io  <- d$A2 * S2
  psi2_hat <- tryCatch(
    solve(crossprod(X2_io), crossprod(X2_io, d$R2)),
    error = function(e) qr.solve(X2_io, d$R2)
  )
  g2_hat    <- as.vector(S2 %*% psi2_hat)
  Y1_regret <- d$R2 - d$A2 * g2_hat + abs(g2_hat)
  S1    <- cbind(1, d$O1)
  X1_io <- d$A1 * S1
  coef1 <- tryCatch(
    solve(crossprod(X1_io), crossprod(X1_io, Y1_regret)),
    error = function(e) qr.solve(X1_io, Y1_regret)
  )
  list(psi1 = coef1, g2_hat = g2_hat, g2_true = d$g2)
}

# IO-QL Oracle: Stage 2를 진짜 g2로 대체 (추정 오차 제거)
fit_ioql_oracle_p <- function(d) {
  g2_hat    <- d$g2   # oracle: 추정 없이 진짜 값 사용
  Y1_oracle <- d$R2 - d$A2 * g2_hat + abs(g2_hat)
  S1    <- cbind(1, d$O1)
  X1_io <- d$A1 * S1
  coef1 <- tryCatch(
    solve(crossprod(X1_io), crossprod(X1_io, Y1_oracle)),
    error = function(e) qr.solve(X1_io, Y1_oracle)
  )
  coef1
}

# 테스트셋은 p별로 사전 생성해서 고정 — 루프 안에서 set.seed 호출 금지
# (루프 내 set.seed는 이후 generate_data_p의 RNG 상태를 고정시켜 Var ≈ 0을 만듦)
make_test_set_p <- function(p) {
  set.seed(SEED + 99)
  O1_t <- matrix(rnorm(N_TEST_DIM * p), N_TEST_DIM, p)
  cbind(1, O1_t)
}

compute_acc_p <- function(psi1_hat, psi1_true, S1_test) {
  d_hat  <- sign(as.vector(S1_test %*% psi1_hat));  d_hat[d_hat  == 0] <- 1L
  d_true <- sign(as.vector(S1_test %*% psi1_true)); d_true[d_true == 0] <- 1L
  mean(d_hat == d_true)
}

compute_boundary_err_p <- function(psi1_hat, psi1_true, S1_test, margin_q = 0.20) {
  score_true <- as.vector(S1_test %*% psi1_true)
  margin     <- abs(score_true) / sqrt(sum(psi1_true[-1]^2) + 1e-8)
  d_hat      <- sign(as.vector(S1_test %*% psi1_hat)); d_hat[d_hat == 0] <- 1L
  d_true     <- sign(score_true); d_true[d_true == 0] <- 1L
  near_mask  <- margin <= quantile(margin, margin_q)
  mean(d_hat[near_mask] != d_true[near_mask])
}

# ════════════════════════════════════════════════════════════
# 메인 시뮬레이션: p별 반복 실험
# ════════════════════════════════════════════════════════════
cat("p 차원 분석 시작...\n")
cat(sprintf("  p ∈ {%s}, n=%d, γ=%.1f, B=%d\n",
    paste(P_GRID, collapse = ", "), N_FIXED, GAMMA_DIM, N_REPS_DIM))

results_dim <- lapply(P_GRID, function(p) {
  cat(sprintf("\n  p=%d 처리 중...\n", p))
  psi1_true_p <- compute_oracle_psi1_p(p)
  psi1_idx    <- 2  # O1₁ 성분 (편향·분산 집중 위치)
  S1_test     <- make_test_set_p(p)  # 고정 테스트셋 (루프 밖에서 한 번만)

  psi1_full  <- numeric(N_REPS_DIM)
  psi1_ora   <- numeric(N_REPS_DIM)
  mse_g2     <- numeric(N_REPS_DIM)
  acc_full   <- numeric(N_REPS_DIM)
  acc_ora    <- numeric(N_REPS_DIM)
  berr_full  <- numeric(N_REPS_DIM)
  berr_ora   <- numeric(N_REPS_DIM)

  for (b in seq_len(N_REPS_DIM)) {
    d <- generate_data_p(N_FIXED, p)

    # IO-QL Full
    res <- tryCatch(fit_ioql_p(d), error = function(e) NULL)
    if (!is.null(res)) {
      psi1_full[b] <- res$psi1[psi1_idx]
      mse_g2[b]    <- mean((res$g2_hat - res$g2_true)^2)
      acc_full[b]  <- compute_acc_p(res$psi1, psi1_true_p, S1_test)
      berr_full[b] <- compute_boundary_err_p(res$psi1, psi1_true_p, S1_test)
    } else {
      psi1_full[b] <- NA; mse_g2[b] <- NA
      acc_full[b] <- NA; berr_full[b] <- NA
    }

    # IO-QL Oracle
    psi1_o <- tryCatch(fit_ioql_oracle_p(d), error = function(e) NULL)
    if (!is.null(psi1_o)) {
      psi1_ora[b] <- psi1_o[psi1_idx]
      acc_ora[b]  <- compute_acc_p(psi1_o, psi1_true_p, S1_test)
      berr_ora[b] <- compute_boundary_err_p(psi1_o, psi1_true_p, S1_test)
    } else {
      psi1_ora[b] <- NA; acc_ora[b] <- NA; berr_ora[b] <- NA
    }
  }

  data.frame(
    p              = p,
    psi1_true      = psi1_true_p[psi1_idx],
    # IO-QL Full
    var_full       = var(psi1_full, na.rm = TRUE),
    bias_full      = abs(mean(psi1_full, na.rm = TRUE) - psi1_true_p[psi1_idx]),
    mse_g2_mean    = mean(mse_g2, na.rm = TRUE),
    acc_full_mean  = mean(acc_full, na.rm = TRUE),
    berr_full_mean = mean(berr_full, na.rm = TRUE),
    # IO-QL Oracle
    var_ora        = var(psi1_ora, na.rm = TRUE),
    bias_ora       = abs(mean(psi1_ora, na.rm = TRUE) - psi1_true_p[psi1_idx]),
    acc_ora_mean   = mean(acc_ora, na.rm = TRUE),
    berr_ora_mean  = mean(berr_ora, na.rm = TRUE)
  )
})
results_dim <- do.call(rbind, results_dim)

# ─── Full-Oracle 갭 ──────────────────────────────────────────
# 추정 오차가 Var(ψ̂₁)에 추가하는 기여분
results_dim$gap_var <- results_dim$var_full - results_dim$var_ora

cat("\n\n═══════════════════════════════════════════════════════════\n")
cat("결과 요약: p 차원 vs 추정 오류\n")
cat("═══════════════════════════════════════════════════════════\n")
cat(sprintf("%-6s  %-10s  %-10s  %-10s  %-10s  %-10s\n",
    "p", "MSE(ĝ₂)", "Var(Full)", "Var(Oracle)", "Gap(F-O)", "경계오류율"))
cat(strrep("-", 66), "\n")
for (i in seq_len(nrow(results_dim))) {
  r <- results_dim[i, ]
  cat(sprintf("%-6d  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %-10.4f\n",
      r$p, r$mse_g2_mean, r$var_full, r$var_ora,
      r$gap_var, r$berr_full_mean))
}

# 갭 배율 (p=2 기준)
gap_base <- results_dim$gap_var[results_dim$p == 2]
cat(sprintf("\nFull-Oracle 갭 배율 (p=2 기준):\n"))
for (i in seq_len(nrow(results_dim))) {
  r <- results_dim[i, ]
  cat(sprintf("  p=%-3d: %.1f배\n", r$p, r$gap_var / gap_base))
}

# ════════════════════════════════════════════════════════════
# 시각화
# ════════════════════════════════════════════════════════════
theme_set(theme_bw(base_size = 12))
BLUE  <- "#2c6fad"
RED   <- "#c0392b"
GREEN <- "#27ae60"
GRAY  <- "#7f8c8d"

# 공통 x 스케일
x_scale <- scale_x_continuous(breaks = P_GRID)

# 4A. Stage 2 추정 오차: MSE(ĝ₂) vs p
p_mse <- ggplot(results_dim, aes(x = p, y = mse_g2_mean)) +
  geom_line(color = RED, linewidth = 1.2) +
  geom_point(color = RED, size = 3) +
  x_scale +
  labs(
    title = "4A. Stage 2 추정 오차 vs 공변량 차원",
    subtitle = "MSE(ĝ₂) = E[(ĝ₂ - g₂)²] — 잡음 변수 p-1개 증가로 OLS 정밀도 하락",
    x = "공변량 차원 p", y = "MSE(ĝ₂)"
  )

# 4B. Var(ψ̂₁₁) vs p: Full vs Oracle
var_long <- results_dim |>
  select(p, var_full, var_ora) |>
  pivot_longer(-p, names_to = "method", values_to = "variance") |>
  mutate(method = recode(method,
    var_full = "IO-QL Full (OLS Stage 2)",
    var_ora  = "IO-QL Oracle (Stage 2 oracle)"))

p_var <- ggplot(var_long, aes(x = p, y = variance, color = method, linetype = method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("IO-QL Full (OLS Stage 2)" = BLUE,
                                "IO-QL Oracle (Stage 2 oracle)" = GREEN)) +
  scale_linetype_manual(values = c("IO-QL Full (OLS Stage 2)" = "solid",
                                   "IO-QL Oracle (Stage 2 oracle)" = "dashed")) +
  x_scale +
  labs(
    title = "4B. Var(ψ̂₁₁) vs 공변량 차원",
    subtitle = "Oracle은 p 무관 — Full과의 갭이 추정 오차 기여분",
    x = "공변량 차원 p", y = "Var(ψ̂₁₁)",
    color = NULL, linetype = NULL
  ) +
  theme(legend.position = "bottom")

# 4C. Full-Oracle 갭 vs p (추정 오차 기여분)
p_gap <- ggplot(results_dim, aes(x = p, y = gap_var)) +
  geom_line(color = RED, linewidth = 1.2) +
  geom_point(color = RED, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = GRAY) +
  x_scale +
  labs(
    title = "4C. 추정 오차 기여 분산 (Full-Oracle 갭) vs p",
    subtitle = "p 증가 → OLS 추정 오차가 Var(ψ̂₁₁)에 추가 기여",
    x = "공변량 차원 p", y = "Var(Full) - Var(Oracle)"
  )

# 4D. 경계 오류율 vs p
berr_long <- results_dim |>
  select(p, berr_full_mean, berr_ora_mean) |>
  pivot_longer(-p, names_to = "method", values_to = "berr") |>
  mutate(method = recode(method,
    berr_full_mean = "IO-QL Full",
    berr_ora_mean  = "IO-QL Oracle"))

p_berr <- ggplot(berr_long, aes(x = p, y = berr, color = method, linetype = method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("IO-QL Full" = BLUE, "IO-QL Oracle" = GREEN)) +
  scale_linetype_manual(values = c("IO-QL Full" = "solid", "IO-QL Oracle" = "dashed")) +
  x_scale +
  labs(
    title = "4D. 경계 근처 오류율 vs 공변량 차원",
    subtitle = "결정경계 근처 20% 포인트 기준 — Full은 p에 따라 악화",
    x = "공변량 차원 p", y = "경계 근처 오류율",
    color = NULL, linetype = NULL
  ) +
  theme(legend.position = "bottom")

# 레이아웃 조합
p_combined <- (p_mse | p_var) / (p_gap | p_berr)

ggsave(
  "/Users/songsusie/dtr_simulation/plots/04_dimension.png",
  plot   = p_combined,
  width  = 14,
  height = 10,
  dpi    = 150
)
cat("\n그래프 저장: plots/04_dimension.png\n")
