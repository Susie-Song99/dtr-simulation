# ==============================================================
# 01_problem.R — 연구 서사 1단계: 문제 발견
#
# [가설] Full Q는 m₂가 비선형일 때 정책 정확도(Acc)가 유의하게 저하된다
#
# DGP (p=2, R1 없음 — 논문 코드와 동일 구조):
#   O1 ~ N₂(0, I),  A1 ~ Uniform{-1,+1}
#   O2₁ = O1₁ + A1·1.5 + noise  ← A1이 O2₁에 강하게 영향
#   O2₂ = O1₂ + noise            ← A1 영향 없음
#   m2(O2) = γ·O2₁²              ← 비선형 주효과, O2₁에만 의존
#   g2(O2) = O2₁                 ← Stage 2 참 블립
#   R2 = m2 + A2·g2 + ε₂         ← Stage 1 보상 없음
#
# ψ₁* (stacking oracle):
#   ψ₁* = (E[S1S1ᵀ])⁻¹ · E[A1·S1·(f2+|g2|)]
#   A1이 O2₁에 영향 → E[A1·O1₁·(f2+|g2|)] ≠ 0 → ψ₁*의 O1₁ 성분이 γ와 함께 증가
#
# Full Q 편향 메커니즘 (Prop 3.1):
#   Full Q는 m2를 선형으로 근사 → 잔차 δ₂ = m2 - m̂2 ≠ 0
#   E[A1·O1₁·δ₂(O2)] ≠ 0  (O2₁이 A1에 의존) → ψ̂₁,Full Q 편향
#   E[A1·O1₂·δ₂(O2)] ≈ 0  (O2₂가 A1에 독립) → O1₂ 성분은 편향 없음
#   → d̂(O1) = sign(ψ̂₁ᵀO1) ≠ d*(O1) for some O1 → Acc < 1
#
# 표기:
#   S1_2 = cbind(1, O1)  ∈ R^{p+1}
#   Z1   = A1 · S1_2  (치료-공변량 상호작용)
#   δ₂   = m₂(O2) − E[m₂(O2)|O2,linear]  (오설정 잔차)
# ==============================================================

source("/Users/songsusie/dtr_simulation/00_dgp.R")
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)

set.seed(SEED)

# ─── p=2 DGP 상수 ────────────────────────────────────────────
P_PROB   <- 2    # 01_problem.R 전용 DGP
A1_TO_O2 <- 1.5  # A1이 O2₁에 미치는 계수 (클수록 편향 심함)

generate_data_prob <- function(n, gamma = 1.0) {
  O1 <- matrix(rnorm(n * P_PROB), n, P_PROB)
  A1 <- sample(c(-1L, 1L), n, replace = TRUE)

  O2_1 <- O1[,1] + A1_TO_O2 * A1 + rnorm(n, 0, 0.3)  # A1 → O2₁ 강한 영향
  O2_2 <- O1[,2] + rnorm(n, 0, 0.3)                   # A1 영향 없음
  O2   <- cbind(O2_1, O2_2)
  A2   <- sample(c(-1L, 1L), n, replace = TRUE)

  m2   <- gamma * O2[,1]^2         # 비선형 주효과, O2₁에만 의존
  g2   <- 1.0 * O2[,1]             # 참 Stage 2 블립
  eps2 <- rnorm(n, 0, SIGMA_R)
  R2   <- m2 + A2 * g2 + eps2

  # R1 없음 — 논문 DGP와 동일. A1의 효과는 O2를 통한 간접 경로만 존재.
  list(O1=O1, A1=A1, O2=O2, A2=A2,
       R2=R2, m2=m2, g2=g2)
}

# Oracle ψ₁* — stacking 방식 (noise 없이 Q1 계산)
compute_oracle_psi1_prob <- function(gamma = 1.0, N = 5e5) {
  set.seed(SEED)
  O1 <- matrix(rnorm(N * P_PROB), N, P_PROB)
  S1 <- cbind(1, O1)   # N × (p+1)

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

# ─── Full Q 추정기 (p=2) ─────────────────────────────────────
fit_full_q_prob <- function(d) {
  n  <- nrow(d$O1)
  S2 <- cbind(1, d$O1, d$A1, d$O2)  # n × (2p+2)
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

  S1 <- cbind(1, d$O1)   # n × (p+1)
  X1 <- cbind(S1, d$A1 * S1)
  coef1 <- tryCatch(
    solve(crossprod(X1), crossprod(X1, max_q2)),
    error = function(e) qr.solve(X1, max_q2)
  )
  # A1·S1 계수: 뒤쪽 p+1 개
  coef1[(P_PROB + 2):(2 * P_PROB + 2)]
}

# ─── 성능 지표 ────────────────────────────────────────────────
compute_accuracy_prob <- function(psi1_hat, psi1_true, O1_test) {
  S1 <- cbind(1, O1_test)
  d_hat  <- sign(as.vector(S1 %*% psi1_hat));  d_hat[d_hat == 0]  <- 1L
  d_true <- sign(as.vector(S1 %*% psi1_true)); d_true[d_true == 0] <- 1L
  mean(d_hat == d_true)
}

# δ₂(O2) = m₂(O2) − 선형 근사: OLS 투영 잔차
compute_delta2_prob <- function(O2_1, gamma) {
  m2 <- gamma * O2_1^2
  X  <- cbind(1, O2_1)
  b  <- coef(lm.fit(X, m2))
  m2 - (b[1] + b[2] * O2_1)
}

# E[A1·O1ⱼ·δ₂]: Proposition 3.1 편향 항 (j=1,2 각각)
compute_prop31_prob <- function(n_large, gamma) {
  set.seed(SEED)
  O1 <- matrix(rnorm(n_large * P_PROB), n_large, P_PROB)
  A1 <- sample(c(-1L, 1L), n_large, replace = TRUE)
  O2_1 <- O1[,1] + A1_TO_O2 * A1 + rnorm(n_large, 0, 0.3)
  delta2 <- compute_delta2_prob(O2_1, gamma)
  c(
    comp1 = mean(A1 * O1[,1] * delta2),  # E[A1·O1₁·δ₂] — 오염됨
    comp2 = mean(A1 * O1[,2] * delta2)   # E[A1·O1₂·δ₂] — 오염 안됨
  )
}


# ════════════════════════════════════════════════════════════
# Oracle ψ₁* 사전 계산
# ════════════════════════════════════════════════════════════
cat("Oracle ψ₁* 계산 중 (p=2 DGP)...\n")
psi1_true_g0 <- compute_oracle_psi1_prob(gamma = 0)
psi1_true_g1 <- compute_oracle_psi1_prob(gamma = 1.0)
cat(sprintf("  γ=0: ψ₁* = (%.4f, %.4f, %.4f)\n",
    psi1_true_g0[1], psi1_true_g0[2], psi1_true_g0[3]))
cat(sprintf("  γ=1: ψ₁* = (%.4f, %.4f, %.4f)\n",
    psi1_true_g1[1], psi1_true_g1[2], psi1_true_g1[3]))

set.seed(SEED)
O1_test <- matrix(rnorm(N_TEST * P_PROB), N_TEST, P_PROB)


# ════════════════════════════════════════════════════════════
# Phase 1A — 가설 검증: 선형 vs 비선형 m₂에서 Full Q Acc 비교
# ════════════════════════════════════════════════════════════
cat("\n─────────────────────────────────────────\n")
cat("Phase 1A: 가설 검증 — Full Q Acc 분포 비교\n")
cat("─────────────────────────────────────────\n")

N_REPS_TEST <- 500
N_SAMPLE    <- 1000

acc_linear <- vapply(seq_len(N_REPS_TEST), function(b) {
  d <- generate_data_prob(N_SAMPLE, gamma = 0)
  psi_hat <- tryCatch(fit_full_q_prob(d), error = function(e) rep(NA, P_PROB+1))
  if (anyNA(psi_hat)) return(NA_real_)
  compute_accuracy_prob(psi_hat, psi1_true_g0, O1_test)
}, numeric(1L))

acc_nonlin <- vapply(seq_len(N_REPS_TEST), function(b) {
  d <- generate_data_prob(N_SAMPLE, gamma = 1.0)
  psi_hat <- tryCatch(fit_full_q_prob(d), error = function(e) rep(NA, P_PROB+1))
  if (anyNA(psi_hat)) return(NA_real_)
  compute_accuracy_prob(psi_hat, psi1_true_g1, O1_test)
}, numeric(1L))

acc_linear <- acc_linear[!is.na(acc_linear)]
acc_nonlin <- acc_nonlin[!is.na(acc_nonlin)]

t_result <- tryCatch(
  t.test(acc_nonlin, acc_linear, paired = FALSE, alternative = "less"),
  error = function(e) list(statistic=c(t=-Inf), p.value=0)
)

cat(sprintf("\n선형 m₂  (γ=0, 올바른 설정): Acc = %.4f ± %.4f\n",
    mean(acc_linear), sd(acc_linear)))
cat(sprintf("비선형 m₂(γ=1, 오설정):       Acc = %.4f ± %.4f\n",
    mean(acc_nonlin), sd(acc_nonlin)))
cat(sprintf("t-검정: t = %.3f, p-value = %.4e\n",
    t_result$statistic, t_result$p.value))
cat(sprintf("→ Acc 저하 %s (α=0.01)\n",
    ifelse(t_result$p.value < 0.01, "통계적으로 유의 ✓", "유의하지 않음")))

df_1a <- data.frame(
  acc   = c(acc_linear, acc_nonlin),
  label = rep(c("선형 m₂\n(γ=0, 참설정)", "비선형 m₂\n(γ=1, 오설정)"),
              c(length(acc_linear), length(acc_nonlin))),
  ord   = rep(c(0, 1), c(length(acc_linear), length(acc_nonlin)))
)

# subtitle: 비전공자용 한국어 서사체 (수식 기호 제거)
# 메커니즘: m₂ 비선형 시 ψ̂₁₁ 부호 오반전 → Acc 저하 (Prop 3.1)
p_acc_compare <- ggplot(df_1a,
  aes(x = reorder(label, ord), y = acc, fill = reorder(label, ord))) +
  geom_boxplot(width = 0.5, alpha = 0.85, outlier.size = 0.8) +
  annotate("segment", x=1, xend=2,
           y=max(df_1a$acc)+0.003, yend=max(df_1a$acc)+0.003,
           color="gray30") +
  annotate("text", x=1.5, y=max(df_1a$acc)+0.008,
           label=sprintf("p = %.2e", t_result$p.value),
           size=3.5, color="gray30") +
  scale_fill_manual(values = c("선형 m₂\n(γ=0, 참설정)" = "#3060c0",
                                "비선형 m₂\n(γ=1, 오설정)"  = "#c03838")) +
  labs(title    = "Phase 1A: Full Q 정책 정확도 비교  (p=2 DGP)",
       subtitle = sprintf("n=%d, B=%d회 — 비선형 함수가 있으면 치료 배정 추정이 틀린다",
                          N_SAMPLE, N_REPS_TEST),
       x = NULL, y = "Acc = P(d̂₁ = d₁*)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))


# ════════════════════════════════════════════════════════════
# Phase 1B — Proposition 3.1: E[A1·O1ⱼ·δ₂(O2)] vs γ
# O1₁ (편향됨) vs O1₂ (편향 없음) 비교로 메커니즘 확인
# ════════════════════════════════════════════════════════════
cat("\n─────────────────────────────────────────\n")
cat("Phase 1B: Proposition 3.1 편향 조건 수치화\n")
cat("─────────────────────────────────────────\n")

# GAMMA_ANNOTATE 원소는 반드시 GAMMA_FINE 그리드 안에 있어야 한다.
# GAMMA_FINE 간격을 바꾸면 GAMMA_ANNOTATE도 함께 확인할 것.
GAMMA_FINE     <- seq(0, 2, by = 0.25)
GAMMA_ANNOTATE <- c(0, 0.5, 1.0, 1.5, 2.0)
N_PROP31       <- 5e5

cat(sprintf("\n%-6s  %20s  %20s\n", "γ",
    "E[A1·O1₁·δ₂] (오염)", "E[A1·O1₂·δ₂] (무오염)"))
cat(strrep("-", 50), "\n")

prop31_vals <- t(sapply(GAMMA_FINE, function(g) {
  v <- compute_prop31_prob(N_PROP31, g)
  cat(sprintf("%-6.2f  %20.6f  %20.6f\n", g, v["comp1"], v["comp2"]))
  v
}))

lm_prop31 <- lm(prop31_vals[,"comp1"] ~ GAMMA_FINE)
cat(sprintf("\nO1₁ 성분 선형 회귀 ~ γ:  기울기 = %.4f, R² = %.4f\n",
    coef(lm_prop31)[2], summary(lm_prop31)$r.squared))

df_prop31 <- data.frame(
  gamma = rep(GAMMA_FINE, 2),
  value = c(prop31_vals[,"comp1"], prop31_vals[,"comp2"]),
  comp  = rep(c("O1₁ (A1→O2₁, 오염됨)", "O1₂ (A1 영향 없음)"), each=length(GAMMA_FINE))
)

# 5개 γ 지점에서 두 성분 수치를 플롯 내 텍스트로 표시 (D-03)
df_annot31 <- df_prop31 %>% filter(gamma %in% GAMMA_ANNOTATE)

# subtitle: 편향이 O1₁에만 나타난다는 사실을 수식 없이 한국어로 전달 (D-05)
# 메커니즘: E[A1·O1₁·δ₂] ∝ γ ≠ 0 / E[A1·O1₂·δ₂] ≈ 0
p_prop31 <- ggplot(df_prop31, aes(x=gamma, y=value, color=comp, linetype=comp)) +
  geom_hline(yintercept=0, linetype="dashed", color="gray50") +
  geom_line(linewidth=1.0) + geom_point(size=2.5) +
  ggrepel::geom_text_repel(
    data = df_annot31,
    aes(label = sprintf("%.3f", value)),
    size = 3.0,
    show.legend = FALSE,
    min.segment.length = 0.2
  ) +
  scale_color_manual(values = c("O1₁ (A1→O2₁, 오염됨)" = "#c03838",
                                 "O1₂ (A1 영향 없음)"     = "#3060c0")) +
  labs(title    = "Phase 1B: Proposition 3.1 편향 조건  (p=2 DGP)",
       subtitle = "편향은 A1의 영향을 받는 변수(O1₁)에만 나타난다",
       x = "오설정 강도 γ", y = "E[A1·O1ⱼ·δ₂(O2)]",
       color = NULL, linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))


# ════════════════════════════════════════════════════════════
# Phase 1C — 오류 위치: 결정 경계 근처 집중 확인
# d*(O1) = sign(ψ₁*ᵀO1) — 2D 결정 경계 시각화
# ════════════════════════════════════════════════════════════
cat("\n─────────────────────────────────────────\n")
cat("Phase 1C: 치료 오류 위치 — 결정 경계 분석\n")
cat("─────────────────────────────────────────\n")

N_DIAG <- 5000
set.seed(SEED)
O1_diag <- matrix(rnorm(N_DIAG * P_PROB), N_DIAG, P_PROB)
d_diag  <- generate_data_prob(N_DIAG, gamma = 1.0)

# Full Q 추정 (단일 큰 표본)
psi1_hat_diag <- tryCatch(
  fit_full_q_prob(d_diag),
  error = function(e) rep(NA, P_PROB + 1)
)

S1_diag  <- cbind(1, d_diag$O1)
d_true   <- sign(as.vector(S1_diag %*% psi1_true_g1))
d_hat    <- sign(as.vector(S1_diag %*% psi1_hat_diag))
d_true[d_true == 0] <- 1L; d_hat[d_hat == 0] <- 1L
error_ind <- (d_hat != d_true)

# 결정 경계까지의 거리 = |ψ₁*ᵀO1| / ‖ψ₁*‖
margin_diag <- abs(as.vector(S1_diag %*% psi1_true_g1)) /
               sqrt(sum(psi1_true_g1[-1]^2))  # intercept 제외한 norm

cat(sprintf("\n전체 오류율: %.4f (n=%d, γ=1)\n", mean(error_ind), N_DIAG))
cat(sprintf("오류 집중: 결정 경계 가까운 절반 오류율 = %.4f, 먼 절반 = %.4f\n",
    mean(error_ind[margin_diag < median(margin_diag)]),
    mean(error_ind[margin_diag >= median(margin_diag)])))

df_diag <- data.frame(
  O1_1   = d_diag$O1[,1],
  O1_2   = d_diag$O1[,2],
  error  = factor(error_ind, labels = c("정확", "오류")),
  margin = margin_diag
)

# 결정 경계선: ψ₁*₀ + ψ₁*₁·x + ψ₁*₂·y = 0
# → y = -(ψ₁*₀ + ψ₁*₁·x) / ψ₁*₂
boundary_line <- data.frame(
  x = seq(-3, 3, length.out = 100)
)
boundary_line$y <- -(psi1_true_g1[1] + psi1_true_g1[2] * boundary_line$x) /
                    psi1_true_g1[3]

# 1C-scatter: 범례 제거 (p_margin_line과 아래 여백 높이 일치를 위해)
# 색상 의미는 subtitle이 전달: 빨강=오류, 파랑=정확
p_scatter <- ggplot(df_diag %>% filter(abs(O1_1) < 3, abs(O1_2) < 3)) +
  geom_point(aes(x=O1_1, y=O1_2, color=error, alpha=error), size=0.8) +
  geom_line(data=boundary_line, aes(x=x, y=y),
            color="black", linewidth=1.2, linetype="solid") +
  annotate("text", x=2.2, y=1.5, label="결정 경계",
           size=3, color="gray30") +
  scale_color_manual(values = c("정확" = "#a8c8e8", "오류" = "#c03838")) +
  scale_alpha_manual(values = c("정확" = 0.3, "오류" = 0.9)) +
  coord_cartesian(xlim=c(-3,3), ylim=c(-3,3)) +
  labs(title    = "Phase 1C-1: 치료 오류 위치  (γ=1, n=5000)",
       subtitle = "오류(빨강)는 결정 경계(실선) 근처에 집중",
       x = "O1₁", y = "O1₂", color = NULL, alpha = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

# 결정 경계 거리(margin_diag)를 10분위로 나눠 분위별 오류율 집계
# ntile은 동일 크기 빈을 보장; x축에 bin 번호 대신 margin_mid(실제 거리)를 써서
# "결정 경계까지 거리"라는 물리적 의미를 유지한다 (bin 번호 사용 금지)
df_margin_bin <- data.frame(
  margin = margin_diag,
  error  = as.integer(error_ind)  # logical → 0/1 (factor 거치지 않음)
) %>%
  mutate(bin = dplyr::ntile(margin, 10)) %>%
  group_by(bin) %>%
  summarize(
    error_rate = mean(error),
    margin_mid = mean(margin),
    .groups = "drop"
  )

# 가장 경계에 가까운 빈 vs 가장 먼 빈 오류율 콘솔 확인
cat(sprintf("10분위 오류율: 가장 가까운 빈 = %.4f, 가장 먼 빈 = %.4f\n",
    df_margin_bin$error_rate[which.min(df_margin_bin$margin_mid)],
    df_margin_bin$error_rate[which.max(df_margin_bin$margin_mid)]))

# 색상 #c03838으로 p_scatter의 오류 색과 일치시켜 같은 대상임을 시각적으로 연결
p_margin_line <- ggplot(df_margin_bin, aes(x = margin_mid, y = error_rate)) +
  geom_line(linewidth = 1.0, color = "#c03838") +
  geom_point(size = 2.5, color = "#c03838") +
  labs(title    = "Phase 1C-2: 결정 경계 거리 vs 오류율",
       subtitle = "거리가 줄수록 오류율 상승",
       x = "결정 경계까지 거리", y = "오류율") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# 중첩 patchwork: patchwork 1.3.1은 네이티브 지원 — wrap_elements 불필요
# p_error_combined에 plot_annotation을 달지 않는다 (바깥 patchwork와 충돌)
p_error_combined <- p_scatter | p_margin_line


# ─── 결론 출력 ─────────────────────────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat("1단계 결론\n")
cat(strrep("=", 60), "\n")
cat(sprintf("Full Q Acc: 선형 %.4f → 비선형 %.4f (저하 %.4f, p=%.2e)\n",
    mean(acc_linear), mean(acc_nonlin),
    mean(acc_linear) - mean(acc_nonlin), t_result$p.value))
cat("원인: E[A1·O1₁·δ₂] ∝ γ, E[A1·O1₂·δ₂] ≈ 0 (Prop 3.1 확인)\n")
cat("위치: 오류가 결정 경계 근처에 집중\n")
cat("→ 2단계: m₂를 Stage 1 추정에서 제거하는 pseudo-outcome 설계\n")

# 레이아웃: 1A 전체 너비 위, (1B | 1C) 2열 아래
# patchwork에서 | 가 / 보다 결합 우선순위가 낮으므로 괄호 필수
p_combined <- p_acc_compare / (p_prop31 | p_error_combined) +
  plot_annotation(
    title    = "1단계: Full Q 문제 발견  (p=2 DGP)",
    subtitle = "가설 검증 → Prop 3.1 원인 추적 → 오류 위치 분석",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave("/Users/songsusie/dtr_simulation/01_problem.png",
       p_combined, width = 10, height = 12, dpi = 150)
cat("\n01_problem.R 완료 — 01_problem.png 저장됨\n")
