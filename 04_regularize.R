# ==============================================================
# 04_regularize.R — 연구 서사 4단계: 후보 해결책 비교
#
# 확인된 원인:
#   - IO-QL Stage 1 Var(ψ̂₁₁) = 0.158 (M̃ 구조적 노이즈)
#   - ψ₁* ≈ (0, 3.85, 0) sparse → OLS가 O1₂에 불필요한 분산 할당
#   - Full-Oracle 갭 = 0.033 (Stage 2 추정 오차)
#
# 후보:
#   Class A — Stage 1 정규화 (ψ₁* sparsity 활용)
#     A1: Ridge        — L2, 균등 수축
#     A2: LASSO        — L1, sparse ψ₁*에 적합
#     A3: Adaptive LASSO — Oracle property, bias 최소화
#     A4: Elastic Net  — L1+L2 결합
#
#   Class B — Ỹ₁ 분산 직접 감소 (근본 원인 공략)
#     B1: Cross-fitting   — Stage 2/1 독립 샘플, 갭 0.033 제거
#     B2: Augmented+CF    — m̂₂ 추정 후 M̃ 직접 감소
#
# 비교 기준:
#   |Bias(ψ̂₁₁)|, Var(ψ̂₁₁), Acc(n=500), 경계오류율
# ==============================================================

source("/Users/songsusie/dtr_simulation/00_dgp.R")
library(glmnet)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

set.seed(SEED)

# ─── DGP 상수 (01~03과 동일) ──────────────────────────────────
P_PROB   <- 2
A1_TO_O2 <- 1.5

generate_data_prob <- function(n, gamma = 1.0) {
  O1   <- matrix(rnorm(n * P_PROB), n, P_PROB)
  A1   <- sample(c(-1L, 1L), n, replace = TRUE)
  O2_1 <- O1[, 1] + A1_TO_O2 * A1 + rnorm(n, 0, 0.3)
  O2_2 <- O1[, 2] + rnorm(n, 0, 0.3)
  O2   <- cbind(O2_1, O2_2)
  A2   <- sample(c(-1L, 1L), n, replace = TRUE)
  m2   <- gamma * O2[, 1]^2
  g2   <- 1.0 * O2[, 1]
  eps2 <- rnorm(n, 0, SIGMA_R)
  R2   <- m2 + A2 * g2 + eps2
  list(O1=O1, A1=A1, O2=O2, A2=A2, R2=R2, m2=m2, g2=g2, eps2=eps2)
}

compute_oracle_psi1_prob <- function(gamma = 1.0, N = 5e5) {
  set.seed(SEED)
  O1 <- matrix(rnorm(N * P_PROB), N, P_PROB)
  S1 <- cbind(1, O1)
  Q1_given <- function(a1_val) {
    A1v  <- rep(a1_val, N)
    O2_1 <- O1[, 1] + A1_TO_O2 * A1v + rnorm(N, 0, 0.3)
    gamma * O2_1^2 + abs(O2_1)
  }
  Qp <- Q1_given(+1); Qm <- Q1_given(-1)
  S_st <- rbind(S1, S1); A_st <- c(rep(1L, N), rep(-1L, N))
  STS  <- crossprod(S_st) / (2 * N)
  STQ  <- crossprod(A_st * S_st, c(Qp, Qm)) / (2 * N)
  as.vector(solve(STS + 1e-8 * diag(P_PROB + 1), STQ))
}

cat("Oracle ψ₁* 계산 중...\n")
psi1_true <- compute_oracle_psi1_prob(gamma = 1.0)
cat(sprintf("  ψ₁* = (%.4f, %.4f, %.4f)\n",
    psi1_true[1], psi1_true[2], psi1_true[3]))

PSI1_IDX <- 2  # O1₁ 성분
N_REPS   <- 500
N_TRAIN  <- 500
N_TEST   <- 5000
K_FOLD   <- 5  # Cross-fitting 폴드 수

# ─── 고정 테스트셋 ────────────────────────────────────────────
set.seed(SEED + 99)
O1_test <- matrix(rnorm(N_TEST * P_PROB), N_TEST, P_PROB)
S1_test <- cbind(1, O1_test)
d_true_test  <- sign(as.vector(S1_test %*% psi1_true))
d_true_test[d_true_test == 0] <- 1L
margin_test  <- abs(as.vector(S1_test %*% psi1_true)) /
                sqrt(sum(psi1_true[-1]^2))
near_mask    <- margin_test < median(margin_test)  # 경계 근처 절반

compute_acc <- function(psi1_hat) {
  d_hat <- sign(as.vector(S1_test %*% psi1_hat))
  d_hat[d_hat == 0] <- 1L
  mean(d_hat == d_true_test)
}
compute_near_err <- function(psi1_hat) {
  d_hat <- sign(as.vector(S1_test %*% psi1_hat))
  d_hat[d_hat == 0] <- 1L
  mean((d_hat != d_true_test)[near_mask])
}

# ─── 기준선: IO-QL OLS (03과 동일) ────────────────────────────
fit_ioql_ols <- function(d) {
  S2   <- cbind(1, d$O2); X2_io <- d$A2 * S2
  psi2 <- tryCatch(solve(crossprod(X2_io), crossprod(X2_io, d$R2)),
                   error = function(e) qr.solve(X2_io, d$R2))
  g2h  <- as.vector(S2 %*% psi2)
  Y1   <- d$R2 - d$A2 * g2h + abs(g2h)
  S1   <- cbind(1, d$O1); X1 <- d$A1 * S1
  cf1  <- tryCatch(solve(crossprod(X1), crossprod(X1, Y1)),
                   error = function(e) qr.solve(X1, Y1))
  cf1
}

# ─── Class A: Stage 1 정규화 ──────────────────────────────────
# X1 = A1*(1, O1₁, O1₂), intercept=FALSE (A1 자체가 intercept 역할)
# λ: 5-fold CV로 선택

fit_ioql_ridge <- function(d) {
  S2   <- cbind(1, d$O2); X2_io <- d$A2 * S2
  psi2 <- tryCatch(solve(crossprod(X2_io), crossprod(X2_io, d$R2)),
                   error = function(e) qr.solve(X2_io, d$R2))
  g2h  <- as.vector(S2 %*% psi2)
  Y1   <- d$R2 - d$A2 * g2h + abs(g2h)
  S1   <- cbind(1, d$O1); X1 <- d$A1 * S1
  cv   <- cv.glmnet(X1, Y1, alpha = 0, intercept = FALSE, nfolds = 5)
  cf1  <- as.vector(coef(cv, s = "lambda.min"))[-1]  # intercept 제거
  cf1
}

fit_ioql_lasso <- function(d) {
  S2   <- cbind(1, d$O2); X2_io <- d$A2 * S2
  psi2 <- tryCatch(solve(crossprod(X2_io), crossprod(X2_io, d$R2)),
                   error = function(e) qr.solve(X2_io, d$R2))
  g2h  <- as.vector(S2 %*% psi2)
  Y1   <- d$R2 - d$A2 * g2h + abs(g2h)
  S1   <- cbind(1, d$O1); X1 <- d$A1 * S1
  cv   <- cv.glmnet(X1, Y1, alpha = 1, intercept = FALSE, nfolds = 5)
  cf1  <- as.vector(coef(cv, s = "lambda.min"))[-1]
  cf1
}

fit_ioql_adalasso <- function(d) {
  # Adaptive LASSO: 초기 Ridge 계수로 penalty.factor 계산
  S2   <- cbind(1, d$O2); X2_io <- d$A2 * S2
  psi2 <- tryCatch(solve(crossprod(X2_io), crossprod(X2_io, d$R2)),
                   error = function(e) qr.solve(X2_io, d$R2))
  g2h  <- as.vector(S2 %*% psi2)
  Y1   <- d$R2 - d$A2 * g2h + abs(g2h)
  S1   <- cbind(1, d$O1); X1 <- d$A1 * S1
  # Step 1: Ridge로 초기 추정
  cv_ridge <- cv.glmnet(X1, Y1, alpha = 0, intercept = FALSE, nfolds = 5)
  init_cf  <- as.vector(coef(cv_ridge, s = "lambda.min"))[-1]
  # Step 2: Adaptive weights (진짜 0 계수 보호: 작은 값 → 큰 penalty)
  ada_w <- 1 / (abs(init_cf) + 1e-4)
  # Step 3: 가중 LASSO
  cv_ada <- cv.glmnet(X1, Y1, alpha = 1, intercept = FALSE,
                      penalty.factor = ada_w, nfolds = 5)
  cf1 <- as.vector(coef(cv_ada, s = "lambda.min"))[-1]
  cf1
}

fit_ioql_enet <- function(d) {
  S2   <- cbind(1, d$O2); X2_io <- d$A2 * S2
  psi2 <- tryCatch(solve(crossprod(X2_io), crossprod(X2_io, d$R2)),
                   error = function(e) qr.solve(X2_io, d$R2))
  g2h  <- as.vector(S2 %*% psi2)
  Y1   <- d$R2 - d$A2 * g2h + abs(g2h)
  S1   <- cbind(1, d$O1); X1 <- d$A1 * S1
  cv   <- cv.glmnet(X1, Y1, alpha = 0.5, intercept = FALSE, nfolds = 5)
  cf1  <- as.vector(coef(cv, s = "lambda.min"))[-1]
  cf1
}

# ─── Class B: Ỹ₁ 분산 직접 감소 ──────────────────────────────

# B1: Cross-fitting (K-fold)
# 폴드별로 Stage 2 추정 → 독립 폴드에서 Ỹ₁ 계산 → 전체 ψ̂₁ 추정
fit_ioql_cf <- function(d) {
  n    <- nrow(d$O1)
  fids <- sample(rep(1:K_FOLD, length.out = n))
  Y1_cf <- numeric(n)

  for (k in seq_len(K_FOLD)) {
    train_idx <- which(fids != k)
    test_idx  <- which(fids == k)
    # Stage 2: train fold에서 ĝ₂ 추정
    S2_tr <- cbind(1, d$O2[train_idx, ])
    X2_tr <- d$A2[train_idx] * S2_tr
    psi2  <- tryCatch(
      solve(crossprod(X2_tr), crossprod(X2_tr, d$R2[train_idx])),
      error = function(e) qr.solve(X2_tr, d$R2[train_idx]))
    # test fold에 적용
    S2_te <- cbind(1, d$O2[test_idx, ])
    g2h_te <- as.vector(S2_te %*% psi2)
    Y1_cf[test_idx] <- d$R2[test_idx] - d$A2[test_idx] * g2h_te + abs(g2h_te)
  }
  # Stage 1: 전체 데이터로 OLS
  S1 <- cbind(1, d$O1); X1 <- d$A1 * S1
  cf1 <- tryCatch(solve(crossprod(X1), crossprod(X1, Y1_cf)),
                  error = function(e) qr.solve(X1, Y1_cf))
  cf1
}

# B2: Augmented + Cross-fitting
# 폴드별로 Stage 2 ĝ₂ AND m̂₂ 추정 → Ỹ₁_aug = Ỹ₁ - m̂₂ + Ē[m̂₂]
# m̂₂를 독립 폴드에서 추정해야 E[A1·S1·(m₂-m̂₂)]=0 보장 (Bias=0 유지)
fit_ioql_aug_cf <- function(d) {
  n    <- nrow(d$O1)
  fids <- sample(rep(1:K_FOLD, length.out = n))
  Y1_aug <- numeric(n)
  m2hat_all <- numeric(n)

  for (k in seq_len(K_FOLD)) {
    train_idx <- which(fids != k)
    test_idx  <- which(fids == k)
    # Stage 2: train fold에서 ĝ₂ 추정 (IO-QL)
    S2_tr <- cbind(1, d$O2[train_idx, ])
    X2_tr <- d$A2[train_idx] * S2_tr
    psi2  <- tryCatch(
      solve(crossprod(X2_tr), crossprod(X2_tr, d$R2[train_idx])),
      error = function(e) qr.solve(X2_tr, d$R2[train_idx]))
    S2_te  <- cbind(1, d$O2[test_idx, ])
    g2h_te <- as.vector(S2_te %*% psi2)

    # m₂ 추정: train fold에서 R2 ~ S2 회귀 (주효과만)
    # m2 = γ·O2₁² → 다항 특성 (1, O2₁, O2₁²) 포함해서 비선형 근사
    O2_tr1    <- d$O2[train_idx, 1]
    X_m2_tr   <- cbind(1, O2_tr1, O2_tr1^2)   # 다항 2차
    # A2 중립화: E[R2|O2] = m2, A2 평균=0이므로 OLS로 m₂ 추정 가능
    m2_coef   <- tryCatch(
      solve(crossprod(X_m2_tr), crossprod(X_m2_tr, d$R2[train_idx])),
      error = function(e) qr.solve(X_m2_tr, d$R2[train_idx]))
    O2_te1    <- d$O2[test_idx, 1]
    X_m2_te   <- cbind(1, O2_te1, O2_te1^2)
    m2hat_te  <- as.vector(X_m2_te %*% m2_coef)
    m2hat_all[test_idx] <- m2hat_te

    Y1_aug[test_idx] <- d$R2[test_idx] - d$A2[test_idx] * g2h_te +
                        abs(g2h_te) - m2hat_te
  }
  # 상수 보정: Ē[m̂₂] 더해서 E[Ỹ₁_aug] = E[Ỹ₁] 수준 맞춤
  Y1_aug <- Y1_aug + mean(m2hat_all)

  # Stage 1: OLS
  S1 <- cbind(1, d$O1); X1 <- d$A1 * S1
  cf1 <- tryCatch(solve(crossprod(X1), crossprod(X1, Y1_aug)),
                  error = function(e) qr.solve(X1, Y1_aug))
  cf1
}

# ════════════════════════════════════════════════════════════
# 메인 시뮬레이션
# ════════════════════════════════════════════════════════════
METHODS <- c("OLS", "Ridge", "LASSO", "AdaLASSO", "ElasticNet",
             "CF", "Aug+CF")

cat(sprintf("\n후보 비교 시뮬레이션 (n=%d, B=%d, γ=1.0)\n", N_TRAIN, N_REPS))
cat("방법:", paste(METHODS, collapse=", "), "\n\n")

set.seed(SEED)
results_raw <- lapply(METHODS, function(m) numeric(N_REPS))
names(results_raw) <- METHODS

acc_raw  <- lapply(METHODS, function(m) numeric(N_REPS))
names(acc_raw) <- METHODS
nerr_raw <- lapply(METHODS, function(m) numeric(N_REPS))
names(nerr_raw) <- METHODS

for (b in seq_len(N_REPS)) {
  if (b %% 100 == 0) cat(sprintf("  반복 %d/%d...\n", b, N_REPS))
  d <- generate_data_prob(N_TRAIN, gamma = 1.0)

  fits <- list(
    OLS        = tryCatch(fit_ioql_ols(d),      error = function(e) NULL),
    Ridge      = tryCatch(fit_ioql_ridge(d),    error = function(e) NULL),
    LASSO      = tryCatch(fit_ioql_lasso(d),    error = function(e) NULL),
    AdaLASSO   = tryCatch(fit_ioql_adalasso(d), error = function(e) NULL),
    ElasticNet = tryCatch(fit_ioql_enet(d),     error = function(e) NULL),
    CF         = tryCatch(fit_ioql_cf(d),       error = function(e) NULL),
    `Aug+CF`   = tryCatch(fit_ioql_aug_cf(d),   error = function(e) NULL)
  )

  for (m in METHODS) {
    cf1 <- fits[[m]]
    if (!is.null(cf1) && !anyNA(cf1)) {
      results_raw[[m]][b] <- cf1[PSI1_IDX]
      acc_raw[[m]][b]     <- compute_acc(cf1)
      nerr_raw[[m]][b]    <- compute_near_err(cf1)
    } else {
      results_raw[[m]][b] <- NA
      acc_raw[[m]][b]     <- NA
      nerr_raw[[m]][b]    <- NA
    }
  }
}

# ─── 결과 집계 ─────────────────────────────────────────────────
summary_df <- do.call(rbind, lapply(METHODS, function(m) {
  est  <- results_raw[[m]][!is.na(results_raw[[m]])]
  acc  <- acc_raw[[m]][!is.na(acc_raw[[m]])]
  nerr <- nerr_raw[[m]][!is.na(nerr_raw[[m]])]
  data.frame(
    method   = m,
    bias     = mean(est) - psi1_true[PSI1_IDX],
    abs_bias = abs(mean(est) - psi1_true[PSI1_IDX]),
    var_psi  = var(est),
    sd_psi   = sd(est),
    acc      = mean(acc),
    near_err = mean(nerr),
    n_valid  = length(est)
  )
}))

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("후보 비교 결과 (n=500, γ=1.0, B=500)\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat(sprintf("%-12s  %8s  %8s  %8s  %6s  %8s\n",
    "방법", "|Bias|", "Var(ψ̂)", "SD(ψ̂)", "Acc", "경계오류율"))
cat(strrep("-", 60), "\n")
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf("%-12s  %8.4f  %8.4f  %8.4f  %6.4f  %8.4f\n",
      r$method, r$abs_bias, r$var_psi, r$sd_psi, r$acc, r$near_err))
}

# ─── 기준(OLS) 대비 개선율 ────────────────────────────────────
ols_row <- summary_df[summary_df$method == "OLS", ]
cat("\n[OLS 대비 개선율]\n")
cat(sprintf("%-12s  %10s  %10s  %10s\n", "방법", "Var 감소", "Acc 향상", "경계오류 감소"))
cat(strrep("-", 50), "\n")
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  if (r$method == "OLS") next
  cat(sprintf("%-12s  %9.1f%%  %9.1f%%  %9.1f%%\n",
      r$method,
      100 * (1 - r$var_psi / ols_row$var_psi),
      100 * (r$acc - ols_row$acc) / ols_row$acc,
      100 * (1 - r$near_err / ols_row$near_err)))
}

# ════════════════════════════════════════════════════════════
# 시각화
# ════════════════════════════════════════════════════════════
theme_set(theme_bw(base_size = 11))

# 색상: Class별로 구분
METHOD_COLORS <- c(
  OLS        = "#7f8c8d",
  Ridge      = "#2980b9",
  LASSO      = "#1abc9c",
  AdaLASSO   = "#27ae60",
  ElasticNet = "#3498db",
  CF         = "#e67e22",
  `Aug+CF`   = "#c0392b"
)
CLASS_LABEL <- c(
  OLS="기준", Ridge="A1", LASSO="A2",
  AdaLASSO="A3", ElasticNet="A4", CF="B1", `Aug+CF`="B2"
)

summary_df$class <- CLASS_LABEL[summary_df$method]
summary_df$method_f <- factor(summary_df$method, levels = METHODS)

# 4A. Bias-Variance plot
p_bv <- ggplot(summary_df,
    aes(x = var_psi, y = abs_bias, color = method_f, label = class)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3.5, fontface = "bold", box.padding = 0.4) +
  geom_vline(xintercept = 0.122, linetype="dashed", color="gray50") +
  annotate("text", x=0.122, y=max(summary_df$abs_bias)*0.95,
           label="Oracle 하한\n(0.122)", hjust=-0.05, size=3, color="gray40") +
  scale_color_manual(values = METHOD_COLORS) +
  labs(title    = "4A. Bias-Variance Tradeoff (n=500, γ=1)",
       subtitle = "왼쪽·아래로 갈수록 우수 | 점선 = Oracle 분산 하한",
       x = "Var(ψ̂₁₁)", y = "|Bias(ψ̂₁₁)|",
       color = NULL) +
  theme(legend.position = "none")

# 4B. Acc 비교
p_acc <- ggplot(summary_df,
    aes(x = method_f, y = acc, fill = method_f)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = summary_df$acc[summary_df$method == "OLS"],
             linetype = "dashed", color = "gray50") +
  geom_text(aes(label = sprintf("%.4f", acc)),
            vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = METHOD_COLORS) +
  scale_y_continuous(limits = c(0.85, 1.0)) +
  labs(title = "4B. 정책 정확도 비교 (n=500)",
       subtitle = "점선 = OLS 기준선",
       x = NULL, y = "Acc", fill = NULL) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

# 4C. 경계오류율 비교
p_nerr <- ggplot(summary_df,
    aes(x = method_f, y = near_err, fill = method_f)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = summary_df$near_err[summary_df$method == "OLS"],
             linetype = "dashed", color = "gray50") +
  geom_text(aes(label = sprintf("%.4f", near_err)),
            vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = METHOD_COLORS) +
  labs(title = "4C. 경계 근처 오류율 비교 (n=500)",
       subtitle = "결정경계 근처 50% — 낮을수록 우수",
       x = NULL, y = "경계오류율", fill = NULL) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

# 4D. Var vs Acc scatter
p_vacc <- ggplot(summary_df,
    aes(x = var_psi, y = acc, color = method_f, label = class)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3.5, fontface = "bold", box.padding = 0.4) +
  scale_color_manual(values = METHOD_COLORS) +
  labs(title    = "4D. Var(ψ̂₁₁) vs 정책 정확도",
       subtitle = "오른쪽 위가 이상적 (Var 작고 Acc 높음)",
       x = "Var(ψ̂₁₁)", y = "Acc", color = NULL) +
  theme(legend.position = "none")

p_combined <- (p_bv | p_vacc) / (p_acc | p_nerr) +
  plot_annotation(
    title    = "4단계: 후보 해결책 비교",
    subtitle = "Class A(정규화) vs Class B(분산 직접 감소) vs OLS 기준",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave("/Users/songsusie/dtr_simulation/plots/04_regularize.png",
       p_combined, width = 14, height = 10, dpi = 150)
cat("\n그래프 저장: plots/04_regularize.png\n")
cat("04_regularize.R 완료\n")
