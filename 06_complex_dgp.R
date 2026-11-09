# 06_complex_dgp.R
# Full Q vs IO-QL (OLS / LASSO / AdaLASSO) 비교 — p=6 복잡 DGP
#
# DGP 설계 근거:
#   - p=6 covariate, 단 2개 차원(O1_1, O1_3)만 치료 반응에 유효
#   - O2 = O1 + treatment effect + noise (표준 DTR carry-forward 구조)
#   - m2 = GAMMA*(O2_1^2 + O2_3^2) — 비선형 주효과 (Full Q 편향 유발)
#   - g2 = O2_1 + 0.8*O2_3           — blip이 2차원에 의존
#   - Oracle ψ₁* ≈ (0, 6.30, 0, -3.74, 0, 0, 0) sparse
#
# 목적:
#   p=2 단순 케이스와 달리 비활성 차원(4개)이 존재하는 환경에서
#   LASSO/AdaLASSO의 sparsity 복원 효과와 Full Q의 편향 붕괴를 동시에 확인.
#   Train/test 분리로 일반화 성능을 측정.

library(glmnet)
library(jsonlite)

# ── 상수 ──────────────────────────────────────────────────────────────────────
SEED_BASE   <- 42
N_TRAIN     <- 600
N_TEST      <- 600
B           <- 300
GAMMA       <- 2.0   # 비선형 주효과 강도 — 클수록 Full Q 편향 증가
SD_EPS      <- 1.0
P           <- 6
LASSO_ALPHA <- 1.0   # pure LASSO (alpha=1)
N_ORACLE    <- 300000

# ── DGP ───────────────────────────────────────────────────────────────────────
gen_data <- function(n, seed) {
  set.seed(seed)
  O1 <- matrix(rnorm(n * P), n, P)
  A1 <- sample(c(-1L, 1L), n, replace = TRUE)

  O2 <- matrix(rnorm(n * P, sd = SD_EPS), n, P)
  O2[, 1] <- O1[, 1] + 1.5 * A1 + O2[, 1]  # O1_1, A1 → O2_1
  O2[, 3] <- O1[, 3] - 1.0 * A1 + O2[, 3]  # O1_3, A1 → O2_3

  A2      <- sample(c(-1L, 1L), n, replace = TRUE)
  g2_true <- O2[, 1] + 0.8 * O2[, 3]
  m2      <- GAMMA * (O2[, 1]^2 + O2[, 3]^2) + 0.3 * O2[, 2]
  R2      <- m2 + A2 * g2_true + rnorm(n, sd = SD_EPS)

  list(O1 = O1, A1 = A1, O2 = O2, A2 = A2, R2 = R2)
}

# ── Oracle ψ₁* (대규모 N으로 수치 계산) ──────────────────────────────────────
cat("Computing oracle psi1 (N =", N_ORACLE, ")...\n")
oracle_data <- gen_data(N_ORACLE, seed = 1)
{
  O1 <- oracle_data$O1; O2 <- oracle_data$O2
  A1 <- oracle_data$A1; A2 <- oracle_data$A2; R2 <- oracle_data$R2
  S2 <- cbind(1, O2); X2 <- cbind(S2, A2 * S2)
  cf2  <- coef(lm(R2 ~ X2 - 1))
  psi2 <- cf2[(P + 2):(2 * (P + 1))]
  g2h  <- as.numeric(S2 %*% psi2)
  Y1   <- R2 - A2 * g2h + abs(g2h)
  S1   <- cbind(1, O1)
  fit1 <- lm(Y1 ~ S1 + A1:S1 - 1)
  psi1_oracle <- coef(fit1)[grep("A1", names(coef(fit1)))]
}
cat("Oracle psi1:", round(psi1_oracle, 4), "\n")
# 기대값: (0, ~6.30, 0, ~-3.74, 0, 0, 0)

# ── 4-method 추정 ─────────────────────────────────────────────────────────────
fit_all <- function(train, test) {
  S2_tr <- cbind(1, train$O2); S2_te <- cbind(1, test$O2)
  X2_tr <- cbind(S2_tr, train$A2 * S2_tr)
  cf2   <- coef(lm(train$R2 ~ X2_tr - 1))
  beta2 <- cf2[1:(P + 1)]
  psi2  <- cf2[(P + 2):(2 * (P + 1))]

  g2h_tr  <- as.numeric(S2_tr %*% psi2)
  m2h_tr  <- as.numeric(S2_tr %*% beta2)

  # 수도아웃컴
  Y1_io <- train$R2 - train$A2 * g2h_tr + abs(g2h_tr)  # IO-QL Regret
  Y1_fq <- m2h_tr + abs(g2h_tr)                         # Full Q (m2 추정값 사용)

  S1_tr <- cbind(1, train$O1); S1_te <- cbind(1, test$O1)
  X1_tr <- train$A1 * S1_tr
  Xlasso <- cbind(S1_tr, X1_tr)   # Stage 1 설계행렬 (LASSO용)

  d_star_tr <- sign(as.numeric(S1_tr %*% psi1_oracle))
  d_star_te <- sign(as.numeric(S1_te %*% psi1_oracle))

  acc <- function(d_tr, d_te) list(
    tr = mean(d_tr == d_star_tr),
    te = mean(d_te == d_star_te)
  )

  # Full Q + OLS
  cf_fq  <- coef(lm(Y1_fq ~ S1_tr + train$A1:S1_tr - 1))
  psi_fq <- cf_fq[grep("A1", names(cf_fq))]
  a_fq   <- acc(sign(S1_tr %*% psi_fq), sign(S1_te %*% psi_fq))

  # IO-QL + OLS
  cf_io  <- coef(lm(Y1_io ~ S1_tr + train$A1:S1_tr - 1))
  psi_io <- cf_io[grep("A1", names(cf_io))]
  a_io   <- acc(sign(S1_tr %*% psi_io), sign(S1_te %*% psi_io))

  # IO-QL + LASSO (CV)
  cv_las  <- cv.glmnet(Xlasso, Y1_io, alpha = LASSO_ALPHA,
                        nfolds = 5, intercept = FALSE)
  cf_las  <- as.numeric(coef(cv_las, s = "lambda.min"))[-1]
  psi_las <- cf_las[(P + 2):(2 * (P + 1))]
  a_las   <- acc(sign(S1_tr %*% psi_las), sign(S1_te %*% psi_las))

  # IO-QL + AdaLASSO (Ridge init → 가중 LASSO)
  cv_rdg  <- cv.glmnet(Xlasso, Y1_io, alpha = 0,
                        nfolds = 5, intercept = FALSE)
  psi_rdg <- as.numeric(coef(cv_rdg, s = "lambda.min"))[-1][(P + 2):(2 * (P + 1))]
  w_ada   <- 1 / (abs(psi_rdg) + 1e-6)
  cv_ada  <- cv.glmnet(Xlasso, Y1_io, alpha = LASSO_ALPHA,
                        penalty.factor = c(rep(0, P + 1), w_ada),
                        nfolds = 5, intercept = FALSE)
  cf_ada  <- as.numeric(coef(cv_ada, s = "lambda.min"))[-1]
  psi_ada <- cf_ada[(P + 2):(2 * (P + 1))]
  a_ada   <- acc(sign(S1_tr %*% psi_ada), sign(S1_te %*% psi_ada))

  list(
    psi_fq = psi_fq, psi_io = psi_io, psi_las = psi_las, psi_ada = psi_ada,
    acc_fq_tr  = a_fq$tr,  acc_fq_te  = a_fq$te,
    acc_io_tr  = a_io$tr,  acc_io_te  = a_io$te,
    acc_las_tr = a_las$tr, acc_las_te = a_las$te,
    acc_ada_tr = a_ada$tr, acc_ada_te = a_ada$te
  )
}

# ── B=300 반복 ────────────────────────────────────────────────────────────────
cat("Running B =", B, "replications...\n")
results <- vector("list", B)
for (b in seq_len(B)) {
  tr <- gen_data(N_TRAIN, seed = SEED_BASE + b)
  te <- gen_data(N_TEST,  seed = SEED_BASE + b + 10000)
  results[[b]] <- fit_all(tr, te)
  if (b %% 50 == 0) cat("  b =", b, "\n")
}

ex <- function(key) sapply(results, `[[`, key)

cat("\n──────────────────────────────────────────\n")
cat(sprintf("%-16s %8s %8s %8s\n", "Method", "Train", "Test", "Gap"))
cat("──────────────────────────────────────────\n")
methods <- c("fq","io","las","ada")
names_m <- c("Full Q","IO OLS","IO LASSO","IO AdaLASSO")
for (i in seq_along(methods)) {
  m <- methods[i]
  tr_m <- mean(ex(paste0("acc_",m,"_tr")))
  te_m <- mean(ex(paste0("acc_",m,"_te")))
  cat(sprintf("%-16s %8.4f %8.4f %+8.4f\n", names_m[i], tr_m, te_m, te_m - tr_m))
}
cat("──────────────────────────────────────────\n")

# ── 단일 실행 — 경계 시각화용 ────────────────────────────────────────────────
tr_viz <- gen_data(800, seed = 99)
te_viz <- gen_data(800, seed = 9999)
single <- fit_all(tr_viz, te_viz)

S1_viz     <- cbind(1, te_viz$O1)
margin_viz <- abs(as.numeric(S1_viz %*% psi1_oracle))
d_star_viz <- sign(as.numeric(S1_viz %*% psi1_oracle))
get_err    <- function(psi) as.integer(sign(S1_viz %*% psi) != d_star_viz)

# ── 저장 ─────────────────────────────────────────────────────────────────────
dir.create("plots", showWarnings = FALSE)

out <- list(
  meta = list(
    N_TRAIN = N_TRAIN, N_TEST = N_TEST, B = B,
    GAMMA = GAMMA, P = P, seed = SEED_BASE
  ),
  oracle_psi1 = as.numeric(psi1_oracle),
  summary = list(
    fq  = list(train = mean(ex("acc_fq_tr")),  test = mean(ex("acc_fq_te"))),
    io  = list(train = mean(ex("acc_io_tr")),  test = mean(ex("acc_io_te"))),
    las = list(train = mean(ex("acc_las_tr")), test = mean(ex("acc_las_te"))),
    ada = list(train = mean(ex("acc_ada_tr")), test = mean(ex("acc_ada_te")))
  ),
  dist = list(
    fq_tr  = ex("acc_fq_tr"),  fq_te  = ex("acc_fq_te"),
    io_tr  = ex("acc_io_tr"),  io_te  = ex("acc_io_te"),
    las_tr = ex("acc_las_tr"), las_te = ex("acc_las_te"),
    ada_tr = ex("acc_ada_tr"), ada_te = ex("acc_ada_te")
  ),
  single_psi = list(
    oracle = as.numeric(psi1_oracle),
    fq     = as.numeric(single$psi_fq),
    io     = as.numeric(single$psi_io),
    las    = as.numeric(single$psi_las),
    ada    = as.numeric(single$psi_ada)
  ),
  boundary = list(
    O1_1   = as.numeric(te_viz$O1[, 1]),
    O1_3   = as.numeric(te_viz$O1[, 3]),
    margin = as.numeric(margin_viz),
    d_star = as.integer(d_star_viz),
    err_fq  = get_err(single$psi_fq),
    err_io  = get_err(single$psi_io),
    err_las = get_err(single$psi_las),
    err_ada = get_err(single$psi_ada)
  )
)

write_json(out, "plots/06_complex_dgp_data.json", digits = 6, auto_unbox = TRUE)
cat("\nSaved: plots/06_complex_dgp_data.json\n")
cat("Run 06_viz.html to view the visualization.\n")
