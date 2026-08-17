# ============================================================
# 00_dgp.R — Data Generating Process + Oracle ψ₁* 계산
#
# 역할:
#   - 모든 실험에서 공통으로 쓰는 DGP 정의
#   - oracle ψ₁*를 수치적으로 계산 (PSI1_TRUE = 1.0 이 아님)
#   - 오설정 강도 γ를 파라미터로 받아 유연하게 설계
#
# DGP:
#   Stage 1: S1 ~ N(0,1),  A1 ~ {-1,+1}
#            m1(S1) = S1²  (비선형 주효과)
#            g1(S1) = ψ1_blip · S1
#            R1 = m1(S1) + A1·g1(S1) + ε1
#
#   Stage 2: S2 = S1 + 0.5·A1 + ν   ← A1이 S2에 영향 (전이 메커니즘)
#            m2(S2) = γ·S2²          ← γ=0: 선형(오설정 없음), γ>0: 비선형
#            g2(S2) = ψ2_blip · S2
#            R2 = m2(S2) + A2·g2(S2) + ε2
#
# 오설정:
#   Full Q는 Stage 2를 선형으로 근사 → δ₂ = m2(S2) - S2ᵀβ₂* ≠ 0 (γ≠0일 때)
#   γ = 0이면 m2가 선형 → 오설정 없음 → Full Q도 unbiased
#   γ 클수록 오설정 심해짐 → 편향 증가
# ============================================================

# ─── 공통 상수 ───────────────────────────────────────────────
SEED        <- 42
PSI1_BLIP   <- 1.0      # R1 안의 블립 계수 (DGP 설계값)
PSI2_BLIP   <- 1.0      # R2 안의 블립 계수 (DGP 설계값)
SIGMA_R     <- 1.0      # 보상 잡음 표준편차
SIGMA_TRANS <- 0.5      # 상태 전이 잡음 표준편차
N_REPS      <- 1000     # 반복 횟수
N_GRID      <- c(100, 200, 500, 1000, 2000, 5000)
N_ORACLE    <- 1e6      # oracle ψ₁* 수치 계산용 샘플 수
N_TEST      <- 1e4      # 정책 평가용 테스트 샘플 수
GAMMA_GRID  <- c(0, 0.5, 1.0, 1.5, 2.0)  # 오설정 강도


# ─── DGP ─────────────────────────────────────────────────────
# gamma: m2의 비선형성 강도. 0이면 m2=0 (선형), 클수록 오설정 심함
generate_data <- function(n, gamma = 1.0) {
  # Stage 1
  S1   <- rnorm(n)
  A1   <- sample(c(-1L, 1L), n, replace = TRUE)
  eps1 <- rnorm(n, 0, SIGMA_R)
  m1   <- S1^2                          # 비선형 주효과
  R1   <- m1 + A1 * PSI1_BLIP * S1 + eps1

  # Stage 2
  # S2가 A1에 의존 → |g2(S2)|가 A1과 얽힘 → Tian 방식 편향 발생 원인
  S2   <- S1 + 0.5 * A1 + rnorm(n, 0, SIGMA_TRANS)
  A2   <- sample(c(-1L, 1L), n, replace = TRUE)
  eps2 <- rnorm(n, 0, SIGMA_R)
  m2   <- gamma * S2^2                  # 비선형 주효과. γ=0이면 m2=0
  R2   <- m2 + A2 * PSI2_BLIP * S2 + eps2

  list(S1=S1, A1=A1, eps1=eps1, R1=R1,
       S2=S2, A2=A2, eps2=eps2, R2=R2,
       m1=m1, m2=m2,
       g1_true = PSI1_BLIP * S1,
       g2_true = PSI2_BLIP * S2)
}


# ─── Oracle ψ₁* 수치 계산 ────────────────────────────────────
# Q1*(S1,A1) = m1 + A1·g1 + E[m2(S2)+|g2(S2)| | S1,A1]
# 이 함수의 A1·S1 계수가 진짜 ψ₁*
# S2가 A1에 의존하므로 E[m2+|g2||S1,A1]도 A1에 의존 → ψ1* ≠ PSI1_BLIP
#
# 수치 계산: 대용량 샘플에서 oracle 수도 아웃컴으로 회귀
compute_oracle_psi1 <- function(gamma = 1.0, n = N_ORACLE) {
  set.seed(SEED)
  S1 <- rnorm(n)
  A1 <- sample(c(-1L, 1L), n, replace = TRUE)
  S2 <- S1 + 0.5 * A1 + rnorm(n, 0, SIGMA_TRANS)

  # oracle 수도 아웃컴: 참값 m2, g2 직접 사용 (추정 없음)
  m2_true  <- gamma * S2^2
  g2_true  <- PSI2_BLIP * S2
  Y1_oracle <- (S1^2 + A1 * PSI1_BLIP * S1) + m2_true + abs(g2_true)

  # Q1* 안의 A1·S1 계수 추출
  X1 <- cbind(1, S1, A1 * S1)
  coef(lm.fit(X1, Y1_oracle))[3]
}


# ─── 오설정 잔차 δ₂ 계산 ──────────────────────────────────────
# δ₂(S2) = m2(S2) - S2ᵀβ₂*
# Full Q는 m2를 선형(β₀+β₁S2)으로 근사 → δ₂가 클수록 편향 큼
compute_delta2 <- function(S2, m2_true) {
  # β₂* = OLS projection of m2 onto span(1, S2)
  X  <- cbind(1, S2)
  beta2_star <- coef(lm.fit(X, m2_true))
  m2_hat     <- beta2_star[1] + beta2_star[2] * S2
  delta2     <- m2_true - m2_hat    # 오설정 잔차
  list(delta2 = delta2, beta2_star = beta2_star)
}


# ─── 사전 계산: gamma별 oracle ψ₁* ───────────────────────────
cat("Oracle ψ₁* 수치 계산 중...\n")
oracle_psi1_by_gamma <- sapply(GAMMA_GRID, function(g) {
  psi <- compute_oracle_psi1(gamma = g)
  cat(sprintf("  γ=%.1f → ψ₁* = %.4f\n", g, psi))
  psi
})
names(oracle_psi1_by_gamma) <- as.character(GAMMA_GRID)

cat("\n[확인] γ=0 (오설정 없음)일 때 ψ₁* ≈", oracle_psi1_by_gamma["0"], "\n")
cat("[확인] γ=1 (오설정 있음)일 때 ψ₁* ≈", oracle_psi1_by_gamma["1"], "\n")
cat("       → S2가 A1에 의존하므로 두 값이 다름\n\n")

# 결과 저장
saveRDS(list(
  oracle_psi1  = oracle_psi1_by_gamma,
  gamma_grid   = GAMMA_GRID
), "/Users/songsusie/dtr_simulation/oracle_psi1.rds")

cat("00_dgp.R 완료 — oracle_psi1.rds 저장됨\n")
