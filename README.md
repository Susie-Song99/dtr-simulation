# DTR Simulation: IO-QL Variance Problem and AdaLASSO Solution

**Dynamic Treatment Regime (DTR) 시뮬레이션 포트폴리오**

2단계 순차 치료 결정 문제에서 Full Q-Learning의 편향 문제를 발견하고,  
Interaction-Only Q-Learning(IO-QL)으로 재설계한 뒤, 남은 분산 문제의 원인을 추적해  
Adaptive LASSO로 최종 해결하는 연구 서사를 재현한다.

---

## 연구 서사 (12단계)

| 단계 | 발견 | 핵심 수치 |
|------|------|-----------|
| 1 | Full Q-Learning: γ 증가 → 정책 정확도 하락 | γ=0: Acc=0.981 → γ=2: Acc=0.751 |
| 2 | 원인: E[A1·O1₁·δ₂] ≠ 0 (Prop 3.1) — 이론 확인 | OLS slope β̂=3.001, R²=1.000 |
| 3 | 오류 위치: 결정경계 근처에 집중 | 경계 근처 오류율 0.149, 먼 곳 0.000 |
| 4 | IO-QL 수도아웃컴 재설계: Ỹ₁ = R2 − A2ĝ₂ + \|ĝ₂\| | E[A2·S2·Ỹ₁]=0 — bias=0 이론 성립 |
| 5 | Tian 추정기 Bias=3.0, Regret 추정기 Bias≈0 | n=500, B=1000: Tian=3.02, Regret=0.03 |
| 6 | 새 문제: IO-QL Var(ψ̂₁₁)=0.158, Acc=96.5% | Full Q Acc=0.851 → IO-QL OLS Acc=0.965 |
| 7 | IO-QL 오류도 경계 집중 — 문제 구조 동일 | 경계 근처 오류율 0.071 |
| 8 | 분산 원천 분해: M̃ 지배적 (60%) | Var(Ỹ₁)=20.28, Var(M̃)=12.20 |
| 9 | Stage 진단: Stage 2 blip 편향≈0, Stage 1 Var=0.158 | \|Bias(ψ̂₂₁)\|=0.003, Var(ψ̂₁₁)=0.158 |
| 10 | Ablation: M̃ 제거 → Bias=3.0 복귀, ψ₁* sparse | ψ₁*=(0, 3.85, 0) — O1₂ 불필요 |
| 11 | 6개 후보 비교: Ridge/LASSO/AdaLASSO/ElasticNet/CF/Aug+CF | n=500, B=500, γ=1.0 |
| 12 | **최종 선택: AdaLASSO** | Acc=0.984, 경계오류율=0.033 (OLS 대비 54% 감소) |

---

## 최종 비교 결과

| 방법 | \|Bias\| | Var | Acc | 경계오류율 | 비고 |
|------|---------|-----|-----|-----------|------|
| OLS | 0.000 | 0.167 | 0.965 | 0.071 | 기준 |
| Ridge | 0.211 | 0.141 | 0.965 | 0.071 | Bias 도입, Acc 개선 없음 |
| LASSO | 0.143 | 0.168 | 0.977 | 0.047 | sparsity 부분 회복 |
| **AdaLASSO** | **0.027** | **0.167** | **0.984** | **0.033** | ✅ 최종 채택 |
| ElasticNet | 0.141 | 0.172 | 0.973 | 0.055 | LASSO보다 열위 |
| CF | 0.026 | 0.169 | 0.965 | 0.071 | 분산 구조 미해결 |
| Aug+CF | 2.993 | 0.069 | 0.865 | 0.245 | M̃ 제거 → Bias 복귀 |

**AdaLASSO 선택 근거:**
- Oracle property → Bias 최소화 (0.027)
- ψ₁*=(0, 3.85, 0) sparse 구조 정확히 복원 → 경계 안정화
- Var 감소가 아닌 sparsity 복원이 실제 개선 메커니즘
- Acc 0.984 (전체 최고), 경계오류 54% 감소

---

## 파일 구조

```
dtr_simulation/
├── 00_dgp.R              # 공유 상수, DGP 함수, oracle ψ₁* 계산
├── 01_problem.R          # Full Q-Learning 편향 발견 (단계 1~3)
├── 02_candidates.R       # IO-QL 수도아웃컴 재설계 (단계 4~5)
├── 03_variance.R         # 분산 원인 추적 (단계 6~10)
│                         #   3A: IO-QL vs Oracle 분산 비교
│                         #   3B: 경계 오류 패턴
│                         #   3C: 분산 원천 분해 (M̃, ε₂, G̃)
│                         #   3D: 경계 집중 3D 시각화
│                         #   3E: Ablation — M̃ 제거 실험
│                         #   3F: Stage 2 + Stage 1 전면 진단
├── 04_regularize.R       # 후보 비교 + AdaLASSO 채택 (단계 11~12)
├── 05_summary.R          # 12단계 연구 서사 요약 + 최종 플롯
├── 06_complex_dgp.R      # 복잡 환경 확장 (p=6, train/test split)
├── run_all.R             # 전체 실행 스크립트
└── plots/
    ├── 04_regularize.png
    ├── 05_summary.png
    ├── 06_complex_dgp_data.json   # B=300 시뮬레이션 결과
    └── 06_complex_dgp_viz.html    # 4-method 비교 인터랙티브 시각화
```

---

## 실행 방법

```r
# 개별 실행 (순서대로)
Rscript --no-init-file 00_dgp.R
Rscript --no-init-file 01_problem.R
Rscript --no-init-file 02_candidates.R
Rscript --no-init-file 03_variance.R
Rscript --no-init-file 04_regularize.R
Rscript --no-init-file 05_summary.R
```

```r
# 전체 한 번에
Rscript --no-init-file run_all.R
```

**`--no-init-file` 필수**: `.Rprofile`에 별도 스크립트가 포함되어 있어 충돌 방지를 위해 필요.

---

## DGP 설정

- **Stage 2**: O2₁ = O1₁ + 1.5·A1 + ε,  m₂ = γ·O2₁²,  g₂ = O2₁
- **Stage 1**: R1 = 0, 치료 A_t ∈ {-1, +1}
- **Oracle ψ₁*** = (0, 3.85, 0) — S1 = (1, O1₁, O1₂) 중 O1₁만 유효
- **시뮬레이션 설정**: n=500, B=500~1000, SEED=42

---

## 주요 발견

1. **Full Q-Learning 편향**: m₂ 비선형 시 Stage 1 추정에 bias 전파 → 정책 정확도 붕괴
2. **IO-QL 수도아웃컴**: Regret 추정기가 bias=0 보장 (Tian 추정기 대비 우위)
3. **M̃ 지배 구조**: Var(Ỹ₁)의 60%가 m₂ 주효과 잔차 → bias=0 조건과 trade-off
4. **sparsity가 핵심**: Var 감소로는 경계오류 해결 불가 — ψ₁* sparse 구조 복원이 실제 메커니즘
5. **AdaLASSO**: 2단계(Ridge 초기화 → 가중 LASSO)로 Oracle property 달성

---

---

## 복잡 환경 확장 실험 (06_complex_dgp.R)

p=2 단순 DGP를 p=6으로 확장하고 train/test를 분리한 추가 검증.

### DGP 설계

- p=6, 활성 차원: O1₁, O1₃ (나머지 4개는 비활성 노이즈)
- O2₁ = O1₁ + 1.5·A1 + ε, O2₃ = O1₃ − 1.0·A1 + ε (O1 → O2 carry-forward)
- g₂ = O2₁ + 0.8·O2₃, m₂ = γ(O2₁² + O2₃²) + 0.3·O2₂  (γ=2.0)
- **Oracle ψ₁*** ≈ (0, 6.30, 0, −3.74, 0, 0, 0)

### 결과 (N_train=N_test=600, B=300)

| 방법 | Train Acc | Test Acc | ψ̂₁ sparsity |
|------|-----------|----------|-------------|
| Full Q + OLS | 59.4% | 59.0% | 비활성 계수 오염, 방향 역전 |
| IO-QL + OLS | 94.8% | 94.8% | 모든 계수 비영 |
| IO-QL + LASSO | 94.0% | 93.9% | 비활성 계수 부분 복원 |
| IO-QL + AdaLASSO | 93.9% | 93.9% | 비활성 계수 완전 복원 |

### 주요 발견

- **Full Q**: γ=2.0 비선형 주효과 → 수도아웃컴 오염 → ψ̂₁ 방향 역전 → 59% (무작위)
- **IO-QL vs Full Q**: +35.8%p — 주효과 편향 제거 효과 압도적
- **LASSO/AdaLASSO**: 정확도는 OLS와 유사하나 비활성 차원 계수를 0으로 복원 → 경계가 올바른 (O1₁, O1₃) 2D 부분공간에 집중
- **Train-test 갭**: 4개 방법 모두 ≈0 → p=6, N=600에서 과적합 없음; 개선은 bias 제거에 의한 것

---

## 환경

```
R >= 4.0
glmnet, ggplot2, dplyr, tidyr, patchwork, jsonlite
```
