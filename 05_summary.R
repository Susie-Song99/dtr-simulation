# ==============================================================
# 05_summary.R — 연구 서사 전체 요약
#
# 12단계 연구 여정:
#   1. Full Q-Learning 편향 현상 발견 (01)
#   2. 편향 원인 이론 확인 — Prop 3.1 (01)
#   3. 오류 위치 분석 — 결정경계 집중 (01)
#   4. IO-QL 수도아웃컴 재설계 (02)
#   5. Tian vs Regret 편향 비교 (02)
#   6. IO-QL 분산 과잉 발견 (03)
#   7. 오류 위치 분석 — IO-QL도 경계 집중 (03)
#   8. 분산 원천 분해 — M̃ 지배적 (03)
#   9. 단계별 전면 진단 — Stage 2는 무죄, Stage 1이 문제 (03)
#  10. Ablation — M̃ 제거 불가, ψ₁* sparse 구조 확인 (03)
#  11. 후보 비교 — A1~A4, B1, B2 (04)
#  12. 최종 선택 — AdaLASSO (04)
# ==============================================================

source("/Users/songsusie/dtr_simulation/00_dgp.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

cat("══════════════════════════════════════════════════════════════\n")
cat("DTR 시뮬레이션 — 연구 서사 전체 요약\n")
cat("══════════════════════════════════════════════════════════════\n\n")

# ════════════════════════════════════════════════════════════
# 1. 연구 흐름 요약 테이블
# ════════════════════════════════════════════════════════════
cat("─────────────────────────────────────────────────────────────\n")
cat("단계별 발견 및 근거\n")
cat("─────────────────────────────────────────────────────────────\n")

narrative <- data.frame(
  step = 1:12,
  finding = c(
    "Full Q-Learning: γ 증가 → 정책 정확도 하락",
    "원인: E[A1·O1₁·δ₂] ≠ 0 (Prop 3.1) — 이론 slope=3.0, 실측=3.001",
    "오류 위치: 결정경계 근처에 집중, 먼 곳 오류=0",
    "IO-QL 수도아웃컴 재설계: Regret Ỹ₁ = R2-A2ĝ₂+|ĝ₂|",
    "Tian |Bias|=3.0, Regret |Bias|≈0 — 편향 문제 해결",
    "새로운 문제: IO-QL Var(ψ̂₁₁)=0.158, n=500 Acc=96.5%",
    "IO-QL 오류도 결정경계 집중 — 경계 근처 오류율 0.071",
    "분산 원천: Var(Ỹ₁)=20.28, M̃=12.20(60%) 지배",
    "Stage 진단: Stage 2 blip 편향≈0, Stage 1 Var=0.158이 문제",
    "Ablation: M̃ 제거 시 Bias=3.0 복귀, ψ₁*=(0,3.85,0) sparse",
    "후보 비교: Ridge/LASSO/AdaLASSO/ElasticNet/CF/Aug+CF",
    "최종 선택: AdaLASSO — Acc=0.984, 경계오류율=0.033"
  ),
  evidence = c(
    "γ=0: Acc=0.981, γ=2: Acc=0.751 (B=1000, n=500)",
    "OLS slope β̂=3.001, R²=1.000 (N=1e5)",
    "경계 근처 오류율 0.149, 먼 곳 0.000",
    "E[A2·S2·Ỹ₁]=0 보장 → bias=0 이론 성립",
    "n=500 B=1000: Tian Bias=3.02, Regret Bias=0.03",
    "B=500: Full Q Acc=0.851, IO-QL Acc=0.965",
    "10분위 경계 근처 0.071 vs 먼 곳 0.000",
    "Var(M̃)=12.20, Var(ε₂)=1.01, Var(G̃)=0.88",
    "Stage 2: |Bias(ψ̂₂₁)|=0.003, Stage 1: |Bias(ψ̂₁₁)|=0.057",
    "No-M: Bias=2.999, Var=0.009 (Tian 복귀)",
    "n=500, B=500, γ=1.0, 6개 후보",
    "|Bias|=0.027, Acc=0.984, 경계오류 54% 감소 (vs OLS)"
  )
)

for (i in 1:nrow(narrative)) {
  cat(sprintf("[%2d] %s\n     → %s\n\n",
      narrative$step[i], narrative$finding[i], narrative$evidence[i]))
}

# ════════════════════════════════════════════════════════════
# 2. 핵심 수치 요약
# ════════════════════════════════════════════════════════════
cat("─────────────────────────────────────────────────────────────\n")
cat("핵심 수치 요약\n")
cat("─────────────────────────────────────────────────────────────\n")

cat("\n[문제: Full Q vs IO-QL @ n=500, γ=1]\n")
cat(sprintf("  Full Q:   |Bias|=3.84, Acc=0.851\n"))
cat(sprintf("  IO-QL:    |Bias|=0.057, Acc=0.965  ← Bias 해결, 분산 문제 등장\n"))

cat("\n[원인: IO-QL 분산 구조]\n")
cat(sprintf("  Var(Ỹ₁) = 20.28\n"))
cat(sprintf("    M̃(m₂ 주효과 잔차) = 12.20  (60%%) — 구조적, Bias=0 위해 필수\n"))
cat(sprintf("    ε₂(보상 잡음)       =  1.01  ( 5%%)\n"))
cat(sprintf("    G̃(블립 잔차)        =  0.88  ( 4%%)\n"))
cat(sprintf("  Oracle Var(ψ̂₁₁) = 0.122  ← 현 구조 이론 하한\n"))

cat("\n[단계 진단 @ n=500, γ=1]\n")
cat(sprintf("  Stage 2: |Bias(ψ̂₂₁)| = 0.003  → 문제 없음\n"))
cat(sprintf("  Stage 1: Var(ψ̂₁₁)   = 0.158  → 여기가 문제\n"))
cat(sprintf("  ψ₁* = (0, 3.85, 0)  → sparse, O1₂ 불필요\n"))

cat("\n[최종 비교 @ n=500, B=500]\n")
cat(sprintf("  %-12s  |Bias|  Var     Acc    경계오류\n", "방법"))
cat(sprintf("  %-12s  %.3f   %.3f   %.3f   %.3f\n", "OLS",       0.000, 0.167, 0.965, 0.071))
cat(sprintf("  %-12s  %.3f   %.3f   %.3f   %.3f\n", "Ridge",     0.211, 0.141, 0.965, 0.071))
cat(sprintf("  %-12s  %.3f   %.3f   %.3f   %.3f\n", "LASSO",     0.143, 0.168, 0.977, 0.047))
cat(sprintf("  %-12s  %.3f   %.3f   %.3f   %.3f  ← 최종 선택\n", "AdaLASSO", 0.027, 0.167, 0.984, 0.033))
cat(sprintf("  %-12s  %.3f   %.3f   %.3f   %.3f\n", "ElasticNet",0.141, 0.172, 0.973, 0.055))
cat(sprintf("  %-12s  %.3f   %.3f   %.3f   %.3f\n", "CF",        0.026, 0.169, 0.965, 0.071))
cat(sprintf("  %-12s  %.3f   %.3f   %.3f   %.3f\n", "Aug+CF",    2.993, 0.069, 0.865, 0.245))

cat("\n[AdaLASSO 선택 근거]\n")
cat("  1. |Bias| = 0.027 — Oracle property로 bias 최소화 (OLS 0.000과 유사)\n")
cat("  2. Acc = 0.984   — 전체 최고 (OLS 대비 +2.0%p)\n")
cat("  3. 경계오류율 = 0.033 — OLS 대비 54% 감소\n")
cat("  4. 원인 기반: ψ₁* sparse 구조에서 O1₂=0 정확히 복원\n")
cat("     → Var 감소가 아닌 sparsity 복원이 경계 안정화의 실제 메커니즘\n")

# ════════════════════════════════════════════════════════════
# 3. 시각화 — 연구 서사 요약 플롯
# ════════════════════════════════════════════════════════════
theme_set(theme_bw(base_size = 11))
BLUE  <- "#2c6fad"; RED <- "#c0392b"; GREEN <- "#27ae60"; GRAY <- "#7f8c8d"

# 3-1. 방법론 비교 (최종)
summary_final <- data.frame(
  method    = factor(c("OLS","Ridge","LASSO","AdaLASSO","ElasticNet","CF","Aug+CF"),
                     levels=c("OLS","Ridge","LASSO","AdaLASSO","ElasticNet","CF","Aug+CF")),
  abs_bias  = c(0.000, 0.211, 0.143, 0.027, 0.141, 0.026, 2.993),
  var_psi   = c(0.167, 0.141, 0.168, 0.167, 0.172, 0.169, 0.069),
  acc       = c(0.965, 0.965, 0.977, 0.984, 0.973, 0.965, 0.865),
  near_err  = c(0.071, 0.071, 0.047, 0.033, 0.055, 0.071, 0.245),
  class     = c("기준","A1","A2","A3","A4","B1","B2")
)
METHOD_COLORS <- c(
  OLS="#7f8c8d", Ridge="#2980b9", LASSO="#1abc9c",
  AdaLASSO="#27ae60", ElasticNet="#3498db",
  CF="#e67e22", `Aug+CF`="#c0392b"
)

p_final_acc <- ggplot(summary_final, aes(x=method, y=acc, fill=method)) +
  geom_col(width=0.6) +
  geom_text(aes(label=sprintf("%.3f", acc)), vjust=-0.4, size=3.2) +
  scale_fill_manual(values=METHOD_COLORS) +
  scale_y_continuous(limits=c(0.75, 1.02)) +
  geom_hline(yintercept=0.965, linetype="dashed", color=GRAY) +
  annotate("text", x=0.7, y=0.968, label="OLS 기준", size=3, color=GRAY) +
  labs(title="최종 비교: 정책 정확도 (n=500, γ=1)",
       x=NULL, y="Acc") +
  theme(legend.position="none",
        axis.text.x=element_text(angle=20, hjust=1),
        plot.title=element_text(face="bold"))

p_final_nerr <- ggplot(summary_final, aes(x=method, y=near_err, fill=method)) +
  geom_col(width=0.6) +
  geom_text(aes(label=sprintf("%.3f", near_err)), vjust=-0.4, size=3.2) +
  scale_fill_manual(values=METHOD_COLORS) +
  geom_hline(yintercept=0.071, linetype="dashed", color=GRAY) +
  annotate("text", x=0.7, y=0.075, label="OLS 기준", size=3, color=GRAY) +
  labs(title="최종 비교: 경계 근처 오류율 (n=500, γ=1)",
       subtitle="낮을수록 우수",
       x=NULL, y="경계오류율") +
  theme(legend.position="none",
        axis.text.x=element_text(angle=20, hjust=1),
        plot.title=element_text(face="bold"))

# 3-2. 연구 여정: 방법론별 Acc vs 단계
journey <- data.frame(
  stage  = c(1, 2, 3),
  label  = c("Full Q\n(편향)", "IO-QL OLS\n(분산)", "IO-QL AdaLASSO\n(최종)"),
  acc    = c(0.851, 0.965, 0.984),
  berr   = c(0.149, 0.071, 0.033)
)
p_journey_acc <- ggplot(journey, aes(x=stage, y=acc)) +
  geom_line(color=BLUE, linewidth=1.2) +
  geom_point(size=4, color=BLUE) +
  geom_text(aes(label=sprintf("%.3f", acc)), vjust=-1.0, size=3.5) +
  scale_x_continuous(breaks=1:3, labels=journey$label) +
  scale_y_continuous(limits=c(0.8, 1.02)) +
  labs(title="연구 여정: Acc 개선",
       x=NULL, y="Acc (n=500, γ=1)") +
  theme(plot.title=element_text(face="bold"))

p_journey_berr <- ggplot(journey, aes(x=stage, y=berr)) +
  geom_line(color=RED, linewidth=1.2) +
  geom_point(size=4, color=RED) +
  geom_text(aes(label=sprintf("%.3f", berr)), vjust=-1.0, size=3.5) +
  scale_x_continuous(breaks=1:3, labels=journey$label) +
  scale_y_continuous(limits=c(0, 0.20)) +
  labs(title="연구 여정: 경계오류율 감소",
       x=NULL, y="경계 근처 오류율") +
  theme(plot.title=element_text(face="bold"))

p_summary <- (p_journey_acc | p_journey_berr) / (p_final_acc | p_final_nerr) +
  plot_annotation(
    title    = "DTR 시뮬레이션 — 연구 서사 요약",
    subtitle = "Full Q 편향 → IO-QL 재설계 → 분산 원인 추적 → AdaLASSO 최종 채택",
    theme    = theme(plot.title=element_text(face="bold", size=14))
  )

ggsave("/Users/songsusie/dtr_simulation/plots/05_summary.png",
       p_summary, width=14, height=10, dpi=150)

cat("\n그래프 저장: plots/05_summary.png\n")
cat("\n══════════════════════════════════════════════════════════════\n")
cat("05_summary.R 완료\n")
cat("══════════════════════════════════════════════════════════════\n")
