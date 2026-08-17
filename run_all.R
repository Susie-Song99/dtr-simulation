# ============================================================
# run_all.R — 전체 시뮬레이션 순서 실행
#
# 실행 순서:
#   00_dgp.R        → oracle_psi1.rds 생성
#   01_problem.R    → Phase 1~3: 문제 확인 + 원인 추적
#   02_candidates.R → Phase 4: 수도 아웃컴 후보 비교
#   03_variance.R   → Phase 5~6: 분산 문제 + 원천 분해
#   04_regularize.R → Phase 7: 분산 감소 방법
#   05_summary.R    → Phase 8: 최종 비교 + 결론
#
# 사용법:
#   Rscript run_all.R
#   또는 RStudio에서 source("run_all.R")
# ============================================================

BASE_DIR <- "/Users/songsusie/dtr_simulation"

scripts <- c(
  "00_dgp.R",
  "01_problem.R",
  "02_candidates.R",
  "03_variance.R",
  "04_regularize.R",
  "05_summary.R"
)

total   <- length(scripts)
timings <- numeric(total)

cat(strrep("=", 60), "\n")
cat("DTR 시뮬레이션 전체 실행\n")
cat(strrep("=", 60), "\n\n")

for (i in seq_along(scripts)) {
  script_path <- file.path(BASE_DIR, scripts[i])
  cat(sprintf("[%d/%d] %s\n", i, total, scripts[i]))
  cat(strrep("-", 40), "\n")

  t_start <- proc.time()["elapsed"]
  tryCatch(
    source(script_path),
    error = function(e) {
      cat(sprintf("  ✗ 오류 발생: %s\n", conditionMessage(e)))
    }
  )
  t_end <- proc.time()["elapsed"]

  timings[i] <- t_end - t_start
  cat(sprintf("\n  완료: %.1f초\n\n", timings[i]))
}

cat(strrep("=", 60), "\n")
cat("전체 완료\n")
cat(strrep("=", 60), "\n")
cat(sprintf("총 소요 시간: %.1f분 (%.0f초)\n",
    sum(timings)/60, sum(timings)))
cat("\n생성된 파일:\n")
for (i in seq_along(scripts)) {
  img <- gsub("\\.R$", ".png", scripts[i])
  img_path <- file.path(BASE_DIR, img)
  if (file.exists(img_path)) {
    cat(sprintf("  ✓ %s\n", img))
  } else if (scripts[i] == "00_dgp.R") {
    rds_path <- file.path(BASE_DIR, "oracle_psi1.rds")
    cat(sprintf("  %s oracle_psi1.rds\n",
        ifelse(file.exists(rds_path), "✓", "✗")))
  }
}
