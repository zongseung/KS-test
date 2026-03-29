# Peter Hall의 Edgeworth Expansion + Bootstrap 연구 조사

> Paper B에서 참조하는 Hall(1992, 1993), Kolassa(1995) 연구의 핵심 내용과 코드 반영 분석
> 조사일: 2026-03-30

---

## 1. Hall (1992) "The Bootstrap and Edgeworth Expansion" — 핵심 이론

### 1.1 핵심 아이디어: Asymptotic Refinement

| 방법 | CDF 오차 | 비고 |
|------|---------|------|
| Chi-square/Normal 근사 | O(n^{-1/2}) | 1차 점근 |
| Edgeworth 전개 | O(n^{-1}) | 2차 점근 (skewness/kurtosis 보정) |
| **Bootstrap (pivotal statistic)** | **O(n^{-1})** | Edgeworth와 동일 차수 달성 |
| Bootstrap (양측 검정) | O(n^{-2}) | 추가 1차수 개선 |

**핵심**: Bootstrap은 Edgeworth 보정 항을 **암묵적으로 추정**한다. κ3, κ4를 직접 계산하지 않아도,
bootstrap 분포가 자동으로 이 보정을 반영한다.

### 1.2 Pivotal Statistic 조건

Bootstrap이 asymptotic refinement를 달성하려면 통계량이 **asymptotically pivotal**이어야 함:
- 귀무분포가 미지 모수에 의존하지 않아야 함
- **KW H 통계량은 이 조건 충족** — H0 하에서 표본 크기만으로 분포 결정

### 1.3 Cumulant Ordering

Hall이 확립한 핵심 정리: 점근적으로 정규인 통계량의 p번째 cumulant는 O(n^{-(p-2)/2}).
- κ1 = O(1), κ2 = O(1) (mean, variance — 스케일링 후)
- κ3 = O(n^{-1/2}) (skewness)
- κ4 = O(n^{-1}) (kurtosis)
- 고차 cumulant는 더 빠르게 감소 → Edgeworth 급수가 유효한 점근 급수

---

## 2. Hall (1993a,b) — Studentized Statistics + Weak Assumptions

### Hall (1993) Annals of Statistics
- Studentized 통계량에 대해 **최소한의 모멘트 조건**으로 Edgeworth 전개 유효성 증명
- Rank statistics는 bounded → 모든 모멘트 존재 → **조건 자동 충족**
- KW 통계량에 Edgeworth 보정 적용의 이론적 정당화

### Hall (1993) JRSS-B
- 비모수 곡선 추정에서 bootstrap 신뢰대역
- 직접적 관련성은 낮으나, bootstrap + Edgeworth 결합의 방법론적 프레임워크 제공

---

## 3. Kolassa (1995) — Rank-Sum 통계량의 Edgeworth 근사

### 3.1 핵심 발견

> **Rank statistics의 표준 Edgeworth 급수는 격자 보정(lattice correction) 없이도 O(1/n) 정확도 달성**

- 일반적으로 이산 격자 분포에는 Sheppard 보정 + 연속성 보정이 필요
- 그러나 rank statistics는 격자 간격이 n 증가에 따라 자동으로 축소
- 따라서 **표준 Edgeworth 공식을 그대로 적용 가능**

### 3.2 KW 구현에의 시사점
- 현재 `edgeworth.py`가 표준 공식을 사용하는 것은 **이론적으로 정당**
- 격자 보정을 추가할 필요 없음 (중간~대표본)
- Saddlepoint의 연속성 보정(SDC1/SDC2)은 소표본에서만 효과적 — Kolassa의 결과와 일관

---

## 4. Paper B에서의 Asymptotic Cumulant 공식 (식 13-14)

### 4.1 공식의 의미

```
κ3 = (Σ 1/ni² - 1/N²) / (Σ 1/ni - 1/N)^{3/2}    (식 13)
κ4 = (Σ 1/ni³ - 1/N³) / (Σ 1/ni - 1/N)²           (식 14)
```

이 공식은 **standardized H 통계량**의 asymptotic skewness/kurtosis:
- 그룹 간 표본 크기 불균형이 분포 형태에 미치는 영향을 포착
- 균형 설계 (모든 ni 동일) → κ3, κ4가 단순화 → chi-square 근사 최적
- 불균형 설계 → κ3, κ4가 커짐 → chi-square 근사 부정확 → Edgeworth 보정 필수

### 4.2 Hall/Kolassa 기반 vs Exact cumulant 기반

| 방법 | κ3, κ4 원천 | 정확도 | 계산 비용 |
|------|------------|--------|----------|
| Exact (열거) | 정확한 분포에서 직접 계산 | **정확** | N ≤ 15만 가능 |
| Simulation | Monte Carlo 추정 | MC 오차 존재 | O(sim × N) |
| **Hall/Kolassa (식 13-14)** | **Closed-form asymptotic** | O(1/N) 오차 | **O(k)** — 즉시 |
| Chi-square 기반 (현재 코드) | sqrt(8/df), 12/df | O(1) — 매우 부정확 | O(1) |

**결론**: 식 13-14는 simulation과 chi-square 사이의 **최적 지점** — 계산 비용 거의 0, 정확도는 simulation에 근접.

---

## 5. Bootstrap Calibration은 KW에 필요한가?

### 5.1 이론적으로는 유효
- KW H는 asymptotically pivotal → bootstrap refinement 가능
- Permutation distribution = bootstrap distribution under H0
- Hall의 이론에 의해 O(n^{-1}) 정확도 달성 가능

### 5.2 실용적으로는 **이미 구현되어 있음**

현재 코드의 아키텍처가 이미 bootstrap calibration과 동등한 역할을 수행:

| Bootstrap 역할 | 현재 코드의 대응 |
|---------------|----------------|
| Permutation bootstrap under H0 | `ExactDistribution` (소표본) |
| Monte Carlo bootstrap | `_compute_moments_from_simulation()` |
| Asymptotic cumulant 추정 | 식 13-14 추가 예정 (G6) |

### 5.3 결론

> **Bootstrap calibration을 별도 방법으로 추가할 필요 없음**.
> 현재 아키텍처(exact → simulation → asymptotic)가 이미 Hall의 이론적 프레임워크와 일치.
>
> 단, 식 13-14의 asymptotic cumulant을 구현하면 **simulation 없이도** bootstrap 수준의 정확도를
> closed-form으로 달성할 수 있으므로 이것이 가장 실용적인 개선.

---

## 6. Edgeworth vs Saddlepoint 관계

Hall의 프레임워크에서 두 근사의 관계:

| 특성 | Edgeworth | Saddlepoint |
|------|-----------|-------------|
| 오차 유형 | O(n^{-1}) **절대** 오차 | O(n^{-1}) **상대** 오차 |
| 꼬리 정확도 | 중심부 좋음, 꼬리 불안정 | **꼬리에서 더 정확** |
| 이론적 관계 | Saddlepoint의 2차 전개 = Edgeworth | Edgeworth의 exponential tilting = Saddlepoint |
| KW 코드 역할 | 벤치마크 (ED, GC-A) | 주 근사 (SD1, SD2, SDC1, SDC2) |

현재 코드의 설계가 이 이론과 **정확히 일치**: saddlepoint을 주 방법으로, Edgeworth를 벤치마크로 사용.

---

## 7. 코드 반영 시사점

### 변경 필요 없는 항목
- `edgeworth.py`의 chi-square + Laguerre 전개 구조 → **이론적으로 올바름** (Kolassa 1995)
- `saddlepoint.py`의 Lugannani-Rice 공식 → Hall의 Edgeworth-Saddlepoint 관계와 일치
- Exact/Simulation/Asymptotic 3단계 구조 → Bootstrap calibration과 동등

### 반영 권장 항목

| 우선순위 | 항목 | 근거 |
|----------|------|------|
| **P2** | 식 13-14 asymptotic cumulant 구현 (G6) | Hall/Kolassa의 핵심 기여. Simulation 대체 가능 |
| P3 | `edgeworth.py` 안정성 판정 개선 | Hall 이론 기반 — 보정항 크기 모니터링 |
| P4 | 문서에 Hall/Kolassa 이론적 근거 명시 | 현재 구현의 이론적 정당성 기록 |

### 신규 구현 불필요 항목
- Bootstrap calibration 별도 모듈 → **불필요** (현행 simulation이 동등)
- Lattice correction for Edgeworth → **불필요** (Kolassa 1995 결과)
- Bartlett correction → 이론적으로 가능하나 실용적 이점 미미

---

## 8. 요약

| 질문 | 답변 |
|------|------|
| Hall의 bootstrap이 뭔가? | Edgeworth 보정을 암묵적으로 수행하는 resampling 기법 |
| KW에 적용 가능한가? | 예 — H가 asymptotically pivotal이므로 |
| 별도 구현이 필요한가? | **아니오** — 현재 코드의 exact/simulation이 이미 동등한 역할 |
| 그럼 뭘 반영해야 하나? | **식 13-14 (asymptotic κ3, κ4) 구현이 핵심** — 이것이 Hall/Kolassa 연구의 실용적 결실 |
| 현재 코드 구조가 맞나? | **예** — chi-square + Laguerre Edgeworth, Lugannani-Rice saddlepoint 모두 이론적으로 정당 |
