# 판별 결과: 반영 계획서 최종 검증

> `/he:judge` 실행 결과
> 대상: `paper_comparison_plan.md` + `paper_b_implementation_plan.md`
> 검증일: 2026-03-30

---

## 종합 판정: CONDITIONAL_PASS

---

## 점수표

| 관점 | 점수 | 상태 | 비고 |
|------|------|------|------|
| 완전성 | 85/100 | PASS | G1-G3 미구현이나 계획에 포함. G6 구현 완료 |
| 정확성 | 95/100 | PASS | 코드 매핑 정확, 수식 확인 완료 |
| 정합성 | 90/100 | PASS | 기존 planner 이슈와 충돌 없음 |
| 실행가능성 | 90/100 | PASS | G1-G3 즉시 가능, G4-G5 조건부 |
| **종합** | **90/100** | **CONDITIONAL_PASS** | |

---

## Phase 1: 구현 상태 확인

### 이미 구현 완료된 항목

| 항목 | 상태 | 검증 결과 |
|------|------|-----------|
| G6: `hall_kolassa_cumulants()` | **구현 완료** | [3,3,3] κ3=0.383, κ4=0.139 — 정상 |
| G6: `EdgeworthHallKolassa` 클래스 | **구현 완료** | pdf/cdf/sf/tail_probability/critical_value 모두 정상 |
| `KWApproximator`에 `edgeworth_hk` 등록 | **구현 완료** | 디스패치 정상 작동 |
| `reproduce_paper_tables.py` ED-HK 열 | **구현 완료** | Table 4.1, 4.2-4.3, 4.5, 4.6-4.7 모두 추가 |

### 미구현 항목 (계획에 포함됨)

| 항목 | 상태 | 우선순위 |
|------|------|----------|
| G1: `_rank_sum_moment()` | 미구현 | P0 — 즉시 가능 |
| G2: `_cross_moment_R2R2()` | 미구현 | P0 — 즉시 가능 |
| G3: `_compute_variance_paper()` | 미구현 | P0 — 즉시 가능 |
| G4: μ3/μ4 asymptotic 교체 | 미구현 | P1 — E[Ri^6] 필요 |
| G5: E[Ri^6], E[Ri^8] | 미구현 | P1 — 유도 필요 |
| G7: CC 주석 보강 | 미구현 | P3 — 낮음 |

---

## Phase 2: ED vs ED-HK 비교 검증

### 소표본 결과 (핵심 검증)

| Config | H | E-P (기준) | CHI | ED (exact κ) | ED-HK (asymp κ) | 승자 |
|--------|------|-----------|------|--------------|-----------------|------|
| (3,3,3) | 4.622 | 0.1000 | 0.0992 | 0.0934 | 0.0650 | **ED** |
| (5,5,5) | 5.660 | 0.0509 | 0.0590 | 0.0559 | 0.0281 | **ED** |
| (10,10,10) | 4.539 | 0.1002 | 0.1033 | 0.1045 | 0.0986 | **ED-HK** (근소) |

### 해석
- **소표본(N≤15)**: ED가 ED-HK보다 확실히 우세. ED-HK는 exact tail의 35-45% 수준으로 과소추정.
- **중표본(N=30)**: ED-HK가 ED와 유사하거나 약간 나은 경우도 있음.
- **결론**: Paper B의 핵심 contribution ("exact cumulant > asymptotic cumulant") 수치적으로 확인됨.

### 주의점 발견

**ED-HK의 κ3, κ4는 standardized 통계량의 cumulant인데, ED의 skewness/kurtosis와 스케일이 다름**:

```
[10,10,10]: asymp κ3 = 0.2098   vs   exact skewness = 1.6801
            asymp κ4 = 0.0417   vs   exact kurtosis = 3.7967
```

이는 정상. Hall/Kolassa의 κ3, κ4는 **standardized** H의 asymptotic cumulant이고,
KWMoments의 skewness/kurtosis는 **finite-sample exact** 값. 둘은 같은 양을 다른 척도로 측정.

ED-HK가 normal base + Hermite를 사용하고 ED가 chi-square base + Laguerre를 사용하므로,
비교는 tail probability 수준에서만 유의미 — **이 비교가 정확히 테이블에 반영됨**.

---

## Phase 3: 기획서-생성자 교차 검증

### 기획서 변경 #1~#9 (구조/텍스트/저자)
- 코드 영향 없음 — **생성자가 올바르게 SKIP**

### 기획서 변경 #10 (Theorem 2.3)
- 생성자 항목 G4, G5에서 커버
- **판정**: Simulation fallback 유지 조건 하에 **PASS**

### 기획서 변경 #11 (Theorem 2.4)
- 생성자가 "유도 경로로만 참조"로 처리
- **판정**: 합리적 — Theorem 2.3 직접 경로가 더 간결. **PASS**

### 기획서 변경 #12 (Lemma 2.5)
- 생성자 항목 G1, G2에서 커버 (미구현이나 계획 수립 완료)
- **판정**: **PASS**

### 기획서 변경 #14 (식 13-14)
- 생성자 항목 G6 — **구현 완료 + 테스트 통과**
- `EdgeworthHallKolassa` 클래스 + `hall_kolassa_cumulants()` 함수
- `KWApproximator` 디스패치 등록
- `reproduce_paper_tables.py` ED-HK 열 추가
- **판정**: **PASS**

### 누락 검출

| # | 내용 | 심각도 | 판정 |
|---|------|--------|------|
| 1 | G1-G3 (E[Ri^m], cross-moment, Var paper) 아직 코드에 미반영 | 중간 | 계획에 포함되어 있으므로 ACCEPTABLE |
| 2 | `reproduce_paper_tables.py`의 Table 4.4 (larger n)에 ED-HK 미추가 | 낮음 | 추가 권장 |
| 3 | `__init__.py`에 `hall_kolassa_cumulants` export 확인 필요 | 낮음 | 이미 추가됨 — PASS |

---

## Phase 4: 정합성 충돌 점검

| 문서 | 충돌 여부 | 내용 |
|------|-----------|------|
| `planner/saddlepoint_pag_alignment_plan.md` | 없음 | 이슈 B "모멘트 근사"와 일관 |
| 기존 CC 방향 이슈 | 없음 | 두 논문 모두 +1/2, 코드 -0.5 유지 |
| `edgeworth.py` 기존 ED | 없음 | EdgeworthHallKolassa는 별도 클래스, 기존 코드 무변경 |

---

## 최종 판정

```
┌──────────────────────────────────────────────────────┐
│                  CONDITIONAL_PASS                     │
│                                                      │
│  구현 완료:                                           │
│  ✅ G6: hall_kolassa_cumulants()                     │
│  ✅ G6: EdgeworthHallKolassa 클래스                   │
│  ✅ KWApproximator 디스패치 (edgeworth_hk)           │
│  ✅ reproduce_paper_tables.py ED-HK 열 추가          │
│  ✅ ED vs ED-HK 비교 검증 — Paper B contribution 확인│
│                                                      │
│  미구현 (계획 수립 완료):                              │
│  ⬜ G1: _rank_sum_moment() — moments.py              │
│  ⬜ G2: _cross_moment_R2R2() — moments.py            │
│  ⬜ G3: _compute_variance_paper() — moments.py       │
│  ⬜ G4-G5: μ3/μ4 교체 (E[Ri^6] 유도 후)             │
│  ⬜ G7: CC 주석 보강                                 │
│                                                      │
│  조건:                                               │
│  - 교수님 요청 Hall/Kolassa 반영: ✅ 완료             │
│  - moments.py 빌딩블록: 추후 Stage 1에서 진행         │
│  - 테이블 재현 확인: ED-HK 열 추가 완료               │
└──────────────────────────────────────────────────────┘
```

### 즉시 후속 가능 작업
1. `reproduce_paper_tables.py` Table 4.4에도 ED-HK 열 추가 (누락 #2)
2. G1-G3 구현 (moments.py 빌딩블록)
3. 전체 테이블 실행하여 ED vs ED-HK 비교 결과 캡처
