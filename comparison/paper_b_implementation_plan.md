# 반영 계획서 + 판별 보고서: Paper B (final_revision) → 코드

> `/he:generate` + `/he:judge` 통합 문서
> 기획서: `comparison/paper_comparison_plan.md` 기반
> 분석일: 2026-03-30

---

## Phase 1: 코드 현황 분석 (Generate)

### 1.1 moments.py 현황

코드는 3가지 경로로 모멘트를 계산:

| 경로 | 조건 | 방식 | Paper B 공식 사용 여부 |
|------|------|------|---------------------|
| **Exact** | N<=15 (k<=3), N<=13 (k>=4) | ExactDistribution 열거 | N/A (정확) |
| **Simulation** | N>한계, `use_simulation=True` | Monte Carlo 50,000회 | N/A (시뮬) |
| **Asymptotic** | simulation 비활성화 시 | chi-square 기반 스케일링 | **전혀 사용 안 함** |

#### 1.1.1 Variance 계산 (`_compute_variance_exact`, L227-252)
- **현재 코드**: Wallace(1959) 공식 사용
  ```
  Var(H) = 2(k-1) - (2/5)*A/[N(N+1)] - (6/5)*Σ(1/ni)
  where A = 3k(k-2) + N(2k²-6k+1)
  ```
- **Paper B** (Theorem 2.2, p.9): E[Ri^4], E[Ri^2·Rj^2] 기반 공식
  ```
  Var(H) = (12/(N(N+1)))^2 [Σ E[Ri^4]/ni^2 + 2·Σ_{i<j} E[Ri^2·Rj^2]/(ni·nj)] - (k-1)^2
  ```
- **판정**: Wallace 공식은 동치이므로 **수치적으로 동일**. 변경 불필요. 단, Paper B는 E[Ri^4]과 E[Ri^2·Rj^2]를 명시적으로 제공하므로, 향후 μ3/μ4 계산의 빌딩블록으로 활용 가능.

#### 1.1.2 μ3 계산 (3차 raw moment, L302-309)
- **현재 코드** (asymptotic 경로): chi-square skewness 기반
  ```python
  skew = sqrt(8.0 / df)    # chi-square(df)의 skewness
  m3_central = skew * var^(3/2)
  ```
- **Paper B** (Theorem 2.3, 식 8):
  ```
  μ3 = (12^3 / N^3(N+1)^3) · Σ E[Ri^6]/ni^3 - 3·E(H)·Var(H) - μ1^3
  ```
  여기서 E[Ri^6]은 Lemma 2.5의 E[Ri^m] 공식들로부터 유도 가능 (직접 제시되지는 않지만, 조합론 공식으로 계산 가능).
- **불일치**: **심각**. chi-square skewness는 O(N^0) 근사이고, Paper B 공식은 exact.
- **수치 증거** (기존 planner 문서): simulation skewness 1.7147 vs 코드 2.0000 (chi-square)

#### 1.1.3 μ4 계산 (4차 raw moment, L311-318)
- **현재 코드** (asymptotic 경로): chi-square kurtosis 기반
  ```python
  kurt = 12.0 / df    # chi-square(df)의 excess kurtosis
  m4_central = (3 + kurt) * var^2
  ```
- **Paper B** (Theorem 2.3, 식 9):
  ```
  μ4 = (12^4 / N^4(N+1)^4) · Σ E[Ri^8]/ni^4 - 4·μ1·μ3 - 6·μ1^2·Var(H) - μ1^4
  ```
- **불일치**: **심각**. chi-square kurtosis는 O(N^0) 근사.
- **수치 증거**: simulation excess kurt 3.8499 vs 코드 6.0000 (chi-square)

#### 1.1.4 E[Ri^m] 빌딩블록 (Lemma 2.5)
- **현재 코드**: E[Ri^m] 공식이 **전혀 구현되어 있지 않음**
- **Paper B** (Lemma 2.5, p.11):
  ```
  E[Ri]   = ni·(N+1)/2
  E[Ri^2] = ni·(N-ni)·(N+1)/12 + ni^2·(N+1)^2/4
  E[Ri^3] = ni^2·(N+1)^3/8 + ni·(N-ni)·(N+1)·(N+2)/24
  E[Ri^4] = ni·(N-ni)·(N+1)·(2N+1)/30 + ni^2·(N+1)^2·(N-ni)/12 + ni^4·(N+1)^4/16
  ```
- **판정**: 이 공식들은 μ3, μ4를 asymptotic으로 계산하기 위한 핵심 빌딩블록. **신규 구현 필요**.

#### 1.1.5 E[Ri^2·Rj^2] Cross-moment (Lemma 2.5)
- **현재 코드**: **구현 없음**
- **Paper B** (Lemma 2.5, p.11, 식 7):
  ```
  E[Ri^2·Rj^2] = E[Ri^2]·E[Rj^2] - ni·nj·(N+1)^2·(N-2)·(N-ni-nj) / (12·(N-1)·(N-3))
  ```
- **판정**: **신규 구현 필요**.

#### 1.1.6 Theorem 2.4 (A 기반 μ 계산)
- **현재 코드**: A = Σ Ri^2/ni 개념 **없음**
- **Paper B**: A 기반 binomial expansion으로 μ1~μ4를 E[A^m]으로 환원
- **판정**: Theorem 2.3의 직접 공식이 더 간결하므로, Theorem 2.4는 유도 경로로만 참조. **구현은 Theorem 2.3 경로 우선**.

### 1.2 edgeworth.py 현황

#### 1.2.1 Asymptotic cumulant (식 13-14)
- **현재 코드**: Hall/Kolassa asymptotic κ3, κ4 공식 **없음**
  - 모멘트를 KWMoments에서 가져옴 (exact/simulation/chi-square 경로)
- **Paper B** (식 13-14):
  ```
  κ3_asymp = (Σ 1/ni^2 - 1/N^2) / (Σ 1/ni - 1/N)^{3/2}
  κ4_asymp = (Σ 1/ni^3 - 1/N^3) / (Σ 1/ni - 1/N)^2
  ```
- **판정**: 대표본에서 simulation도 비용이 높을 때 사용할 수 있는 closed-form. **추가 구현 권장** (우선순위 낮음).

#### 1.2.2 Edgeworth 구현 방식
- chi-square 기반 Laguerre 전개 사용 — **Paper B Section 3.3과 일치**
- `cdf_normal_based`는 Hermite 기반 normal Edgeworth (GC-A와 유사) — Paper B도 이를 reference benchmark으로 제시

### 1.3 saddlepoint.py 현황

#### 1.3.1 연속성 보정 방향
- **현재 코드** (L378, L422): `x = x - 0.5`, `v = v - 0.5`
- **코드 주석**: "paper text writes +0.5 but table values match -0.5"
- **Paper B** (p.14): `v + 1/2` (Paper A와 동일, 변경 없음)
- **판정**: 코드가 `-0.5`를 쓰는 근거는 테이블 재현. 두 논문 버전 모두 `+1/2` 표기이므로 **이슈 유지**.
  - 해석: Paper의 `v+1/2`는 연속성 보정 적용 시 **이산 격자에서 다음 격자점**을 뜻하므로, P(H ≥ v)를 P_cont(H ≥ v - 1/2)로 치환하는 것과 수학적으로 동치. 코드의 `-0.5`가 올바를 가능성 높음.

#### 1.3.2 CGF 구현
- ER1, Wang, KT 등 **Paper B Section 3.2와 일치** — 변경 불필요

### 1.4 pam.py 현황
- Paper B Section 3.1과 **완전 일치** — 변경 불필요
- 모멘트 정확도에 의존하므로 moments.py 개선 시 자동으로 혜택

---

## Phase 2: 반영 계획 (Generate Output)

### 2.1 코드 수정 필요 항목

| # | 우선순위 | 변경 근거 | 대상 파일:위치 | 현재 코드 | 수정 방향 |
|---|----------|-----------|---------------|-----------|-----------|
| G1 | **P0** | Lemma 2.5 | `moments.py` (신규 함수) | E[Ri^m] 없음 | `_compute_rank_sum_moment(ni, N, m)` 함수 추가 (m=1..4) |
| G2 | **P0** | 식 7 | `moments.py` (신규 함수) | E[Ri^2·Rj^2] 없음 | `_compute_cross_moment_R2R2(ni, nj, N)` 함수 추가 |
| G3 | **P0** | Theorem 2.2 | `moments.py` (신규 함수) | 직접 구현 없음 | `_compute_variance_paper()` — E[Ri^4], cross-moment 기반 Var(H) |
| G4 | **P1** | Theorem 2.3 | `moments.py:L302-318` | chi-square skew/kurt | μ3, μ4를 Paper B 식 8-9로 교체 |
| G5 | **P1** | Theorem 2.3 확장 | `moments.py` (신규) | E[Ri^6], E[Ri^8] 없음 | 고차 모멘트 E[Ri^6], E[Ri^8] 공식 유도 또는 조합론 계산 |
| G6 | **P2** | 식 13-14 | `edgeworth.py` (신규) | asymptotic κ 없음 | `asymptotic_cumulants(sample_sizes)` 함수 추가 |
| G7 | **P3** | 연속성 보정 | `saddlepoint.py:L378,422` | 주석만 존재 | 주석 보강 (수학적 동치 설명) |

### 2.2 신규 구현 상세

#### G1: `_compute_rank_sum_moment(ni, N, m)` — E[Ri^m]

```python
@staticmethod
def _rank_sum_moment(ni: int, N: int, m: int) -> float:
    """
    Compute E[Ri^m] for rank sum Ri of group with size ni.
    Paper B, Lemma 2.5.
    """
    if m == 1:
        return ni * (N + 1) / 2
    elif m == 2:
        return ni * (N - ni) * (N + 1) / 12 + ni**2 * (N + 1)**2 / 4
    elif m == 3:
        return (ni**2 * (N + 1)**3 / 8
                + ni * (N - ni) * (N + 1) * (N + 2) / 24)
    elif m == 4:
        return (ni * (N - ni) * (N + 1) * (2*N + 1) / 30
                + ni**2 * (N + 1)**2 * (N - ni) / 12
                + ni**4 * (N + 1)**4 / 16)
    else:
        raise NotImplementedError(f"E[Ri^{m}] not implemented for m > 4")
```

#### G2: `_cross_moment_R2R2(ni, nj, N)` — E[Ri^2·Rj^2]

```python
@staticmethod
def _cross_moment_R2R2(ni: int, nj: int, N: int) -> float:
    """
    Compute E[Ri^2 * Rj^2] for i != j.
    Paper B, Lemma 2.5, Eq. (7).
    """
    E_Ri2 = ni * (N - ni) * (N + 1) / 12 + ni**2 * (N + 1)**2 / 4
    E_Rj2 = nj * (N - nj) * (N + 1) / 12 + nj**2 * (N + 1)**2 / 4

    if N <= 3:
        return E_Ri2 * E_Rj2

    correction = (ni * nj * (N + 1)**2 * (N - 2) * (N - ni - nj)
                  / (12 * (N - 1) * (N - 3)))
    return E_Ri2 * E_Rj2 - correction
```

#### G3: `_compute_variance_paper()` — Theorem 2.2 기반

```python
def _compute_variance_paper(self) -> float:
    """
    Var(H) via Paper B Theorem 2.2 using E[Ri^4] and E[Ri^2*Rj^2].

    Var(H) = (12/(N(N+1)))^2 * [Σ E[Ri^4]/ni^2 + 2·Σ_{i<j} E[Ri^2*Rj^2]/(ni*nj)]
             - (k-1)^2
    """
    N, k, n = self.N, self.k, self.sample_sizes
    factor = (12.0 / (N * (N + 1)))**2

    sum_diag = 0.0
    for i in range(k):
        E_Ri4 = self._rank_sum_moment(int(n[i]), N, 4)
        sum_diag += E_Ri4 / n[i]**2

    sum_cross = 0.0
    for i in range(k):
        for j in range(i+1, k):
            E_RiRj = self._cross_moment_R2R2(int(n[i]), int(n[j]), N)
            sum_cross += E_RiRj / (n[i] * n[j])

    return factor * (sum_diag + 2 * sum_cross) - (k - 1)**2
```

#### G4: μ3, μ4 교체 — Theorem 2.3 기반

**핵심 문제**: Theorem 2.3의 μ3 공식에 E[Ri^6]이 필요하고, μ4에 E[Ri^8]이 필요함.
Lemma 2.5는 E[Ri^m]을 m=1..4까지만 제공. E[Ri^6], E[Ri^8]은 논문에 명시되지 않음.

**대안 접근법 (Theorem 2.4 + Lemma 2.5 경로)**:
- A = Σ Ri^2/ni 정의
- μ3 = (12/(N(N+1)))^3 · E[A^3] - 3·3(N+1)·(12/(N(N+1)))^2 · E[A^2] + ...
- E[A^3] = Σ E[Ri^6]/ni^3 + 3·Σ_{i≠j} E[Ri^4·Rj^2]/(ni^2·nj) + 6·Σ_{i<j<l} E[Ri^2·Rj^2·Rl^2]/(ni·nj·nl)
- 이 경로도 E[Ri^6], E[Ri^4·Rj^2], E[Ri^2·Rj^2·Rl^2] 등 고차 교차모멘트 필요

**실현 가능한 접근**:
1. **E[Ri^6], E[Ri^8]를 조합론 공식으로 유도** — sampling without replacement의 고차 모멘트 (Kendall & Stuart, 1979)
2. **Theorem 2.4의 직접 계산** — E[A^m]을 multinomial 전개 + 개별/교차 모멘트로 계산
3. **Simulation fallback 유지** — 현재 코드처럼, 정확한 closed-form이 없으면 simulation

**권장**: 접근법 1과 3의 하이브리드
- E[Ri^m] (m≤4)은 Lemma 2.5로 exact 계산
- Var(H)는 Theorem 2.2 공식으로 exact 계산 (Wallace 대체)
- μ3, μ4는 N이 작으면 exact distribution, 그 외 simulation (현행 유지)
- asymptotic 경로는 chi-square 대신 **Lemma 2.5 기반 partial formula** 사용

#### G6: Asymptotic cumulants (식 13-14)

```python
@staticmethod
def asymptotic_cumulants(sample_sizes):
    """
    Compute asymptotic κ3, κ4 per Hall(1992)/Kolassa(1995).
    Paper B, Equations (13)-(14).
    """
    n = np.array(sample_sizes)
    N = np.sum(n)

    sum_inv = np.sum(1.0 / n) - 1.0 / N
    sum_inv2 = np.sum(1.0 / n**2) - 1.0 / N**2
    sum_inv3 = np.sum(1.0 / n**3) - 1.0 / N**3

    kappa3 = sum_inv2 / sum_inv**(3/2)
    kappa4 = sum_inv3 / sum_inv**2

    return kappa3, kappa4
```

### 2.3 기존 Planner 문서 업데이트

| 문서 | 수정 내용 |
|------|-----------|
| `planner/saddlepoint_pag_alignment_plan.md` | 이슈 B 상태: "Paper B에서 E[Ri^m] 공식 확보됨" 추가 |
| `planner/saddlepoint_pag_alignment_plan.md` | 이슈 A 주석: "두 버전 모두 +1/2, 수학적으로 -0.5와 동치 해석 가능" |
| `planner/saddlepoint_pag_alignment_plan.md` | 단계 3 (모멘트 엔진): "3-1안 채택 — Lemma 2.5 기반" |

### 2.4 테스트 수정 사항

| 항목 | 내용 |
|------|------|
| 신규 | `test_rank_sum_moments.py` — E[Ri^m] 공식을 brute-force 열거와 대조 검증 |
| 신규 | `test_cross_moments.py` — E[Ri^2·Rj^2] 공식을 열거 대조 |
| 수정 | `test_approximations.py` — Var(H) Paper 공식 경로 검증 추가 |
| 신규 | `test_asymptotic_cumulants.py` — 식 13-14를 simulation과 대조 |

---

## Phase 3: 판별 보고서 (Judge Output)

### 3.1 교차 검증 결과

| 기획서 항목 | 생성자 반영 | 판정 | 비고 |
|------------|-----------|------|------|
| #10 Theorem 2.3 (μ3, μ4) | G4, G5 | **CONDITIONAL** | E[Ri^6], E[Ri^8] 공식 미제시 — 추가 유도 필요 |
| #11 Theorem 2.4 (A 기반) | 참조만 | **PASS** | Theorem 2.3 경로 우선, 2.4는 유도 근거로만 사용 |
| #12 Lemma 2.5 (E[Ri^m]) | G1, G2 | **PASS** | m=1..4 + cross-moment 완전 커버 |
| #14 식 13-14 (asymp κ) | G6 | **PASS** | closed-form 함수 추가 |
| CC 방향 | G7 | **PASS** | 주석 보강으로 충분 |

### 3.2 누락 항목 점검

#### 누락 #1: E[Ri^6] 유도 필요 (μ3 계산용)
- **문제**: Theorem 2.3의 μ3에 E[Ri^6]이 필요하지만, Lemma 2.5는 m=4까지만 제공
- **해결 방안**:
  - A) Sampling without replacement의 6차 모멘트를 Kendall & Stuart(1979) 또는 combinatorial formula로 유도
  - B) Paper B Theorem 2.4의 E[A^3] 경로: A^3 전개에서 E[Ri^6] 대신 E[Ri^4]·E[Rj^2] 등 교차모멘트 경로 활용
  - C) Simulation fallback 유지 (현행)
- **권장**: 당장은 C(simulation), 향후 A 또는 B로 전환

#### 누락 #2: E[Ri^4·Rj^2] 고차 교차모멘트
- Theorem 2.4의 E[A^3] 계산에 필요하지만, 논문에 공식 미제시
- **권장**: Theorem 2.3 경로보다 simulation이 현실적

#### 누락 #3: Gram-Charlier (GC-A) 별도 구현 확인
- Paper B의 GC-A는 edgeworth.py의 `cdf_normal_based` 함수가 담당
- 현재 구현은 **Paper B Section 3.4와 일치** — 변경 불필요

### 3.3 구현 정합성 체크

| 코드 파일 | Paper B 섹션 | 현재 정합성 | 반영 후 예상 |
|-----------|-------------|-----------|------------|
| `moments.py` (exact) | - | **정확** (열거 기반) | 변경 없음 |
| `moments.py` (simulation) | - | **정확** (MC) | 변경 없음 |
| `moments.py` (asymptotic) | Thm 2.2, 2.3, Lemma 2.5 | **불일치** (chi-square 근사) | G1-G5 반영 시 **일치** |
| `edgeworth.py` | Section 3.3, 3.4 | **일치** | G6 추가 시 보강 |
| `saddlepoint.py` | Section 3.2 | **일치** (CC 방향 이슈 유지) | G7 주석 보강 |
| `pam.py` | Section 3.1 | **완전 일치** | 변경 없음 |

### 3.4 최종 판정

```
┌─────────────────────────────────────────┐
│         CONDITIONAL_PASS                │
│                                         │
│  사유:                                   │
│  1. G1, G2, G6: 즉시 구현 가능 — PASS    │
│  2. G3: 즉시 구현 가능 — PASS            │
│  3. G4-G5: E[Ri^6]/E[Ri^8] 미확보       │
│     → simulation fallback 유지 조건부 PASS│
│  4. G7: 주석 보강 충분 — PASS            │
│                                         │
│  조건:                                   │
│  - G4의 asymptotic μ3/μ4 교체는           │
│    E[Ri^6] 유도 완료 후 진행              │
│  - 그 전까지 simulation 경로 유지         │
└─────────────────────────────────────────┘
```

---

## Phase 4: 실행 로드맵

### Stage 1: 즉시 실행 가능 (E[Ri^m] 빌딩블록)
1. `moments.py`에 `_rank_sum_moment(ni, N, m)` 추가 (G1)
2. `moments.py`에 `_cross_moment_R2R2(ni, nj, N)` 추가 (G2)
3. `moments.py`에 `_compute_variance_paper()` 추가 및 검증 (G3)
4. 단위테스트: brute-force 열거와 대조

### Stage 2: Asymptotic cumulant 추가
5. `edgeworth.py` 또는 `moments.py`에 `asymptotic_cumulants()` 추가 (G6)
6. 대표본 설계에서 simulation 결과와 대조 검증

### Stage 3: Asymptotic μ3/μ4 교체 (E[Ri^6] 확보 후)
7. E[Ri^6] 공식 유도 (Kendall & Stuart 참조 또는 combinatorial 직접 계산)
8. `_compute_central_moment_3_asymptotic()` 교체 (G4)
9. `_compute_central_moment_4_asymptotic()` 교체 (G4)
10. 회귀 테스트: 모든 테이블 재현 확인

### Stage 4: 문서/주석 정리
11. `saddlepoint.py` CC 주석 보강 (G7)
12. `planner/saddlepoint_pag_alignment_plan.md` 업데이트

---

## Appendix: E[Ri^m] 고차 모멘트 유도 노트

Sampling without replacement에서 Ri = Σ_{j=1}^{ni} X_j, X_j ∈ {1,...,N}일 때:

E[Ri^m]은 power sum의 기대값이며, m ≤ 4까지는 Lemma 2.5에서 확보.
m = 6, 8의 경우:

- **방법 A**: Symmetric function identities. E[Ri^m] = Σ (multinomial products of power sums) × (combinatorial coefficients)
- **방법 B**: David & Barton (1962), Kendall & Stuart (1979)의 표 참조
- **방법 C**: Exact 열거가 가능한 소표본에서 regression으로 공식 계수 추정 후 검증

Paper B 자체도 "E[Ri^6] and E[Ri^8] can be computed using standard formulas for sums of powers of integers and properties of sampling without replacement" (p.10 상단)으로 기술하며 구체적 closed-form은 제시하지 않음.

**결론**: E[Ri^6]은 유도 가능하나 복잡. 당분간 simulation 경로가 실용적.
