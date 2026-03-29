# 논문 비교 기획서: Paper A (revision) → Paper B (final_revision)

> **Paper A**: `Kruskal_Wallis_Test_revisionn.pdf` (29 pages)
> **Paper B**: `Kruskal_Wallis_Test_final_revision.pdf` (30 pages)
> **분석일**: 2026-03-30

---

## 1. 변경 요약

| # | 영역 | 유형 | 변경 내용 요약 | 코드 영향도 |
|---|------|------|---------------|-----------|
| 1 | 제목/포맷 | STRUCTURE | 타이틀 페이지 → 학술지 투고 포맷(Elsevier) 전환 | 없음 |
| 2 | 저자 순서 | AUTHOR | Lee, Murakami, Ha → **Murakami, Lee, Ha** (1저자 변경) | 없음 |
| 3 | 저자 소속 표기 | AUTHOR | a,b,c 소속 재배열 (Murakami가 a→첫째, Ha가 a,c→b,c) | 없음 |
| 4 | Abstract | TEXT | "based on moments with polynomial correction terms" 표현 추가; Hall(1992,1993), Kolassa(1995) Edgeworth 비교 언급 **신규 추가** | 없음 |
| 5 | Section 2 제목 | STRUCTURE | "Exact and Asymptotic Theory" → "Kruskal-Wallis Statistic" (간결화) | 없음 |
| 6 | Section 2.1 | TEXT | PGF 서술 보강 — P(z) 정의를 coefficient counting 관점으로 재작성, Diaconis/Gangolli(1977), Anderson(2001), Stanley(1999) 참고문헌 추가 | 없음 |
| 7 | Section 2.2 | TEXT | Iman et al.(1975) 재귀식에 Hollander/Wolfe(1999), Conover(1999) 참고 추가 | 없음 |
| 8 | Section 2.3 | STRUCTURE | "Analytic geometry and asymptotic convergence" → "Exact representation"으로 제목 변경; Section 2.3.1 하위절 제거 → 본문 통합 | 없음 |
| 9 | **Theorem 2.1~2.2** | FORMULA | 기존과 동일 (변경 없음) | 없음 |
| 10 | **Theorem 2.3 (μ3, μ4)** | **FORMULA** | **Paper A에 없던 내용 — Paper B에서 새로 추가**: μ3 = (12^3/N^3(N+1)^3) Σ E(Ri^6)/ni^3 - 3E(H)Var(H) - μ1^3, μ4 유사 형태 | **높음** |
| 11 | **Theorem 2.4 (First four raw moments)** | **FORMULA** | **Paper B에서 새로 추가**: A = Σ Ri^2/ni 정의 후 μ1~μ4를 E[A^m] 기반 binomial expansion으로 표현 | **높음** |
| 12 | **Lemma 2.5 (A^m expansion)** | **FORMULA** | **Paper B에서 새로 추가**: A^2, A^3, A^4의 cross-moment 전개식 + E[Ri^m] 공식 (m=1,2,3,4) + E[Ri^2 Rj^2] cross-moment 공식 | **매우 높음** |
| 13 | Section 3.2 (Saddlepoint) | TEXT | 연속성 보정 서술에 "as shown in Butler (2007)" 참조 명시 추가 | 낮음 |
| 14 | Section 3.3 (Edgeworth) | **METHOD** | **Paper B에서 대폭 확장**: Hall(1992,1993), Kolassa(1995)의 asymptotic cumulant 기반 Edgeworth 전개 명시적 추가. **κ3, κ4의 asymptotic 공식 (식 13, 14)** 신규 추가 | **높음** |
| 15 | Section 3.4 (Gram-Charlier) | TEXT | 동일 (변경 없음) | 없음 |
| 16 | Section 4 (Tables) | TABLE | 테이블 번호/구성 동일, **수치값 동일** | 없음 |
| 17 | Section 4 서술 | TEXT | "Tables 4.2-4.3" → 단순 서술 변경 | 없음 |
| 18 | References | REFERENCE | **7개 참고문헌 신규 추가**: Hall(1992), Hall(1993a), Hall(1993b), Kolassa(1995), Bickel & van Zwet(1976), Andrews(1976), + 기타 조합론 관련 | 낮음 |

---

## 2. 핵심 변경 상세

### 변경 #10: Theorem 2.3 — μ3, μ4 raw moment 공식 (신규)

- **Paper A**: 이 정리 자체가 존재하지 않음. Section 2에서 μ1 = k-1, Var(H) 공식만 제시하고, μ3/μ4에 대한 명시적 공식 없음
- **Paper B** (p.9, 식 8-9):
  ```
  μ3 = E(H^3) = (12^3 / N^3(N+1)^3) Σ E(Ri^6)/ni^3 - 3·E(H)·Var(H) - μ1^3
  μ4 = E(H^4) = (12^4 / N^4(N+1)^4) Σ E(Ri^8)/ni^4 - 4·μ1·μ3 - 6·μ1^2·Var(H) - μ1^4
  ```
- **차이**: Paper B는 3차/4차 raw moment를 E(Ri^6), E(Ri^8) 기반으로 **명시적으로** 유도. 이전 버전에서는 이 공식이 논문에 없어 구현 시 추측에 의존해야 했음
- **영향 코드**: `kw_approx/moments.py` — μ3, μ4 계산 로직
- **기존 이슈 연결**: 이슈 B (모멘트 근사가 논문 의도와 다를 가능성)

### 변경 #11: Theorem 2.4 — A 기반 raw moment 전개 (신규)

- **Paper A**: 없음
- **Paper B** (p.10): A := Σ Ri^2/ni 로 정의하고, H = 12/(N(N+1))·A - 3(N+1)의 관계를 이용하여:
  ```
  μ1 = (12/N(N+1)) Σ E[Ri^2]/ni - 3(N+1)
  μ2 = (12/N(N+1))^2 E[A^2] - 2·3(N+1)·(12/N(N+1))·E[A] + [3(N+1)]^2
  μ3 = (12/N(N+1))^3 E[A^3] - 3·3(N+1)·(12/N(N+1))^2·E[A^2]
       + 3·[3(N+1)]^2·(12/N(N+1))·E[A] - [3(N+1)]^3
  μ4 = (유사 binomial expansion)
  ```
- **차이**: binomial expansion을 통해 μ1~μ4를 E[A], E[A^2], E[A^3], E[A^4]로 환원 — 이전에는 이 중간 단계가 논문에 없었음
- **영향 코드**: `kw_approx/moments.py` — 전체 moment 계산 파이프라인

### 변경 #12: Lemma 2.5 — A^m의 multinomial 전개 + Ri moment 공식 (신규, 핵심)

- **Paper A**: 없음
- **Paper B** (p.11): A^m (m=2,3,4)의 cross-product 전개:
  ```
  A^2 = Σ Ri^4/ni^2 + 2·Σ_{i<j} Ri^2·Rj^2/(ni·nj)
  A^3 = Σ Ri^6/ni^3 + 3·Σ_{i≠j} Ri^4·Rj^2/(ni^2·nj) + 6·Σ_{i<j<l} Ri^2·Rj^2·Rl^2/(ni·nj·nl)
  A^4 = (유사, 4중 cross-product까지)
  ```
  그리고 Ri의 개별 moment:
  ```
  E[Ri]   = ni·(N+1)/2
  E[Ri^2] = ni·(N-ni)(N+1)/12 + ni^2·(N+1)^2/4
  E[Ri^3] = ni^2·(N+1)^3/8 + ni·(N-ni)(N+1)(N+2)/24
  E[Ri^4] = ni·(N-ni)(N+1)(2N+1)/30 + ni^2·(N+1)^2·(N-ni)/12 + ni^4·(N+1)^4/16
  ```
  Cross-moment:
  ```
  E[Ri^2·Rj^2] = E[Ri^2]·E[Rj^2] - ni·nj·(N+1)^2·(N-2)(N-ni-nj) / (12·(N-1)·(N-3))
  ```
- **차이**: 이것이 **가장 중대한 추가**. 코드 구현의 moment 엔진이 의존하는 핵심 수식들이 Paper A에는 전혀 없었고, Paper B에서 처음으로 명시됨
- **영향 코드**: `kw_approx/moments.py` 전체
- **기존 이슈 연결**: 이슈 B (모멘트 근사), 이슈 C (테스트 기준 확립)

### 변경 #14: Edgeworth expansion에 asymptotic cumulant 공식 추가

- **Paper A** (p.16): Edgeworth 전개는 일반 형태만 제시. κ3, κ4의 구체적 공식 없음
- **Paper B** (p.16, 식 13-14): Hall/Kolassa의 asymptotic cumulant 공식 명시:
  ```
  κ3 = (Σ 1/ni^2 - 1/N^2) / (Σ 1/ni - 1/N)^{3/2}    (식 13)
  κ4 = (Σ 1/ni^3 - 1/N^3) / (Σ 1/ni - 1/N)^2          (식 14)
  ```
- **차이**: 이 공식은 standardized H의 asymptotic skewness/kurtosis를 ni와 N만으로 closed-form으로 표현. Paper A에서는 이 공식이 없어 exact cumulant만 사용 가능했음
- **영향 코드**: `kw_approx/edgeworth.py` — asymptotic cumulant 경로 추가 가능
- **기존 이슈 연결**: 이슈 B와 간접 관련 (대표본 근사의 이론적 근거)

---

## 3. 새로 추가된 내용 (Paper B only)

| # | 내용 | 위치 | 코드 영향 |
|---|------|------|---------|
| 1 | **Theorem 2.3**: μ3, μ4 raw moment 공식 (E[Ri^6], E[Ri^8] 기반) | p.9 | `moments.py` |
| 2 | **Theorem 2.4**: A 기반 μ1~μ4 binomial expansion | p.10 | `moments.py` |
| 3 | **Lemma 2.5**: A^m 전개 + E[Ri^m] (m=1..4) + E[Ri^2·Rj^2] cross-moment | p.11 | `moments.py` |
| 4 | **식 13-14**: Asymptotic κ3, κ4 closed-form (Hall/Kolassa) | p.16 | `edgeworth.py` |
| 5 | Hall(1992,1993), Kolassa(1995) Edgeworth 비교 프레이밍 | Abstract, Section 3.3 | 없음 |
| 6 | Diaconis/Gangolli(1977), Anderson(2001), Stanley(1999), Andrews(1976) 등 조합론 참고문헌 | Section 2.1, References | 없음 |
| 7 | Bickel & van Zwet(1976) 참고문헌 | References | 없음 |
| 8 | Var(H) exact formula with E[Ri^4], E[Ri^2·Rj^2] 명시 | p.9 (식 6-7) | `moments.py` |

---

## 4. 삭제/축소된 내용 (Paper A only)

| # | 내용 | Paper A 위치 | 비고 |
|---|------|-------------|------|
| 1 | Section 2.3.1 "Exact finite-sample representation" 하위절 구조 | p.9 | Paper B에서 본문 통합, 내용은 보존 |
| 2 | "This perspective reveals a rich geometric and algebraic structure..." 문단 | p.3 | Paper B에서 축약 |
| 3 | 타이틀 페이지의 별도 소속/이메일 레이아웃 | p.1 | Elsevier 포맷으로 대체 |

> **실질적 삭제는 없음** — Paper A의 모든 핵심 내용은 Paper B에 보존되고, 추가만 이루어짐.

---

## 5. 연속성 보정 관련 확인

- **Paper A** (p.14): `v + 1/2` (연속성 보정에서 `w_cc = sgn(t̂)√(2(t̂(v+1/2) - K(t̂)))`)
- **Paper B** (p.14): 동일 — `v + 1/2`
- **코드 현황**: `saddlepoint.py`에서 `v - 0.5` 사용 중
- **결론**: 두 논문 버전 모두 `v + 1/2` 표기. **코드의 `-0.5`와 논문의 `+1/2` 불일치는 여전히 존재**
  - 단, 기존 분석(saddlepoint_pag_alignment_plan.md)에서 `-0.5`가 테이블 수치 재현에 더 근접한다는 결과가 있었음

---

## 6. 테이블 수치 비교

두 논문의 테이블 수치는 **완전히 동일**:

| 테이블 | Paper A | Paper B | 일치 |
|--------|---------|---------|------|
| Table 1 (3,3,3, α=0.10) | 동일 | 동일 | O |
| Table 2 (3그룹 balanced, α=0.10) | 동일 | 동일 | O |
| Table 3 (3그룹 balanced, α=0.05) | 동일 | 동일 | O |
| Table 4-5 (3그룹 larger n) | 동일 | 동일 | O |
| Table 6 (4그룹 (3,2,2,5)) | 동일 | 동일 | O |
| Table 7-8 (4그룹 추가) | 동일 | 동일 | O |

---

## 7. 참고문헌 비교

### Paper B에서 추가된 참고문헌 (7개)
1. **Andrews, G.E. (1976)**. *The Theory of Partitions*. — 조합론 기초
2. **Anderson, C.W. (2001)**. *An Introduction to Multivariate Statistical Analysis*. — 다변량 통계
3. **Bickel, P.J. and van Zwet, W.R. (1976)**. Asymptotic expansions for the power of distribution-free tests. — Edgeworth 이론적 근거
4. **Conover, W.J. (1999)**. *Practical Nonparametric Statistics*. — exact 계산 참조
5. **Diaconis, P. and Gangolli, A. (1977)**. Rectangular arrays with fixed margins. — PGF 이론
6. **Hall, P. (1992)**. *The Bootstrap and Edgeworth Expansion*. — Edgeworth 비교 대상
7. **Hall, P. (1993a)**. On Edgeworth expansion and bootstrap confidence bands. — Edgeworth 비교
8. **Hall, P. (1993b)**. Edgeworth expansions for studentized statistics under weak assumptions. — Edgeworth 비교
9. **Hollander, M. and Wolfe, D.A. (1999)**. *Nonparametric Statistical Methods*. — exact 계산
10. **Kendall, M.G. and Stuart, A. (1979)**. *The Advanced Theory of Statistics*. — moment 공식
11. **Kolassa, J.E. (1995)**. Edgeworth approximations for rank-sum test statistics. — 핵심 비교 대상
12. **Stanley, R.P. (1999)**. *Enumerative Combinatorics, Volume 2*. — q-binomial

### Paper A에서 삭제된 참고문헌
- 없음 (Paper A의 모든 참고문헌은 Paper B에도 존재)

---

## 8. 코드 영향 매핑

### 즉시 반영 필요 (코드 영향도: 높음~매우높음)

| 우선순위 | 변경 | 영향 파일 | 조치 |
|----------|------|---------|------|
| **P1** | Lemma 2.5: E[Ri^3], E[Ri^4] 명시적 공식 | `kw_approx/moments.py` | 기존 구현과 대조 → 불일치 시 교체 |
| **P2** | Lemma 2.5: E[Ri^2·Rj^2] cross-moment | `kw_approx/moments.py` | Var(H) 정확도 검증, 교차 모멘트 구현 확인 |
| **P3** | Theorem 2.3: μ3, μ4 공식 | `kw_approx/moments.py` | 3차/4차 raw moment 계산 경로 검증 |
| **P4** | 식 13-14: Asymptotic κ3, κ4 | `kw_approx/edgeworth.py` | 대표본 경로에 asymptotic cumulant 추가 가능 |

### 참고 확인 (코드 영향도: 낮음)

| 변경 | 영향 파일 | 비고 |
|------|---------|------|
| CC 방향 (v+1/2) | `kw_approx/saddlepoint.py` | 두 버전 모두 동일, 기존 이슈 A 유지 |
| Butler(2007) 참조 명시 | — | 코드 변경 불필요, 문서만 업데이트 |

---

## 9. 기존 이슈와의 연결

| 기존 이슈 | 관련 Paper B 변경 | 해소 여부 |
|----------|-----------------|---------|
| **이슈 A** (CC 부호 불일치) | 변경 없음 (v+1/2 유지) | 미해소 — 여전히 코드(-0.5)와 논문(+1/2) 충돌 |
| **이슈 B** (모멘트 근사 부정확) | **Theorem 2.3, 2.4, Lemma 2.5로 명시적 공식 확보** | **부분 해소** — 구현 정답 기준이 확립됨 |
| **이슈 C** (테스트 느슨) | 수치 테이블 변경 없음 | 미해소 |
| **이슈 D** (Gamma saddlepoint 단순화) | 변경 없음 | 미해소 |
| **이슈 E** (ER1 tail 비단조) | 변경 없음 | 미해소 |
| **이슈 F** (PAG(6) 불안정) | 변경 없음 | 미해소 |
| **이슈 G** (PAG pdf/cdf 불일치) | 변경 없음 | 미해소 |

---

## 10. 검증 기준 (판별자 입력)

- [x] 모든 FORMULA 변경이 식별되었는가 — **Yes** (Theorem 2.3, 2.4, Lemma 2.5, 식 13-14)
- [x] 코드 영향 매핑이 정확한가 — `moments.py`, `edgeworth.py` 중심
- [x] 기존 planner 이슈와의 연결이 완전한가 — 7개 이슈 모두 매핑 완료
- [x] 새로 추가된 내용이 누락 없이 열거되었는가 — 8개 항목
- [x] 삭제된 내용 확인 — 실질 삭제 없음 (구조 변경만)
- [x] 테이블 수치 일치 확인 — 8개 테이블 모두 동일

---

## 11. 요약 및 다음 단계

### 핵심 발견
**Paper B(final_revision)의 가장 중요한 변경은 Section 2에 추가된 3개의 새로운 수학적 결과(Theorem 2.3, 2.4 + Lemma 2.5)이다.** 이 결과들은:

1. **μ3, μ4의 명시적 공식**을 E[Ri^6], E[Ri^8] 기반으로 제시
2. **A = Σ Ri^2/ni 기반 binomial expansion**으로 μ1~μ4 계산의 체계적 경로를 확립
3. **E[Ri^m] (m=1..4)과 E[Ri^2·Rj^2] cross-moment**의 closed-form을 제시

이는 기존 코드의 moment 엔진(`moments.py`)이 의존하는 수식들의 **정답 기준(ground truth)**을 처음으로 논문 수준에서 확보한 것이다.

추가로, **Edgeworth 섹션에 Hall/Kolassa의 asymptotic cumulant (식 13-14)**가 추가되어, 대표본 경로에서의 이론적 근거가 강화되었다.

### 다음 단계
1. **`/he:generate`**: 이 기획서를 기반으로 반영 계획서 작성
   - P1~P3: `moments.py`의 E[Ri^m], cross-moment, μ3/μ4 계산이 Lemma 2.5/Theorem 2.3과 일치하는지 코드 수준 검증
   - P4: `edgeworth.py`에 asymptotic cumulant 경로 추가 여부 결정
2. **`/he:judge`**: 반영 계획서의 완전성/정합성 판정
