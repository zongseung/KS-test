# Murakami, Lee & Ha (Kruskal–Wallis) 재현 스크립트 수정 계획서 (단일 MD)

대상: `examples/reproduce_paper_tables.py` (및 관련 util 함수/클래스)  
목표: Random designs 파트에서 발생하는 **self-calibration(자기보정) 구조**를 제거하고, 논문 의도대로 **Reference(Exact/Simulation) 임계값** 기준으로 모든 근사법 tail을 비교하도록 수정한다.

---

## 0. 배경 및 문제 정의

### 현상
- Random 3/4/5-group designs 출력에서 `PAG(6)`(=pam6)가 **항상 0.100000**으로 찍히는 케이스가 반복적으로 발생.

### 원인(핵심)
- Random designs에서 `H(10%)`를 `pam6.critical_value(0.10)`로 설정함.
- 그러면 같은 `H`에 대해 `pam6.tail_probability(H)`는 정의상 `≈ 0.10`으로 귀결되어, **PAG(6) 비교열이 정보가 없어짐(비교 불가능)**.

### 논문 의도(수정 방향)
- 임계값 `H(α)`는 **reference 분포**로부터 결정되어야 한다.
  - 작은 N: Exact distribution
  - 큰 N: Monte Carlo simulation
- 이후 그 **동일한 H**에서 CHI/SD/ED/PAG 등을 비교한다.

---

## 1. 설계 원칙(Reference 정의)

### 1.1 Reference 임계값 우선순위
Random designs에서 `H(α)`는 아래 규칙으로 정한다.

1) **Exact critical value** (가능할 때)  
2) **Monte Carlo simulation critical value** (Exact 불가/비현실적일 때)  
3) **Chi-square critical value** (MC도 실패/금지 등 최후의 fallback)

> 중요: Random designs에서 “PAM6 임계값”은 reference로 절대 쓰지 않는다.

### 1.2 Exact 가능 여부 기준(논문/실무 기준)
- k=3: N <= 15 정도는 exact 가능
- k=4: N <= 13 정도로 더 엄격
- k>=5: N <= 10 정도로 매우 엄격

※ 위 기준은 현재 스크립트에 있는 limit 로직을 **reference 임계값 산출에도 동일하게 적용**한다.

---

## 2. 유틸/API 구조 개편

현재: `get_reference_probability(approx, h, N, k, ...)`가 “주어진 h에서 tail만” 반환  
필요: Random designs는 “h(임계값) 자체”가 reference여야 하므로 **reference critical value 함수**가 필요

---

---
name: refactor-utils
description: Random designs에서 reference 임계값을 일관되게 얻기 위해 get_reference_critical_value()를 추가하고, 기존 get_reference_probability()의 역할을 명확히 분리합니다.
---

## 2.1 신규 함수 추가: get_reference_critical_value()
**추가**: `get_reference_critical_value(sample_sizes, alpha, n_simulations=10000, seed=None)`

- 입력:
  - `sample_sizes: list[int]`
  - `alpha: float`
  - `n_simulations: int = 10000`
  - `seed: Optional[int]` (재현성 확보용)
- 출력:
  - `h_ref: float` (reference critical value)
  - `ref_alpha: float` (실제 tail; exact면 정확, sim이면 MC 오차 반영)
  - `source: str` in {"exact", "simulation", "chi_square"}

### 의사코드
- k=len(sample_sizes), N=sum(sample_sizes)
- limit_N = (k==3 ? 15 : k==4 ? 13 : 10)
- if N <= limit_N:
  - try exact critical value: `(h_ref, ref_alpha) = ExactDistribution(...).critical_value(alpha)`
  - return (h_ref, ref_alpha, "exact")
- else:
  - try simulation critical value: `(h_ref, ref_alpha) = MonteCarloSimulation(...).critical_value(alpha)`
  - return (h_ref, ref_alpha, "simulation")
- fallback:
  - h_ref = chi2.ppf(1-alpha, df=k-1)
  - ref_alpha = 1 - chi2.cdf(h_ref, df=k-1)
  - return (h_ref, ref_alpha, "chi_square")

## 2.2 기존 함수 정리(선택)
- `get_reference_probability()`는 **paper-table(4.2~4.3 등)**의 “E-P column(특정 H에서 reference tail)” 계산에만 쓰고,
- Random designs 흐름에서는 **사용하지 않도록** 정리한다.

## 2.3 seed 정책(권장)
- Random designs는 출력이 변동되므로, `seed`를 옵션으로 두고
  - 스크립트 상단에서 `SEED = 20260202` 같은 상수로 고정(또는 CLI 인자)
  - `random.seed(SEED)`, `np.random.seed(SEED)`
  - `MonteCarloSimulation(..., seed=SEED)` 적용을 권장한다.

---

## 3. Random designs 출력 로직 수정(핵심 변경)

대상 함수:
- `generate_random_three_group_designs()`
- `generate_random_four_group_designs()`
- `generate_random_k_group_designs()`
- `comprehensive_random_study()` 내부 루프

---

---
name: fix-random-designs
description: Random designs에서 H(10%)를 PAM6 임계값으로 두는 self-calibration을 제거하고, reference(Exact/Simulation) 임계값으로 교체합니다.
---

## 3.1 변경 전(현재)
- H(10%) = pam6.critical_value(0.10) (또는 N 작으면 exact)
- E-P는 exact 또는 simulation tail로 계산
- 결과적으로 PAG(6)=0.100000 고정이 빈번

## 3.2 변경 후(목표)
- `H(10%)`를 항상 `get_reference_critical_value()`로 결정
- 그 `H_ref`에서 모든 근사법 tail을 계산해 비교

### 적용 방식(예시)
- 기존:
  - `cv = approx.critical_value(0.10, 'pam6')`
- 수정:
  - `H_ref, ref_alpha, src = get_reference_critical_value(sample_sizes, 0.10, n_simulations=10000, seed=SEED)`

## 3.3 테이블 컬럼 정의 개선(권장)
Random designs 테이블에서 E-P 컬럼 의미 명확화:

- `E-P` 컬럼 = `ref_alpha`
- simulation 기반이면 `*` 마킹
- `SRC` 컬럼 추가(권장): `exact / sim / chi` 중 하나 출력

예시 행:
- `E-P = 0.104400*` + `SRC=sim`
- `E-P = 0.101587` + `SRC=exact`

## 3.4 적용 범위(빠짐 없이)
- Random 3-group, Random 4-group, Random k-group(>=5), comprehensive study 모두 동일 규칙 적용
- “H(10%) 산출” 라인이 pam6로 되어 있는 모든 곳을 reference로 교체

---

## 4. Monte Carlo 임계값 정의(인덱싱/부등호) 점검

---

---
name: explain-code
description: MonteCarloSimulation.critical_value()의 인덱싱/부등호 정의를 점검하고, 임계값(quantile) 산출이 alpha를 제대로 만족하는지 설명 및 수정 지침을 제공합니다.
---

## 4.1 현재 구현(점검 대상)
```python
idx = int(np.ceil((1 - alpha) * self.n_simulations)) - 1
cv = self._H_values[idx]
actual_alpha = np.mean(self._H_values >= cv)
