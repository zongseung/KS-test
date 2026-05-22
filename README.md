[**English**](README_EN.md) | **Korean**

# kw-approx: Higher Order Asymptotic Approximations for Kruskal-Wallis Statistics

Kruskal-Wallis 검정 통계량의 고차 점근 근사를 구현한 Python 패키지입니다.

## 논문 정보

> **Lee, J.-S., Murakami, H., & Ha, H.-T. (2026).** *Higher Order Asymptotic Approximations of Kruskal-Wallis Statistics Based on Skewness and Kurtosis.* Preprint submitted to Journals.

## 배경

Kruskal-Wallis 검정은 k개 집단의 위치 모수 동일성을 검정하는 비모수적 방법으로, 일원분산분석(One-way ANOVA)의 순위 기반 대안입니다.

### Kruskal-Wallis H 통계량

$$H = \frac{12}{N(N+1)} \sum_{i=1}^{k} \frac{R_i^2}{n_i} - 3(N+1)$$

- $N = \sum_{i=1}^{k} n_i$: 전체 표본 크기
- $n_i$: $i$번째 집단의 표본 크기
- $R_i$: $i$번째 집단의 순위합

전통적으로 $H$의 귀무분포는 자유도 $k-1$인 카이제곱 분포로 근사하지만, 소표본에서는 정확도가 떨어집니다. 이 패키지는 새들포인트, Edgeworth, Gram-Charlier, PAG 근사를 통해 더 정확한 근사를 제공합니다.

## 프로젝트 구조

```
KS-test/
|
|-- kw_approx/                    # 핵심 패키지
|   |-- __init__.py               #   패키지 초기화 (v0.3.0)
|   |-- kruskal_wallis.py         #   H 통계량 계산
|   |-- moments.py                #   모멘트/큐뮬런트 (exact/simulation)
|   |-- saddlepoint.py            #   새들포인트 근사 (ER1, ER2, KT + L-R)
|   |-- edgeworth.py              #   Edgeworth 전개 (chi-sq 기반 Laguerre)
|   |-- gram_charlier.py          #   Gram-Charlier Type A (Hermite)
|   |-- pam.py                    #   PAG(d) 다항식 조정 감마 근사
|   |-- exact.py                  #   정확 분포 (재귀 열거, 소표본)
|   |-- simulation.py             #   Monte Carlo 시뮬레이션
|   +-- approximator.py           #   통합 인터페이스 (KWApproximator)
|
|-- examples/
|   +-- reproduce_paper_tables.py #   논문 Table 4.1~4.7 재현 + 확장
|
|-- tests/
|   +-- test_approximations.py    #   36개 테스트 케이스
|
|
|-- README.md
|-- README_EN.md
|-- CLAUDE.md
|-- pyproject.toml
+-- Kruskal_Wallis_Test.pdf       # 논문 원본
```

## 핵심 파이프라인

```
                         +------------------+
                         |   sample_sizes   |
                         | (n1, n2, ..., nk)|
                         +--------+---------+
                                  |
                                  v
                     +------------+-------------+
                     |        KWMoments          |
                     |  (exact enum / MC sim)    |
                     |                           |
                     |  k<=3, N<=15: exact       |
                     |  k>=4, N<=13: exact       |
                     |  otherwise:   simulation  |
                     +------------+-------------+
                                  |
                                  v
                       +--------------------+
                       | Cumulants k1...k4  |
                       | (exact finite-     |
                       |  sample, NOT       |
                       |  asymptotic chi-sq)|
                       +----+----+----+----+
                            |    |    |    |
            +---------------+    |    |    +----------------+
            |           +--------+    +---------+           |
            v           v                       v           v
     +-----------+ +-----------+         +-----------+ +-----------+
     |Saddlepoint| |Edgeworth  |         |Gram-      | |   PAG(d)  |
     |  ER1/ER2  | |(chi-sq +  |         |Charlier   | | (gamma x  |
     |    KT     | | Laguerre) |         |Type A     | | polynomial|
     +-----------+ +-----------+         |(Hermite)  | |  d=4,6)   |
     |  +- CC    |                       +-----------+ +-----------+
     |  variants |
     +-----+-----+
           |
           v
    +-------------+
    |Lugannani-   |
    |Rice tail    |---------->  P(H >= v)
    |probability  |
    +-------------+
```

## 설치

```bash
pip install numpy scipy
pip install -e .
```

## 빠른 시작

```python
from kw_approx import KWApproximator

# 3개 집단, 각 3명씩
approx = KWApproximator([3, 3, 3])

# H = 4.62에서 꼬리확률 P(H >= 4.62)
p_value = approx.tail_probability(4.62, method='ER1')
print(f"ER1 P-value: {p_value:.6f}")

# 여러 방법 비교
results = approx.compare_methods(4.62)
for method, p in results.items():
    print(f"{method}: {p:.6f}")
```

## 구현된 근사 방법

### 방법 일람 (논문 <-> 코드)

```
+----------+----------------+---------------------+-----------------+
| 논문     | 코드 method    | CGF                 | 비고            |
+----------+----------------+---------------------+-----------------+
| SD1      | 'ER1'          | Easton-Ronchetti 1  | Saddlepoint     |
| SD2(KT)  | 'KT'           | Kakizawa-Taniguchi  | Saddlepoint     |
| SDC1     | 'ER1_cc'       | ER1 + CC            | 연속성 보정     |
| SDC2(KT) | 'KT_cc'        | KT + CC             | 연속성 보정     |
| CHI      | 'chi_square'   | -                   | Baseline        |
| ED       | 'edgeworth'    | -                   | Laguerre 기반   |
| GC-A     | 'gram_charlier'| -                   | Hermite 기반    |
| PAG(4)   | 'pam'          | -                   | 감마 x 다항식   |
| PAG(6)   | 'pam6'         | -                   | 감마 x 다항식   |
| -        | 'exact'        | -                   | 소표본 전용     |
| -        | 'simulation'   | -                   | MC 시뮬레이션   |
+----------+----------------+---------------------+-----------------+
```

> **참고:** 논문 테이블의 SD2/SDC2 컬럼은 Wang damped CGF에서 **Kakizawa-Taniguchi(KT)** CGF로 교체되었습니다 (헤더 표기 `SD2(KT)` / `SDC2(KT)`). Wang은 N≈15 부근에서 K''(t̂)가 0에 가까워져 Lugannani-Rice 보정항이 폭주하는 수치 불안정이 있어, 다항식 CGF로 안정적인 KT로 대체했습니다. Wang CGF 구현은 삭제하지 않고 코드 주석으로 보존했지만, 현재 실행 경로에서는 비활성화되어 있습니다.

### 큐뮬런트 (Cumulants) -- Exact Finite-Sample

논문 Section 4에 따라, 모든 근사 방법은 **정확한 유한표본 큐뮬런트**를 사용합니다:

- 소표본 (k<=3, N<=15 또는 k>=4, N<=13): exact distribution으로 계산
- 대표본: Monte Carlo simulation (50,000 iterations, deterministic seed)

> **참고:** 점근적 chi-sq 큐뮬런트는 소표본에서 큰 오차를 보이므로 사용하지 않습니다.
> 예: (3,3,3)에서 kappa_2: exact 2.72 vs asymptotic 4.0

### Saddlepoint Approximation -- Lugannani-Rice

Daniels (1954)의 새들포인트 밀도 근사:

$$f_{SP}(x) = \left(2\pi K_H^{(2)}(\hat{t})\right)^{-1/2} \exp(K_H(\hat{t}) - x\hat{t})$$

**Lugannani-Rice 꼬리확률 근사:**

$$\Pr(H \geq v) \approx 1 - \Phi(\hat{w}) + \phi(\hat{w})\left(\frac{1}{\hat{u}} - \frac{1}{\hat{w}}\right)$$

여기서:
- $\hat{w} = \text{sgn}(\hat{t})\sqrt{2(\hat{t}v - K_H(\hat{t}))}$
- $\hat{u} = \hat{t}\sqrt{K_H^{(2)}(\hat{t})}$

### CGF 근사 방법

**ER1 (Easton-Ronchetti 1st):**
$$K_H(t) \approx \sum_{i=1}^{4} \frac{\kappa_i t^i}{i!}$$

**ER2 (Easton-Ronchetti 2nd):**
$$K_H(t) \approx \kappa_1 t + \frac{\kappa_2}{2}t^2 + \log\left(1 + \frac{\kappa_3}{6}t^3 + \frac{3\kappa_4}{72}t^4 + \frac{\kappa_3^2}{72}t^6\right)$$

**Wang (damped, currently disabled and preserved in comments):**
$$K_H(t) \approx \kappa_1 t + \frac{\kappa_2}{2}t^2 + \left(\frac{\kappa_3}{6}t^3 + \frac{\kappa_4}{24}t^4\right)\eta_p(t)$$

여기서 $\eta_p(t) = \exp(-\kappa_2 p^2 t^2 / 2)$, $p$는 $K''_W(t;p) \geq 0$을 보장하는 최소 damping parameter.

**K-T (Kakizawa-Taniguchi):**
$$K_H(t) \approx \kappa_1 t + \frac{(1+\kappa_2)}{2}t^2 + \frac{\kappa_3}{6}t^3 + \frac{\kappa_4}{24}t^4$$

### 연속성 보정 (Continuity Correction)

이산 분포의 특성을 보정하기 위해 $v$를 $v - 1/2$로 대체하여 CC 변형 계산.

### Edgeworth Expansion

카이제곱 분포를 기저로 사용하며, 일반화된 Laguerre 다항식을 통해 skewness/kurtosis 보정:

$$F_{ED}(x) = G_{\nu}(x) - g_{\nu}(x)\left[\frac{\gamma_1}{6}L_3^{(\nu/2-1)}(x/2) + \frac{\gamma_2}{24}L_4^{(\nu/2-1)}(x/2) + \frac{\gamma_1^2}{72}L_6^{(\nu/2-1)}(x/2)\right]$$

### Gram-Charlier Type A

정규분포를 기저로 사용하며, Hermite 다항식 $H_3, H_4, H_6$을 통해 보정.

### PAG(d) -- Polynomially Adjusted Gamma

감마 분포에 $d$차 다항식을 곱한 밀도 근사:

$$f_{PAG}(x; d) = \psi(x) \sum_{i=0}^{d} \xi_i x^i$$

처음 $d+1$개 모멘트를 일치시키는 계수 $\xi_0, \ldots, \xi_d$를 moment matrix 역행렬로 결정. $d=4$ 및 $d=6$ 변형 제공.

## 모멘트와 큐뮬런트

귀무가설 하에서 H 통계량의 기본 모멘트:

- **평균**: $E(H) = k - 1$
- **분산**: Wallace (1959) exact formula:

$$\text{Var}(H) = 2(k-1) - \frac{2A_W}{5\,N(N+1)} - \frac{6}{5}\sum_{i=1}^{k}\frac{1}{n_i}$$

여기서 $A_W = 3k(k-2) + N(2k^2 - 6k + 1)$.

**큐뮬런트** (exact finite-sample):
- $\kappa_1 = E(H) = k - 1$
- $\kappa_2 = \text{Var}(H)$
- $\kappa_3 = \mu_3 - 3\mu_2\mu_1 + 2\mu_1^3$
- $\kappa_4 = \mu_4 - 4\mu_3\mu_1 - 3\mu_2^2 + 12\mu_2\mu_1^2 - 6\mu_1^4$

## 상세 사용법

### 개별 근사 클래스 사용

```python
from kw_approx import SaddlepointApproximation, ExactDistribution, MonteCarloSimulation

sample_sizes = [3, 3, 3]

# Saddlepoint with different CGF methods
sp_er1 = SaddlepointApproximation(sample_sizes, cgf_method='ER1')
sp_kt = SaddlepointApproximation(sample_sizes, cgf_method='KT')

print(f"ER1: {sp_er1.tail_probability_lr(4.62):.6f}")
print(f"KT:  {sp_kt.tail_probability_lr(4.62):.6f}")

# Exact (소표본)
exact = ExactDistribution(sample_sizes)
print(f"Exact: {exact.tail_probability(4.62):.6f}")

# Monte Carlo (대표본)
sim = MonteCarloSimulation(sample_sizes, n_simulations=10000, seed=42)
print(f"Simulation: {sim.tail_probability(4.62):.6f}")
```

### 임계값 계산

```python
from kw_approx import KWApproximator

approx = KWApproximator([5, 5, 5])

for alpha in [0.10, 0.05, 0.01]:
    cv_exact = approx.critical_value(alpha, method='exact')
    cv_chi2 = approx.critical_value(alpha, method='chi_square')
    cv_er1 = approx.critical_value(alpha, method='ER1')
    print(f"alpha={alpha}: Exact={cv_exact:.4f}, Chi2={cv_chi2:.4f}, ER1={cv_er1:.4f}")
```

### 실제 데이터에 적용

```python
import numpy as np
from scipy import stats
from kw_approx import KWApproximator

np.random.seed(42)
group1 = np.random.normal(10, 2, 5)
group2 = np.random.normal(12, 2, 5)
group3 = np.random.normal(11, 2, 5)

result = stats.kruskal(group1, group2, group3)
H = result.statistic

approx = KWApproximator([5, 5, 5])
print(f"H statistic: {H:.4f}")
print(f"SciPy p-value: {result.pvalue:.4f}")

for method in ['exact', 'chi_square', 'ER1', 'KT', 'edgeworth', 'gram_charlier', 'pam']:
    p = approx.tail_probability(H, method)
    print(f"{method}: {p:.4f}")
```

## 방법별 권장 사용 상황

```
+---------------+----------+---------------+---------------------------+
| 표본 크기 (N) | 그룹 (k) | 권장 방법     | 비고                      |
+---------------+----------+---------------+---------------------------+
| N <= 15       | k <= 3   | exact         | 정확 분포 계산 가능       |
| N <= 13       | k = 4    | exact         | 정확 분포 계산 가능       |
| 15 < N < 100  | any      | ER1, KT       | 새들포인트 + L-R          |
| N >= 100      | any      | chi_square    | 점근 근사 충분            |
| any           | any      | simulation    | MC 참조값 (느림)          |
+---------------+----------+---------------+---------------------------+
```

## 테스트

```bash
# 전체 테스트 실행 (36개)
python -m pytest tests/test_approximations.py -v

# 단일 테스트 클래스
python -m pytest tests/test_approximations.py::TestKWMoments -v

# 논문 테이블 재현
uv run python examples/reproduce_paper_tables.py
```

## 의존성

- Python >= 3.13
- NumPy >= 2.4.1
- SciPy >= 1.17.0

## 참고문헌

1. Kruskal, W. H. & Wallis, A. (1952). Use of ranks in one-criterion variance analysis. *JASA*, 47, 583-621.
2. Wallace, D. L. (1959). Simplified beta-approximations to the Kruskal-Wallis H test. *JASA*, 54, 225-230.
3. Daniels, H. E. (1954). Saddlepoint approximations in statistics. *Annals of Mathematical Statistics*, 25, 631-650.
4. Iman, R. L., Quade, D., & Alexander, D. A. (1975). Exact probability levels for the Kruskal-Wallis test.
5. Lugannani, R. & Rice, S. O. (1980). Saddlepoint approximation for the distribution of the sum of independent random variables. *Advances in Applied Probability*, 12, 475-490.
6. Easton, G. S. & Ronchetti, E. (1986). General saddlepoint approximations with applications to L statistics. *JASA*, 81, 420-430.
7. Wang, S. (1992). General saddlepoint approximations in the bootstrap. *Statistics & Probability Letters*, 13, 61-66.
8. Ha, H.-T. & Provost, S. B. (2007). A viable alternative to resorting to statistical tables. *Communications in Statistics*, 36, 1135-1151.
9. Hall, P. (1992). *The Bootstrap and Edgeworth Expansion*. Springer.

## 라이선스

MIT License
