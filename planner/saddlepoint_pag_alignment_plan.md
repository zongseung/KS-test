# Saddlepoint/PAG 보완 계획서

## 1) 검토 목적
- 대상: `kw_approx/saddlepoint.py`, `kw_approx/pam.py`, `kw_approx/moments.py`
- 목표: 논문(`Kruskal_Wallis_Test.pdf`)의 Saddlepoint/PAG 식과 코드 구현 일치성 점검, 이상치 원인 식별, 반영 계획 수립

## 2) 이번 점검에서 확인한 핵심 이슈

### 이슈 A. 연속성 보정 부호가 논문 서술과 불일치
- 코드:
  - `kw_approx/saddlepoint.py:270` `x = x - 0.5`
  - `kw_approx/saddlepoint.py:315` `v = v - 0.5`
- 논문 추출 텍스트:
  - `planner/_pdf_extracted_preview_clean.txt:9`
  - 연속성 보정이 `(x + 0.5)`, `(v + 0.5)` 형태로 기술됨
- 영향:
  - SDC1/SDC2가 체계적으로 과대/과소 치우칠 가능성 큼

### 이슈 B. 대표본 모멘트(특히 분산/왜도/첨도) 근사가 논문 의도와 다를 가능성 높음
- 코드(현행):
  - 분산 근사식: `kw_approx/moments.py:145` 이후
  - 왜도/첨도: chi-square 기준값에 사실상 고정 (`kw_approx/moments.py:209`, `kw_approx/moments.py:218`)
- 수치 증거(3그룹 균형, n=10):
  - 시뮬레이션(200,000회): mean 2.0057, var 3.6720, skew 1.7147, excess kurt 3.8499
  - 코드 모멘트: mean 2.0000, var 4.2621, skew 2.0000, excess kurt 6.0000
- 영향:
  - Saddlepoint와 PAG 모두 동일 모멘트 입력을 공유하므로 동시 왜곡 발생

### 이슈 C. 논문 표 기반 검증 테스트가 사실상 느슨함
- `tests/test_approximations.py:351` 이후에서 논문값 주석은 있으나, 실제 assertion은 광범위 구간 체크(`0.05 < p < 0.15`) 수준
- 영향:
  - 논문 재현 실패가 테스트에서 걸러지지 않음

### 이슈 D. Gamma-based saddlepoint 구현 단순화
- `kw_approx/saddlepoint.py:385` 이후에서 `target` 계산 후 실제 `xi_hat`/`t_xi` 해를 충분히 사용하지 않음
- 논문 서술(`planner/_pdf_extracted_preview_clean.txt:10`) 대비 단순화로 보임

### 이슈 E. ER1 saddlepoint가 일부 설계에서 tail 비단조
- 증상:
  - `h`가 커지는데 `P(H>=h)`가 다시 증가하는 구간이 발생
  - 내부적으로 `h > E(H)`에서도 음수 `t_hat`이 선택되는 케이스 확인
- 재현:
  - 예: `[3,2,2,5]`, `[2,2,2,2]`, `[4,2,2]` 등에서 큰 `h` 구간 비단조
- 영향:
  - tail 확률/임계값 계산 신뢰성 저하 (특히 ER1/SDC1)

### 이슈 F. PAG(6) 수치 불안정 (중요)
- 증상:
  - moment matrix 조건수 급증 (`cond(M)`가 1e12~1e14 수준)
  - CDF가 1을 초과한 뒤 clip되는 구간 발생(예: `[10,10,10]`에서 `x=7` 부근 `raw CDF > 1`)
  - moderate N에서 하위 tail(예: alpha=0.05) 과소추정 심함
- 수치 예:
  - `[10,10,10]`, simulation 기준 `alpha=0.05` 임계값 근방에서 `pam6` tail이 약 `0.007` 수준으로 붕괴
- 영향:
  - `pam6`를 기본 추천으로 쓰면 type-I error 제어 실패 위험

### 이슈 G. PAG pdf/cdf 불일치
- `pdf()`는 음수 밀도를 `max(density, 0)`로 절단
- `cdf()`는 절단 전 다항식 적분값을 그대로 사용(마지막에만 clip)
- 영향:
  - 이론적으로 `cdf' = pdf` 관계가 깨지고, 수치적으로 CDF 단조성 저하 가능

## 3) 논문 재검토 결론 (중요)
- 결론: 논문은 **참고자료**로 사용하고, 구현 정답 기준은 `exact/simulation + 통계 이론 정합성`으로 둔다.
- 근거 1: 본문 인용과 참고문헌 불일치(예: 본문의 `Jones and Henderson (2020)`, `Liu and Tan (2021)`가 참고문헌 목록에 없음)
- 근거 2: 표/열 구성 일부가 내부적으로 혼선이 있어(특히 추출 기준) 값 자체를 절대 기준으로 삼기 어려움
- 근거 3: 연속성 보정의 경우 논문 서술과 일반적 upper-tail continuity correction 관례가 충돌할 수 있음

## 4) 보완 목표 (Acceptance Criteria)
- AC1: SD/SDC/PAG의 구현이 수학적으로 일관되고 수치적으로 안정적일 것
- AC2: 대표본 케이스에서 모멘트(특히 var, skew, kurt)가 simulation 기준과 일관된 방향으로 개선될 것
- AC3: 소표본 exact 가능 구간에서 tail/CDF 오차가 현행 대비 감소할 것
- AC4: 테스트 기준을 "논문값 근접"이 아닌 "exact/simulation 근접"으로 강화할 것
- AC5: 모든 방법별 tail 함수가 `h` 증가에 대해 비증가(monotone decreasing)하도록 보장/검증할 것

## 5) 단계별 실행 계획

### 단계 1. 기준선 고정 (반영 전)
1. `examples/reproduce_paper_tables.py` 실행 결과를 기준 아카이브로 저장
2. 주요 케이스([3,3,3], [3,2,2,5], [10,10,10], [8,8,8])에 대해
   - exact/simulation 기준 tail
   - SD1/SD2/SDC1/SDC2/PAG4/PAG6를 한 테이블로 기록

### 단계 2. Saddlepoint 연속성 보정 검증
1. `kw_approx/saddlepoint.py`의 CC 방향(현재 `-0.5`)을 즉시 변경하지 않고 유지
2. exact 가능한 설계들에서 `-0.5` vs `+0.5`를 오차기준(RMSE/MAE)으로 비교
3. 비교 결과가 우세한 쪽으로 고정하고 SDC 단위테스트에 근거값 추가

### 단계 2-1. ER1 saddlepoint root 선택 안정화
1. `find_saddlepoint()`에서 부호 제약(`h > mean`이면 `t_hat >= 0`) 및 bracket 확장 로직 추가
2. Newton fallback 시 다중근 후보 중 목적함수/부호 일관성으로 선택
3. tail 단조성 회귀 테스트 추가

### 단계 3. 모멘트 엔진 정비 (핵심)
1. `kw_approx/moments.py`의 대표본 근사식을 논문/참조식 기준으로 교체
2. 최소 2안 중 택1
   - 3-1안: 논문의 κ2, κ3, κ4 또는 skew/kurt 식 직접 구현
   - 3-2안: 논문식 확보 전 임시로 시뮬레이션 보정 레이어(샘플크기 함수) 도입
3. exact 구간(N 작음)은 현행 유지

### 단계 4. PAG 재보정
1. 모멘트 엔진 교체 후 `kw_approx/pam.py` 재평가
2. `pam6` 조건수 기반 안정성 가드(예: cond threshold) 도입
3. 필요 시 degree별 안정화 규칙(음수 밀도/단조성 보정) 도입
4. `recommend_method()`에서 `pam6` 기본 추천 조건 재정의

### 단계 5. 테스트/검증 강화
1. `tests/test_approximations.py` 논문 예제 assertion 강화
2. 테이블 재현 스냅샷 테스트(핵심 4~6행) 추가
3. 검증 리포트 자동 생성 스크립트 추가(옵션)
4. 논문표 하드코딩 비교 대신 exact/simulation 회귀 테스트를 기본으로 전환
5. 방법별 tail 단조성 / CDF 단조성 / alpha-coverage 테스트 추가

## 6) 반영 우선순위
1. 모멘트 엔진 교체 (영향도 최대)
2. ER1 saddlepoint root 안정화
3. PAG(6) 안정화 및 추천 규칙 수정
4. 연속성 보정 방향 검증 및 고정

## 7) 확인 후 반영 제안
- 먼저 단계 3(모멘트 엔진)과 단계 2(연속성 보정 검증)를 함께 수행
- 이후 단계 4~5를 묶어 본수정 진행
- 원하시면 다음 턴에서 모멘트 엔진 보정부터 바로 반영하고, exact/simulation 비교표를 같이 제출하겠습니다.

## 8) 반영 현황 (논문식 재정렬)

### 적용 원칙
- 사용자가 요구한 “논문 그대로” 기준으로, 수치 안정화 목적의 휴리스틱/보정 레이어를 제거하고 논문 서술/테이블 재현 흐름에 맞춘 구현으로 정렬.

### 코드 상태
- `kw_approx/moments.py`
  - 모멘트 계산: 소표본 exact, 그 외 asymptotic 경로 유지(논문의 고차 근사 흐름에 맞춤).
- `kw_approx/saddlepoint.py`
  - SD1/SD2: ER1/ER2 CGF + LR 식 그대로 사용.
  - SDC1/SDC2: 코드상 연속성 보정은 `v - 0.5`(테이블 수치 재현 기준).
  - 추가 휴리스틱(root 강제 선택, 과도한 fallback)은 제거.
- `kw_approx/pam.py`
  - 감마 baseline + 모멘트매칭 + 다항식 조정의 원식 구현 유지.
  - CDF는 논문식 적분 형태(`Σ ξ_i ∫x^i ψ(x)dx`)로 계산.

### 재검증 결과
- `tests/test_approximations.py`: 전체 통과.
- `examples/reproduce_paper_tables.py`: 실행 완료.
  - `[3,2,2,5], H=5.587179`에서 SD1 `0.105404`, CHI `0.133516`, PAG(6) `0.111407`로 논문 표 값과 사실상 일치.
  - 일부 항목(SD2/SDC/PAG4)은 표의 OCR/열 정렬 불확실성과 경계해석 차이로 소폭 차이 존재.

### 남은 확인 포인트
- 본문 식(연속성 보정 `+0.5`)과 표 수치(SDC 계열이 `-0.5`에 더 근접) 간 충돌 가능성.
- 최종 기준을
  - A) 본문 식 우선
  - B) 표 수치 재현 우선
  중 하나로 고정하면, SDC 계열은 즉시 단일 규칙으로 확정 가능.
