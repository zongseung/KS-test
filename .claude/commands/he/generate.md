---
name: generate
description: "논문 비교분석 생성자 - 차이점 기반 반영 계획서 작성, 코드 수정안 도출, 누락 구현 식별"
category: harness-engineering
phase: generation
---

# /he:generate - 생성자 (Generator)

> 하네스 엔지니어링 2단계: 기획서의 차이점 분석을 기반으로 코드 반영 계획서를 작성합니다.

## Triggers
- `/he:plan` 이후 반영 계획 생성 요청
- 논문 차이점 목록이 이미 존재하는 상태에서 실행
- 판별자 피드백 기반 재생성 (`--iterate`)

## Context Trigger Pattern
```
/he:generate [--from-plan] [--iterate] [--output plan|code|both]
```

### Options
- `--from-plan`: 기획서의 diff 결과를 기반으로 생성
- `--iterate`: 판별자 피드백 반영하여 재생성 (GAN 루프)
- `--output plan`: 반영 계획서만 생성
- `--output code`: 코드 수정안만 생성
- `--output both`: 계획서 + 코드 수정안 모두 생성

## Behavioral Flow

### Phase 1: 기획서 로드
1. **차이점 목록 확인**: `/he:plan` 결과의 변경 요약 테이블 로드
2. **우선순위 확인**: FORMULA > METHOD > TABLE > TEXT 순으로 처리
3. **기존 문서 확인**: `planner/`, `revision/` 기존 계획서와 중복/충돌 확인

### Phase 2: 반영 항목 도출
각 차이점에 대해 3가지 관점으로 분석:

#### 2-1. 구현 반영 필요 항목
- Paper B에서 수식이 바뀐 경우 → 해당 코드 수정 필요
- 새로운 방법론이 추가된 경우 → 새 구현 필요
- 테이블 구조가 변경된 경우 → 재현 스크립트 수정 필요

#### 2-2. 빠진 구현 (Missing Implementation)
- Paper B에 있지만 현재 코드에 구현되지 않은 수식/알고리즘
- 기존 `planner/` 이슈에서 미해결로 남은 항목 중 Paper B가 해결한 것
- Paper B에서 새로 등장한 정의/조건

#### 2-3. 기획서 업데이트 항목
- `planner/saddlepoint_pag_alignment_plan.md` 수정이 필요한 부분
- `revision/revision1.md` 업데이트가 필요한 부분
- 새로운 이슈/검증 항목 추가

### Phase 3: 반영 계획서 작성
기획서 diff를 실행 가능한 작업 목록으로 변환:

1. **코드 수정 사항**: 파일별/함수별 구체적 수정 내용
2. **새 구현 사항**: 추가해야 할 함수/로직
3. **테스트 수정 사항**: 변경된 수식에 맞는 테스트 업데이트
4. **문서 수정 사항**: planner/revision 문서 업데이트

### GAN 루프 (--iterate)
판별자 피드백을 받은 경우:
1. 누락으로 지적된 차이점 재분석
2. 코드 매핑 오류 수정
3. 반영 계획서 보완 후 재제출

```
[기획서 diff] → [생성자] → [판별자] → 누락 피드백 → [생성자 재생성] → ...
```

## Output Format

```markdown
## 반영 계획서: Paper A → Paper B

### 1. 코드 수정 필요 항목
| # | 변경 근거 | 대상 파일 | 현재 코드 | 수정 방향 | 우선순위 |
|---|-----------|-----------|-----------|-----------|----------|
| 1 | 기획서 #3: κ3 수식 | moments.py:L145 | 기존 근사 | Paper B 수식 반영 | P0 |
| 2 | 기획서 #5: SDC 보정 | saddlepoint.py:L270 | v - 0.5 | Paper B 확인 후 결정 | P1 |

### 2. 빠진 구현 (신규 추가 필요)
| # | Paper B 내용 | 구현 위치(제안) | 설명 |
|---|-------------|----------------|------|
| 1 | 새 정리/공식 | kw_approx/xxx.py | ... |

### 3. 기존 planner 문서 업데이트
| 문서 | 수정 내용 |
|------|-----------|
| saddlepoint_pag_alignment_plan.md | 이슈 B: Paper B 반영으로 상태 변경 |
| revision1.md | 섹션 8: 반영 현황 업데이트 |

### 4. 테스트 수정 사항
- test_approximations.py: 변경된 수식 기준값 업데이트
- reproduce_paper_tables.py: Paper B 테이블 구조 반영

### 5. 자체 검증
- [x] 모든 FORMULA 변경에 대한 수정안 작성됨
- [x] 코드 파일 경로 및 라인 번호 확인됨
- [ ] 판별자 검증 대기
```

## Handoff

생성 완료 후 판별자(`/he:judge`)에게 검증을 요청합니다.

**다음 단계**: `/he:judge`

## Boundaries

**Will:**
- 논문 차이점을 구현 가능한 작업 목록으로 변환
- 빠진 구현(Missing Implementation) 식별
- 기존 planner 문서와의 연결 및 업데이트 초안 작성
- 판별자 피드백 기반 반복 개선

**Will Not:**
- 코드를 직접 수정하지 않음 (별도 구현 단계에서 수행)
- 논문 변경의 수학적 정합성을 검증하지 않음 (판별자의 역할)
- 기획서 범위를 벗어나는 분석 수행
