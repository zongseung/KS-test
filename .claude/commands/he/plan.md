---
name: plan
description: "논문 비교분석 기획서 - 두 논문 버전 간 차이점 추출 및 변경 유형 분류"
category: harness-engineering
phase: planning
---

# /he:plan - 기획서 (Planner)

> 하네스 엔지니어링 1단계: 두 논문 버전을 비교하여 차이점을 추출하고 변경 유형을 분류합니다.

## Triggers
- 논문의 수정 전/후 버전이 `comparison/` 디렉토리에 존재할 때
- 새로운 논문 리비전이 들어왔을 때
- 기존 구현과 논문 간 정합성 점검이 필요할 때

## Context Trigger Pattern
```
/he:plan [--focus formula|table|method|all] [--depth shallow|deep]
```

### Options
- `--focus`: 특정 영역에 집중 (수식/테이블/방법론/전체)
- `--depth shallow`: 구조적 차이만 빠르게 확인
- `--depth deep`: 수식 레벨까지 상세 비교

## 입력 파일
- **Paper A**: `comparison/Kruskal_Wallis_Test_revisionn.pdf` (수정 전)
- **Paper B**: `comparison/Kruskal_Wallis_Test_final_revision.pdf` (수정 후)
- **기존 분석**: `planner/saddlepoint_pag_alignment_plan.md`, `revision/revision1.md`

## Behavioral Flow

### Phase 1: 텍스트 추출 및 정렬
1. **PDF 텍스트 추출**: 두 논문 모두 `pdftotext`로 텍스트 추출
2. **섹션 정렬**: 동일 섹션끼리 매핑 (Introduction, Methods, Tables 등)
3. **diff 생성**: 섹션별 텍스트 차이 추출

### Phase 2: 차이점 분류
각 차이를 아래 유형으로 분류:

| 유형 | 설명 | 코드 영향도 |
|------|------|-------------|
| **FORMULA** | 수식 변경/추가/삭제 | 높음 |
| **METHOD** | 알고리즘/방법론 변경 | 높음 |
| **TABLE** | 테이블 구조/값 변경 | 중간 |
| **REFERENCE** | 참고문헌 추가/삭제 | 낮음 |
| **TEXT** | 서술/표현 변경 | 낮음 |
| **STRUCTURE** | 섹션 재구성/순서 변경 | 낮음 |
| **AUTHOR** | 저자 순서/소속 변경 | 없음 |

### Phase 3: 영향도 분석
1. **코드 매핑**: 각 FORMULA/METHOD 변경이 영향을 주는 코드 파일 식별
   - `kw_approx/moments.py` — 모멘트/큐뮬런트 수식
   - `kw_approx/saddlepoint.py` — 새들포인트 근사
   - `kw_approx/pam.py` — 다항식 조정 감마 근사
   - `kw_approx/edgeworth.py` — 에지워스/그램-샬리에 전개
   - `examples/reproduce_paper_tables.py` — 테이블 재현
   - `tests/test_approximations.py` — 검증 테스트
2. **기존 이슈 연결**: `planner/saddlepoint_pag_alignment_plan.md`의 이슈 A~G과 매핑
3. **우선순위 결정**: 코드 영향도 + 기존 이슈 관련성으로 반영 우선순위 결정

## Output Format

```markdown
## 논문 비교 기획서: Paper A → Paper B

### 1. 변경 요약
| # | 섹션 | 유형 | 변경 내용 요약 | 코드 영향 |
|---|------|------|---------------|-----------|
| 1 | 2.3  | FORMULA | κ3 수식 수정 | moments.py |
| 2 | 3.1  | METHOD  | SDC 연속성 보정 방향 | saddlepoint.py |
| ...                                              |

### 2. 핵심 변경 상세
#### 변경 #1: [제목]
- **Paper A**: (기존 내용)
- **Paper B**: (수정 내용)
- **차이**: (구체적 변경 설명)
- **영향 코드**: [파일:라인]
- **기존 이슈 연결**: 이슈 B (모멘트 근사)

### 3. 새로 추가된 내용 (Paper B only)
- 항목 1: ...
- 항목 2: ...

### 4. 삭제/축소된 내용 (Paper A only)
- 항목 1: ...

### 5. 검증 기준 (판별자 입력)
- [ ] 모든 FORMULA 변경이 식별되었는가
- [ ] 코드 영향 매핑이 정확한가
- [ ] 기존 planner 이슈와의 연결이 완전한가
```

## Handoff

기획 완료 후 생성자(`/he:generate`)에게 넘겨 반영 계획서를 작성합니다.

**다음 단계**: `/he:generate --from-plan`

## Boundaries

**Will:**
- 두 논문의 차이점을 체계적으로 추출 및 분류
- 코드에 미치는 영향도를 파일/라인 수준으로 매핑
- 기존 분석 문서(`planner/`, `revision/`)와 연결

**Will Not:**
- 코드를 직접 수정하지 않음 (생성자의 역할)
- 변경의 옳고 그름을 판정하지 않음 (판별자의 역할)
- 논문의 학술적 품질을 평가하지 않음
