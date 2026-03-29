---
name: he
description: "하네스 엔지니어링 디스패처 - 논문 비교분석 파이프라인 (기획서/생성자/판별자)"
category: harness-engineering
phase: orchestration
---

# /he - 하네스 엔지니어링 (Harness Engineering)

> 논문 버전 간 비교분석을 위한 기획서-생성자-판별자 GAN 패턴 프레임워크

## Overview

동일 논문의 수정 전/후 버전을 체계적으로 비교분석하여, 코드 구현에 반영해야 할 변경사항을 도출합니다:

```
[기획서/Planner] → [생성자/Generator] ⇄ [판별자/Discriminator]
  논문 diff 분석      반영 계획서 작성       누락/오류 검증
```

## 대상 논문
- **Paper A (수정 전)**: `comparison/Kruskal_Wallis_Test_revisionn.pdf`
- **Paper B (수정 후)**: `comparison/Kruskal_Wallis_Test_final_revision.pdf`
- **기존 구현 기준**: `Kruskal_Wallis_Test (2).pdf`

## Commands

| 커맨드 | 역할 | 설명 |
|--------|------|------|
| `/he:plan` | 기획서 | 두 논문 간 차이점 추출 및 분석 |
| `/he:generate` | 생성자 | 차이점 기반 반영 계획서/코드 수정안 작성 |
| `/he:judge` | 판별자 | 누락 항목 검증, 구현 정합성 판정 |
| `/he` | 디스패처 | 전체 파이프라인 오케스트레이션 |

## Full Pipeline

```
/he [--auto] [--max-iter 3]
```

### Pipeline Flow

```
1. /he:plan    → 논문 diff 분석
                  - Paper A vs Paper B 텍스트/수식/테이블 차이 추출
                  - 변경 유형 분류 (수식, 방법론, 표현, 구조)
                    ↓
2. /he:generate → 반영 계획서 작성
                  - 코드에 반영해야 할 변경사항 목록
                  - 빠진 구현/새로 추가해야 할 구현 식별
                  - 기존 planner 문서 업데이트 초안
                    ↓
3. /he:judge   → 검증 및 판정
                  - 생성자가 놓친 차이점 없는지 교차 검증
                  - 구현 코드와의 정합성 체크
                  - PASS/FAIL 판정 후 피드백
                    ↓
              ┌─ PASS → 최종 반영 계획서 확정
              ├─ CONDITIONAL_PASS → 사용자 확인 후 반영
              └─ FAIL → /he:generate --iterate (재분석)
```

## Quick Start

### 단계별 실행 (권장)
```bash
# 1단계: 두 논문 차이점 분석
/he:plan

# 2단계: 반영 계획서 생성
/he:generate --from-plan

# 3단계: 누락 검증
/he:judge

# 피드백이 있으면 재생성
/he:generate --iterate
```

## Design Principles

1. **수식 우선**: 텍스트 변경보다 수식/알고리즘 변경을 우선 분석
2. **구현 연동**: 단순 논문 비교가 아닌, 코드 반영 관점에서 분석
3. **누락 방지**: 판별자가 생성자의 분석을 교차 검증하여 빠진 항목 탐지
4. **기존 문서 연계**: `planner/`, `revision/`, `result/`의 기존 분석과 연결

## Version
Harness Engineering v1.0.0 — Paper Comparison Mode
