# 재현 가능한 그림

`genomecurator.treefigure.render()`는 BEAST 분석의 요약 트리를 한 번의 함수 호출로 논문용 그림으로 그리고, `genomecurator.repro.manifest()`는 그림마다 매니페스트를 옆에 써서
같은 입력으로 그림을 다시 만들고 확인할 수 있게 합니다.

```python
from genomecurator import beastprep, nextstrain, treefigure
run = nextstrain.get_run(RUN_ID)                       # BEAST 분석의 클럭 실행
tree = beastprep.annotated_tree(run, "my-analysis")    # 끝마디에 나라·리니지·날짜가 붙은 TreeAnnotator 요약 트리
treefigure.render(tree, "figures/tree", focus={"tipA", "tipB"}, labels={"tipA": "alias A"}, run=run,
                  intervals=[{"label": "MRCA", "lo": 2024.30, "hi": 2024.50, "on": "clade"}], title="...", collapse_min=3)
```

- **그리는 것.** 달력 축 위의 트리. 관심 클레이드가 마지막에 오는 사다리 순서, `collapse_min`개 이상 끝마디를 가진 곁가지는 나라 구성을 적은 삼각형으로 접기, 끝마디마다 리니지 칩과 수집월,
  사후확률 `support` 이상 노드, 관심 클레이드 음영, `intervals`로 준 신용구간(클레이드 MRCA 또는 줄기)을 그립니다.
- **비공개 유전체.** 비공개 유전체(실행 멤버 중 `private:`로 시작하는 것, 또는 `private_ids`의 이름)는 `labels`에 별칭이 있어야 하고, 없으면 그림을 만들지 않습니다.
- **매니페스트.** `<그림>.manifest.json`에 트리와 모든 입력 파일의 SHA-256과 크기, 매개변수, 프로그램 버전과 git 커밋(그 뒤 수정 여부), Python·라이브러리·BEAST·TreeAnnotator 버전,
  쓴 파일마다 SHA-256을 적습니다. `repro.verify(매니페스트)`가 달라진 파일을 알려 줍니다.
- **같은 입력, 같은 파일.** SVG는 같은 입력이면 바이트까지 같습니다(시각 정보 없음, 해시 솔트 고정). PNG와 PDF에는 생성 날짜가 없습니다.
- **설치.** matplotlib은 선택입니다: `pip install genomecurator[figures]`.

매니페스트가 보장하지 않는 것: BEAST 체인은 확률적입니다. 같은 시드, 같은 BEAST·BEAGLE 빌드와 스레드 수면 체인이 되풀이되어야 하지만 CPU나 스레드 수가 다르면 달라질 수 있습니다.
요약(중앙값, HPD)을 비교하고, 다른 시드를 적어도 하나 더 돌리세요.
