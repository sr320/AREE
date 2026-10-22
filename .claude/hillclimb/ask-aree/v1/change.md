Tell the chat when it sees only part of a result, and to stay within what it retrieved.

The baseline failed only on faithfulness, mostly by generalizing from the first rows of
a shortened result to the whole set, plus a few answers describing studies never looked up.

- `chat.tools`: a result that returns fewer rows than its total now carries
  `partial_result`, saying how many rows were left out and that claims about them must
  rest on the totals. Tool descriptions now state how rows are ordered.
- `query.find_replicated_features`: reports `feature_type_counts`, so a count is never
  described as covering feature types it doesn't.
- `chat.prompt`: rules to base whole-set claims on totals, to describe a study only from
  its record, and never to claim a check that wasn't run.

Same model (Opus 5.5, medium effort), questions, judge and rubric as the baseline.
