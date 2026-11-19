Shorter answers with fewer, directly sourced details (on top of v1).

In the baseline and v1, faithfulness fell with answer length (100% under 200 words,
about 55% over 300), and 35-39 of 52 answers ignored the 250-word guideline. Failures
were mostly extra details: hand-counted rows, inferred background, side remarks.

- `chat.prompt`: under 200 words; the direct answer, about five supporting numbers or
  ids at most, and the one caveat that matters most; offer a follow-up instead of
  listing everything; every number is a tool's number or simple arithmetic on it, never
  a hand count of rows; no side remarks about things not queried; tables of about six
  rows at most. Requested lists stay complete when short.
- Keeps v1's partial-result notes, ordering descriptions and feature_type_counts.

Same model (Opus 5.5, medium effort), questions, judge and rubric as the baseline.
