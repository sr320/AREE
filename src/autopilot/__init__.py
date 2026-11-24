"""Klone autopilot: watch AREE study runs on UW Hyak (Klone) and report their state.

Phase 1 is report-only: it reads Slurm, Nextflow and storage state, writes a
ledger, edits a status issue on GitHub and opens alert issues. It never submits
or cancels a job. See docs/klone_autopilot_plan.md.
"""
