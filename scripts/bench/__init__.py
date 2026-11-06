"""Shared evaluation harness for the AI-extension experiments (docs/DESIGN_SPEC.md §4).

registry: encoder features · tasks: targets built with the frozen label rules · folds: outer (fixed) and inner folds,
label-efficiency subsets · metrics: classification / regression scores · stats: DeLong, paired bootstrap, Wilson.
Run scripts from the repo root with `scripts` on the path (the bench_*.py scripts do this themselves).
"""
