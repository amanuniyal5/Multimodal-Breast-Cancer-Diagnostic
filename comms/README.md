# How we communicate

The repo is our shared channel, for us and for our Claude Code sessions.

| Need | Use |
|---|---|
| Status update, result, handoff, quick question | Append an entry to [`LOG.md`](LOG.md) and push |
| What to work on next | `docs/PROJECT_PLAN.md` (no fixed owners; mark items `[~]` before starting) |
| A task not yet in the plan, with a done-condition | GitHub issue (template: *Task*), then add it to the plan |
| A question that blocks you | GitHub issue (template: *Question*) |
| Code or doc changes | `git pull --rebase`, then push to `main`; optional branch + PR for larger code changes |

## LOG.md entry format

Append at the bottom. Never edit or delete earlier entries.

```markdown
## 2026-10-13 14:30 IST · Uniyal (via Claude Code)

**Did:** what changed, with file paths and commit/PR links.
**Results:** numbers with n and CIs, or "none".
**Blocked on:** anything you need, or "nothing".
**Next open items:** ... (from `docs/PROJECT_PLAN.md`; nobody is assigned)
**Questions:** ... (or "none")
```

Write "(via Claude Code)" when your Claude wrote the entry, so the reader knows who to ask about it.
