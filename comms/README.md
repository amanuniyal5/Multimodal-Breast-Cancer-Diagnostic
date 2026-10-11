# How we communicate

The repo is our shared channel, for us and for our Claude Code sessions.

| Need | Use |
|---|---|
| Status update, result, handoff, quick question | Append an entry to [`LOG.md`](LOG.md) and push |
| A task with an owner and a done-condition | GitHub issue (template: *Task*) |
| A question that blocks you | GitHub issue (template: *Question*), mention the other person |
| Code or doc changes | Branch `uniyal/<topic>` or `srivastava/<topic>` → pull request into `main` |

## LOG.md entry format

Append at the bottom. Never edit or delete earlier entries.

```markdown
## 2026-10-13 14:30 IST · Uniyal (via Claude Code)

**Did:** what changed, with file paths and commit/PR links.
**Results:** numbers with n and CIs, or "none".
**Blocked on:** anything you need, or "nothing".
**Next for me:** ...
**Ask for Srivastava:** ... (or "none")
```

Write "(via Claude Code)" when your Claude wrote the entry, so the reader knows who to ask about it.
