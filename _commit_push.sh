#!/usr/bin/env bash
set -e
cd ~/Academia/Projects/Internship_USTC

git add -A
echo "=== staged summary ==="
git status --short | sed -n '1,8p'
echo "  ... ($(git status --short | wc -l) total changed paths)"

git commit -F - <<'MSG'
Reorganize: stelloft as editable package under Code/, untrack Source_Codes

- Move first-party stelloft package from Code/Source_Codes/ to Code/stelloft/
- Add pyproject.toml and install editable so `import stelloft` works from any
  directory or notebook
- Untrack + gitignore Code/Source_Codes/ (local reference copies of third-party
  packages installed separately from ~/Academia/Software)
- Ignore build/output artifacts (*.torosurf, *.stl, *.egg-info, Zone.Identifier)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
MSG

echo "=== commit done ==="
git log --oneline -1
echo "=== push ==="
git push origin main 2>&1 | tail -8
