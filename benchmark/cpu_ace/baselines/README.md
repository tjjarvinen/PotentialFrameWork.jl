Historical baseline: `ea076d3c18d5ae66776f1fa785e729f2eadfba54` (clean tracked tree).

Before editing, `git archive HEAD -o /tmp/pfw-stage1-baseline.tar` produced SHA-256
`9f4e36bbd310404dd1c8d6d9e365d1a1e24d546884e2fab88608db2d3dd519fd`.
Extracted at `/tmp/pfw-stage1-baseline`. Recreate from the pinned commit using
`git archive ea076d3c18d5ae66776f1fa785e729f2eadfba54 | tar -x -C DEST`.
The supplied untracked instructions are preserved in the working checkout.
The commit is the durable immutable source snapshot; no duplicate repository is stored here.
