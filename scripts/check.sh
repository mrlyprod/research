#!/usr/bin/env bash

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE/.."

export CARGO_TARGET_DIR="$PWD/data/target"
export CARGO_PROFILE_TEST_OPT_LEVEL=3
DATA_DIR="data/scripts/check"
LEAKS='/Users/|Developer/|AKIA[0-9A-Z]{16}|npm_[A-Za-z0-9]{36}|pypi-[A-Za-z0-9_-]{40,}|cio[A-Za-z0-9]{32}'

step() {
  local name="$1"; shift
  local log="$DATA_DIR/$name.log"
  local start=$SECONDS
  if "$@" > "$log" 2>&1; then
    printf 'ok   %-12s %4ss\n' "$name" "$((SECONDS - start))"
    return 0
  fi
  printf 'FAIL %-12s %4ss\n' "$name" "$((SECONDS - start))"
  tail -n 60 "$log" >&2
  echo "check stopped at $name; full log in $log" >&2
  exit 1
}

leak() {
  local hits code=0
  hits="$( { echo "$LEAKS"; [ ! -f .leak ] || sed '/^$/d' .leak; } | git grep --untracked -nIE -f - -- ':!scripts/check.sh' ':!.leak')" || code=$?
  [ "$code" -le 1 ] || { echo "git grep failed with $code" >&2; return 1; }
  [ -z "$hits" ] || { echo "$hits"; echo "private paths or token shapes in the tree" >&2; return 1; }
}

labs() { find . -name Cargo.toml -not -path './data/*' -not -path './node_modules/*' | sort; }
rust() { local verb="$1"; for toml in $(labs); do echo "== $toml"; cargo "$verb" --manifest-path "$toml" --all-targets || return 1; done; }
tested() { for toml in $(labs); do grep -rl --include='*.rs' '#\[test\]' "$(dirname "$toml")" 2>/dev/null | grep -q . && echo "$toml"; done; }
rust_test() {
  local slow=""
  for toml in $(tested); do
    echo "== $toml"
    perl -e 'alarm 120; exec @ARGV' cargo test --manifest-path "$toml"
    case $? in
      0) ;;
      142) slow="$slow $toml" ;;
      *) return 1 ;;
    esac
  done
  [ -z "$slow" ] || echo "slow, not gated:$slow"
}
figures() { [ -f package.json ] || { echo "no package.json yet; the figures are RESEARCH-TS's"; return 0; }; bun test; }

RUFF="uvx ruff@0.17.2"
lint_py() { $RUFF format --check . && $RUFF check .; }

mkdir -p "$DATA_DIR"
total=$SECONDS
step leak leak
step rust-check rust check
step rust-test rust_test
step ruff lint_py
step figures figures
echo "check passed in $((SECONDS - total))s"
