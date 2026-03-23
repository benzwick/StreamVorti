#!/bin/bash
# Run all demo examples locally with ParaView output.
#
# Usage:
#   ./scripts/run-demos.sh              # run from repo root (uses build/)
#   ./scripts/run-demos.sh /path/to/build
#
# Results are written to results/ under the repo root.
# Open the .pvd files in ParaView to view the output.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="${1:-$REPO_ROOT/build}"
BINARY="$BUILD_DIR/StreamVorti"

if [ ! -x "$BINARY" ]; then
    echo "Error: $BINARY not found or not executable."
    echo "Build first: cd build && cmake .. && make -j6"
    exit 1
fi

# Each entry: "demo_name|expected_to_fail|failure_reason"
DEMOS=(
    "cavity||"
    "cavity-fdm||"
    "cavity-quick||"
    "cavity-re1000||"
    "cavity-regularized||"
    "cavity-oscillating||"
    "cavity-tall||"
    "double-lid-cavity||"
    "channel||"
    "poiseuille||"
    "backward-step|yes|compound boundary predicates (and/or) not yet supported (issue #30)"
    "cylinder|yes|SDL::DIFFERENCE (CSG subtraction) not implemented (issue #25)"
    "von-karman|yes|SDL::DIFFERENCE (CSG subtraction) not implemented (issue #25)"
)

PASSED=0
FAILED=0
XFAILED=0
XPASSED=0
FAILURES=""

mkdir -p "$REPO_ROOT/results"

for entry in "${DEMOS[@]}"; do
    IFS='|' read -r demo expect_fail reason <<< "$entry"
    DEMO_FILE="$REPO_ROOT/demo/${demo}.lisp"

    if [ ! -f "$DEMO_FILE" ]; then
        echo "SKIP: $DEMO_FILE not found"
        continue
    fi

    echo ""
    echo "========================================"
    echo "Running: $demo"
    if [ -n "$expect_fail" ]; then
        echo "  (expected to fail: $reason)"
    fi
    echo "========================================"

    LOG="$REPO_ROOT/results/${demo}.log"
    if (cd "$REPO_ROOT" && "$BINARY" -f "$DEMO_FILE" -lp "$REPO_ROOT/lisp" -pv) > "$LOG" 2>&1; then
        if [ -n "$expect_fail" ]; then
            echo "XPASS: $demo (expected to fail but passed!)"
            XPASSED=$((XPASSED + 1))
        else
            echo "PASS: $demo"
            PASSED=$((PASSED + 1))
        fi
    else
        if [ -n "$expect_fail" ]; then
            echo "XFAIL: $demo ($reason)"
            XFAILED=$((XFAILED + 1))
        else
            echo "FAIL: $demo (see $LOG)"
            FAILED=$((FAILED + 1))
            FAILURES="$FAILURES $demo"
        fi
    fi
done

echo ""
echo "========================================"
echo "Summary: $PASSED passed, $FAILED failed, $XFAILED expected failures, $XPASSED unexpected passes"
if [ -n "$FAILURES" ]; then
    echo "Unexpected failures:$FAILURES"
fi
echo "========================================"
echo ""
echo "Results in: $REPO_ROOT/results/"
echo "Open .pvd files in ParaView to view output."

exit $FAILED
