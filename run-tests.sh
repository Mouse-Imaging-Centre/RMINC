#!/bin/bash
set -euo pipefail

MODE="${1:-all}"
FILTER="${2:-}"

if [[ "${FILTER}" =~ [^A-Za-z0-9_-] ]]; then
    echo "ERROR: filter contains invalid characters (allowed: alphanumeric, underscore, dash)" >&2
    exit 1
fi

echo "=== Building RMINC from mounted source ==="
cd /RMINC
rm -f RMINC_*.tar.gz
R CMD build .
R CMD INSTALL RMINC_*.tar.gz

run_check() {
    echo "=== Running R CMD check --as-cran ==="
    cd /RMINC
    R CMD check --as-cran --no-install --no-manual RMINC_*.tar.gz
    if [ -f RMINC.Rcheck/00check.log ]; then
        cat RMINC.Rcheck/00check.log
    fi
    local status
    status=$(grep '^Status:' RMINC.Rcheck/00check.log 2>/dev/null || echo "unknown")
    echo "$status"
    if echo "$status" | grep -q "ERROR"; then
        echo "=== R CMD check had ERRORS ==="
        exit 1
    fi
    echo "=== R CMD check passed ==="
}

run_tests() {
    echo "=== Running tests (${FILTER:-all}) ==="
    cd /RMINC
    if [ -n "${FILTER}" ]; then
        Rscript -e "library(RMINC); devtools::test('.', filter = '${FILTER}')"
    else
        Rscript -e "library(RMINC); devtools::test('.')"
    fi
}

case "$MODE" in
    check)
        run_check
        ;;
    test)
        run_tests
        ;;
    all)
        run_check
        run_tests
        ;;
    shell)
        exec bash
        ;;
    *)
        echo "Usage: run-tests.sh [all|check|test|shell] [filter]"
        echo ""
        echo "  all     Run R CMD check + tests (default)"
        echo "  check   Run R CMD check --as-cran only"
        echo "  test    Run tests only"
        echo "  shell   Drop into an interactive bash shell"
        echo ""
        echo "  test <filter>  Run only tests matching <filter> (e.g. mincLm)"
        echo ""
        echo "Mount source with -v \$(pwd):/RMINC to test local changes."
        exit 1
        ;;
esac
