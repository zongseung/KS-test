#!/bin/bash
#
# Validation Script for Kruskal-Wallis Approximation Package
#
# Runs comprehensive verification tests:
# 1. Exact distribution verification (brute force vs recursive)
# 2. Paper table value comparison
# 3. Simulation with large iteration counts
# 4. R cross-validation (optional)
#
# Usage:
#   ./validate.sh           # Run Python verifications
#   ./validate.sh --all     # Run Python + R verifications
#   ./validate.sh --exact   # Only exact distribution tests
#   ./validate.sh --sim     # Only simulation tests
#   ./validate.sh --r       # Only R cross-validation
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_header() {
    echo ""
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
    echo ""
}

print_success() {
    echo -e "${GREEN}[PASS]${NC} $1"
}

print_failure() {
    echo -e "${RED}[FAIL]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

# Check Python environment
check_python() {
    print_header "Checking Python Environment"

    if ! command -v python3 &> /dev/null; then
        print_failure "Python3 not found"
        exit 1
    fi
    print_success "Python3 found: $(python3 --version)"

    # Check required packages
    python3 -c "import numpy" 2>/dev/null && print_success "numpy installed" || print_failure "numpy not installed"
    python3 -c "import scipy" 2>/dev/null && print_success "scipy installed" || print_failure "scipy not installed"

    # Check our package
    cd "$PROJECT_ROOT"
    python3 -c "from kw_approx import ExactDistribution" 2>/dev/null && print_success "kw_approx package available" || print_failure "kw_approx package not available"
}

# Run exact distribution verification
run_exact_verification() {
    print_header "Running Exact Distribution Verification"

    cd "$PROJECT_ROOT"
    python3 "$SCRIPT_DIR/verify_exact.py"

    if [ $? -eq 0 ]; then
        print_success "Exact distribution verification completed"
    else
        print_failure "Exact distribution verification failed"
        return 1
    fi
}

# Run simulation verification
run_simulation_verification() {
    print_header "Running Simulation Verification"

    cd "$PROJECT_ROOT"
    python3 "$SCRIPT_DIR/verify_simulation.py"

    if [ $? -eq 0 ]; then
        print_success "Simulation verification completed"
    else
        print_failure "Simulation verification failed"
        return 1
    fi
}

# Run R cross-validation
run_r_verification() {
    print_header "Running R Cross-Validation"

    if ! command -v Rscript &> /dev/null; then
        print_warning "R not found, skipping R verification"
        return 0
    fi

    print_success "R found: $(R --version | head -1)"

    cd "$SCRIPT_DIR"
    Rscript verify_r_comparison.R

    if [ $? -eq 0 ]; then
        print_success "R cross-validation completed"
    else
        print_warning "R cross-validation had issues"
        return 0  # Don't fail the whole script for R issues
    fi
}

# Generate verification report
generate_report() {
    print_header "Generating Verification Report"

    REPORT_FILE="$PROJECT_ROOT/result/verification_report.txt"

    {
        echo "=============================================="
        echo "KRUSKAL-WALLIS VERIFICATION REPORT"
        echo "Generated: $(date)"
        echo "=============================================="
        echo ""
        echo "1. EXACT DISTRIBUTION VERIFICATION"
        echo "   - Brute force vs recursive algorithm: TESTED"
        echo "   - Paper Table 4.1 (3,3,3): COMPARED"
        echo "   - Paper Table 4.2 (three groups): COMPARED"
        echo "   - Paper Table 4.6 (four groups): COMPARED"
        echo ""
        echo "2. SIMULATION VERIFICATION"
        echo "   - Iterations: 200,000"
        echo "   - Seed: 42 (reproducible)"
        echo "   - Confidence intervals: 95% Wilson score"
        echo "   - Convergence analysis: PERFORMED"
        echo ""
        echo "3. CROSS-VALIDATION"
        echo "   - SciPy kruskal: COMPARED"
        echo "   - R coin package: OPTIONAL"
        echo ""
        echo "=============================================="
    } > "$REPORT_FILE"

    print_success "Report saved to: $REPORT_FILE"
}

# Quick test (small cases only)
run_quick_test() {
    print_header "Running Quick Test"

    cd "$PROJECT_ROOT"
    python3 << 'EOF'
import sys
sys.path.insert(0, '.')

from kw_approx.exact import ExactDistribution

# Test (3,3,3) case
sample_sizes = [3, 3, 3]
exact = ExactDistribution(sample_sizes)

cv, ep = exact.critical_value(0.10)
expected_cv = 4.622222
expected_ep = 0.10

cv_match = abs(cv - expected_cv) < 0.001
ep_match = abs(ep - expected_ep) < 0.001

if cv_match and ep_match:
    print(f"Quick test PASSED: CV={cv:.6f}, EP={ep:.6f}")
    sys.exit(0)
else:
    print(f"Quick test FAILED: CV={cv:.6f} (expected {expected_cv}), EP={ep:.6f} (expected {expected_ep})")
    sys.exit(1)
EOF
}

# Main
main() {
    echo ""
    echo "=============================================="
    echo "KRUSKAL-WALLIS VALIDATION SUITE"
    echo "=============================================="
    echo ""

    case "${1:-all}" in
        --quick)
            run_quick_test
            ;;
        --exact)
            check_python
            run_exact_verification
            ;;
        --sim)
            check_python
            run_simulation_verification
            ;;
        --r)
            run_r_verification
            ;;
        --all)
            check_python
            run_exact_verification
            run_simulation_verification
            run_r_verification
            generate_report
            ;;
        --help|-h)
            echo "Usage: $0 [option]"
            echo ""
            echo "Options:"
            echo "  --quick   Quick test (small cases only)"
            echo "  --exact   Exact distribution tests"
            echo "  --sim     Simulation tests"
            echo "  --r       R cross-validation"
            echo "  --all     All tests (default)"
            echo "  --help    Show this help"
            ;;
        *)
            check_python
            run_exact_verification
            run_simulation_verification
            generate_report
            ;;
    esac

    print_header "Validation Complete"
}

main "$@"
