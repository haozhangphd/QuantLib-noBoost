/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003 RiskMap srl
 Copyright (C) 2012 StatPro Italia srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include "utilities.hpp"
#include <ql/math/solvers1d/brent.hpp>
#include <ql/math/solvers1d/bisection.hpp>
#include <ql/math/solvers1d/falseposition.hpp>
#include <ql/math/solvers1d/ridder.hpp>
#include <ql/math/solvers1d/secant.hpp>
#include <ql/math/solvers1d/newton.hpp>
#include <ql/math/solvers1d/halley.hpp>
#include <ql/math/solvers1d/finitedifferencenewtonsafe.hpp>

using namespace QuantLib;


namespace {

    class F1 {
    public:
        Real operator()(Real x) const { return x * x - 1.0; }

        Real derivative(Real x) const { return 2.0 * x; }

        Real derivative2(Real x) const { return 0.5 / x; }
    };

    class F2 {
    public:
        Real operator()(Real x) const { return 1.0 - x * x; }

        Real derivative(Real x) const { return -2.0 * x; }

        Real derivative2(Real x) const { return 0.5 / x; }
    };

    class F3 {
    public:
        Real operator()(Real x) const { return std::atan(x - 1); }

        Real derivative(Real x) const { return 1.0 / (1.0 + (x - 1.0) * (x - 1.0)); }

        Real derivative2(Real x) const { return -x / (1 + x * x); }
    };

    template<class S, class F>
    void test_not_bracketed(const S &solver, const std::string &name,
                            const F &f, Real guess) {
        Real accuracy[] = {1.0e-4, 1.0e-6, 1.0e-8};
        Real expected = 1.0;
        for (Size i = 0; i < LENGTH(accuracy); i++) {
            Real root = solver.solve(f, accuracy[i], guess, 0.1);
            if (std::fabs(root - expected) > accuracy[i]) {
                FAIL(name << " solver (not bracketed):\n"
                          << "    expected:   " << expected << "\n"
                          << "    calculated: " << root << "\n"
                          << "    accuracy:   " << accuracy[i]);
            }
        }
    }

    template<class S, class F>
    void test_bracketed(const S &solver, const std::string &name,
                        const F &f, Real guess) {
        Real accuracy[] = {1.0e-4, 1.0e-6, 1.0e-8};
        Real expected = 1.0;
        for (Size i = 0; i < LENGTH(accuracy); i++) {
            // guess on the left side of the root, increasing function
            Real root = solver.solve(f, accuracy[i], guess, 0.0, 2.0);
            if (std::fabs(root - expected) > accuracy[i]) {
                FAIL(name << " solver (bracketed):\n"
                          << "    expected:   " << expected << "\n"
                          << "    calculated: " << root << "\n"
                          << "    accuracy:   " << accuracy[i]);
            }
        }
    }

    class Probe {
    public:
        Probe(Real &result, Real offset)
                : result_(result), previous_(result), offset_(offset) {}

        Real operator()(Real x) const {
            result_ = x;
            return previous_ + offset_ - x * x;
        }

        Real derivative(Real x) const { return 2.0 * x; }

        Real derivative2(Real x) const { return 0.5 / x; }

    private:
        Real &result_;
        Real previous_;
        Real offset_;
    };

    template<class S>
    void test_last_call_with_root(const S &solver, const std::string &name,
                                  bool bracketed, Real accuracy) {

        Real mins[] = {3.0, 2.25, 1.5, 1.0};
        Real maxs[] = {7.0, 5.75, 4.5, 3.0};
        Real steps[] = {0.2, 0.2, 0.1, 0.1};
        Real offsets[] = {25.0, 11.0, 5.0, 1.0};
        Real guesses[] = {4.5, 4.5, 2.5, 2.5};
        //Real expected[] = { 5.0, 4.0, 3.0, 2.0 };

        Real argument = 0.0;
        Real result;

        for (Size i = 0; i < 4; ++i) {
            if (bracketed) {
                result = solver.solve(Probe(argument, offsets[i]), accuracy,
                                      guesses[i], mins[i], maxs[i]);
            } else {
                result = solver.solve(Probe(argument, offsets[i]), accuracy,
                                      guesses[i], steps[i]);
            }

            Real error = std::fabs(result - argument);
            // no floating-point comparison: the solver should have
            // called the function with the very same value it's
            // returning
            if (result != argument) {
                FAIL(name << " solver ("
                          << (bracketed ? "" : "not ")
                          << "bracketed):\n"
                          << "    index:   " << i << "\n"
                          << "    expected:   " << result << "\n"
                          << "    calculated: " << argument << "\n"
                          << "    error: " << error);
            }
        }
    }

    template<class S>
    void test_solver(const S &solver, const std::string &name, Real accuracy) {
        // guess on the left side of the root, increasing function
        test_not_bracketed(solver, name, F1(), 0.5);
        test_bracketed(solver, name, F1(), 0.5);
        // guess on the right side of the root, increasing function
        test_not_bracketed(solver, name, F1(), 1.5);
        test_bracketed(solver, name, F1(), 1.5);
        // guess on the left side of the root, decreasing function
        test_not_bracketed(solver, name, F2(), 0.5);
        test_bracketed(solver, name, F2(), 0.5);
        // guess on the right side of the root, decreasing function
        test_not_bracketed(solver, name, F2(), 1.5);
        test_bracketed(solver, name, F2(), 1.5);
        // situation where bisection is used in the finite difference
        // newton solver as the first step and where the initial
        // guess is equal to the next estimate (which causes an infinite
        // derivative if we do not handle this case with special care)
        test_not_bracketed(solver, name, F3(), 1.00001);
        // check that the last function call is made with the root value
        if (accuracy != Null<Real>()) {
            test_last_call_with_root(solver, name, false, accuracy);
            test_last_call_with_root(solver, name, true, accuracy);
        }
    }

}


TEST_CASE("Solver1D_Brent", "[Solver1D]") {
    INFO("Testing Brent solver...");
    test_solver(Brent(), "Brent", 1.0e-6);
}

TEST_CASE("Solver1D_Bisection", "[Solver1D]") {
    INFO("Testing bisection solver...");
    test_solver(Bisection(), "Bisection", 1.0e-6);
}

TEST_CASE("Solver1D_FalsePosition", "[Solver1D]") {
    INFO("Testing false-position solver...");
    test_solver(FalsePosition(), "FalsePosition", 1.0e-6);
}

TEST_CASE("Solver1D_Newton", "[Solver1D]") {
    INFO("Testing Newton solver...");
    test_solver(Newton(), "Newton", 1.0e-12);
}

TEST_CASE("Solver1D_NewtonSafe", "[Solver1D]") {
    INFO("Testing Newton-safe solver...");
    test_solver(NewtonSafe(), "NewtonSafe", 1.0e-9);
}

TEST_CASE("Solver1D_Halley", "[.]") {
    INFO("Testing Halley solver...");
    test_solver(Halley(), "Halley", 1.0e-12);
}

TEST_CASE("Solver1D_HalleySafe", "[.]") {
    INFO("Testing Halley-safe solver...");
    test_solver(HalleySafe(), "HalleySafe", 1.0e-9);
}

TEST_CASE("Solver1D_FiniteDifferenceNewtonSafe", "[Solver1D]") {
    INFO("Testing finite-difference Newton-safe solver...");
    test_solver(FiniteDifferenceNewtonSafe(), "FiniteDifferenceNewtonSafe", Null<Real>());
}

TEST_CASE("Solver1D_Ridder", "[Solver1D]") {
    INFO("Testing Ridder solver...");
    test_solver(Ridder(), "Ridder", 1.0e-6);
}

TEST_CASE("Solver1D_Secant", "[Solver1D]") {
    INFO("Testing secant solver...");
    test_solver(Secant(), "Secant", 1.0e-6);
}

