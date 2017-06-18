/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003 Ferdinando Ametrano
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2015 Peter Caspers

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

#include <ql/math/statistics/incrementalstatistics.hpp>
#include <iomanip>
#include <numeric>

namespace QuantLib {

    Size IncrementalStatistics::samples() const {
        return data_.size();
    }

    Real IncrementalStatistics::weightSum() const {
        return std::accumulate(data_.begin(), data_.end(), 0.0,
                               [](Real x, std::pair<Real, Real> p) { return x + p.second; });
    }

    Real IncrementalStatistics::mean() const {
        QL_REQUIRE(weightSum() > 0.0, "sampleWeight_= 0, unsufficient");
        Real weighted_samples = std::accumulate(data_.begin(), data_.end(), 0.0,
                                                [](Real x, std::pair<Real, Real> p) { return x + p.first * p.second; });
        return weighted_samples / weightSum();
    }

    Real IncrementalStatistics::variance() const {
        QL_REQUIRE(weightSum() > 0.0, "sampleWeight_= 0, unsufficient");
        QL_REQUIRE(samples() > 1, "sample number <= 1, unsufficient");
        Real n = static_cast<Real>(samples());
        Real mn = mean();
        Real co = n / weightSum() / (n - 1.0);
        return co * std::accumulate(data_.begin(), data_.end(), 0.0,
                                   [mn](Real x, std::pair<Real, Real> p) {
                                       return x + pow((p.first - mn), 2.0) * p.second; });

    }

    Real IncrementalStatistics::standardDeviation() const {
        return std::sqrt(variance());
    }

    Real IncrementalStatistics::errorEstimate() const {
        return std::sqrt(variance() / (samples()));
    }

    Real IncrementalStatistics::skewness() const {
        QL_REQUIRE(samples() > 2, "sample number <= 2, unsufficient");
        Real n = static_cast<Real>(samples());
        Real sigma3 = pow(variance(), 1.5);
        Real mn = mean();
        Real co = n * n / weightSum() / sigma3 / (n - 1.0) / (n - 2.0);

        return co * std::accumulate(data_.begin(), data_.end(), 0.0,
                                       [mn](Real x, std::pair<Real, Real> p) {
                                           return x + pow((p.first - mn), 3.0) * p.second; });
    }

    Real IncrementalStatistics::kurtosis() const {
        QL_REQUIRE(samples() > 3,
                   "sample number <= 3, unsufficient");
        Real n = static_cast<Real>(samples());
        Real sigma4 = pow(variance(), 2.0);
        Real mn = mean();
        Real co1 = n * n * (n + 1.0)/ weightSum() / sigma4 /(n - 1.0) / (n - 2.0) / (n - 3.0);
        Real co0 = 3.0 * (n - 1.0) * (n - 1.0) / (n - 2.0) / (n - 3.0);

        return co1 * std::accumulate(data_.begin(), data_.end(), 0.0,
                                                   [mn](Real x, std::pair<Real, Real> p) {
                                                       return x + pow((p.first - mn), 4.0) * p.second;
                                                   }) - co0;
    }

    Real IncrementalStatistics::min() const {
        QL_REQUIRE(samples() > 0, "empty sample set");
        return std::min_element(data_.begin(), data_.end())->first;
    }

    Real IncrementalStatistics::max() const {
        QL_REQUIRE(samples() > 0, "empty sample set");
        return std::max_element(data_.begin(), data_.end())->first;
    }

    Size IncrementalStatistics::downsideSamples() const {
        return downsideData_.size();
    }

    Real IncrementalStatistics::downsideWeightSum() const {
        return std::accumulate(downsideData_.begin(), downsideData_.end(), 0.0,
                               [](Real x, std::pair<Real, Real> p) { return x + p.second; });
    }

    Real IncrementalStatistics::downsideVariance() const {
        QL_REQUIRE(downsideWeightSum() > 0.0, "sampleWeight_= 0, unsufficient");
        QL_REQUIRE(downsideSamples() > 1, "sample number <= 1, unsufficient");
        Real n = static_cast<Real >(downsideSamples());
        Real co = n / downsideWeightSum() / (n-1.0);
        return co * std::accumulate(downsideData_.begin(), downsideData_.end(), 0.0,
                                   [](Real x, std::pair<Real, Real> p) {
                                       return x + pow(p.first, 2.0) * p.second; });

    }

    Real IncrementalStatistics::downsideDeviation() const {
        return sqrt(downsideVariance());
    }

    void IncrementalStatistics::add(Real value, Real valueWeight) {
        QL_REQUIRE(valueWeight >= 0.0, "negative weight (" << valueWeight
                                                           << ") not allowed");
        data_.emplace_back(std::make_pair(value, valueWeight));
        if (value < 0.0)
            downsideData_.emplace_back(std::make_pair(value, valueWeight));
    }

    void IncrementalStatistics::reset() {
        data_.clear();
        downsideData_.clear();
    }

}
