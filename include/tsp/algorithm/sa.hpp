/**
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */

#pragma once

#include <random>

#include "math/matrix.hpp"
#include "tsp/algorithm/algorithm.hpp"

namespace tsp::algorithm
{
class SA : public Algorithm
{
    const float kMinimumTemperature = 0.001;

public:
    SA(math::Matrix<uint32_t> distances, float temperature, float multiplier, uint16_t epoch_size);

public:
    Solution Solve() override;

protected:
    bool AcceptSolution(const Solution& solution, float temperature);

    Solution CreateNeighbourSolution(const Solution& value);

    uint16_t CalculateWeight(const Path& value) const;

    static double GenerateRandomInteger(int32_t min, int32_t max);

private:
    Solution solution_;

    const float temperature_;
    const float multiplier_;
    const uint16_t epoch_size_;
};
} // namespace tsp::algorithm