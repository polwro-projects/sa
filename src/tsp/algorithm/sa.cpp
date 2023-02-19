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

#include "tsp/algorithm/sa.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tsp::algorithm
{
SA::SA(math::Matrix<uint32_t> distances, float temperature, float multiplier, uint16_t epoch_size)
    : Algorithm{ distances }, temperature_{ temperature }, multiplier_{ multiplier }, epoch_size_{ epoch_size }
{
    solution_.path.resize(distances_.Columns());

    generate(solution_.path.begin(), solution_.path.end(), [n = 0]() mutable { return n++; });
    random_shuffle(solution_.path.begin(), solution_.path.end());

    solution_.weight = CalculateWeight(solution_.path);
}

Algorithm::Solution SA::Solve()
{
    // while the stop criterion is not yet satisfied do
    float current_temperature = temperature_;
    while (current_temperature > kMinimumTemperature)
    {
        for (uint16_t epoch = 0; epoch < epoch_size_; epoch++)
        {
            const auto solution = CreateNeighbourSolution(solution_);
            if (!AcceptSolution(solution, current_temperature))
            {
                continue;
            }
            solution_ = solution;
        }

        // Update the temperature
        current_temperature *= multiplier_;
    }

    return solution_;
}

bool SA::AcceptSolution(const Solution& solution, float temperature)
{
    // If the path is smaller
    if (solution_.weight > solution.weight)
    {
        return true;
    }

    // We can still accept the solution, even if it is worse
    if ((int) (exp((solution_.weight - solution.weight) / temperature) * 100) > GenerateRandomInteger(0, 101))
    {
        return true;
    }
    return false;
}

SA::Solution SA::CreateNeighbourSolution(const Solution& value)
{
    const auto start_index = GenerateRandomInteger(0, distances_.Columns());
    const auto end_index = GenerateRandomInteger(0, distances_.Columns());
    if (start_index == end_index)
    {
        return value;
    }

    Solution solution = value;
    std::swap(solution.path[start_index], solution.path[end_index]);

    solution.weight = CalculateWeight(solution.path);
    return solution;
}

uint16_t SA::CalculateWeight(const Path& value) const
{
    uint16_t weight{};
    for (uint8_t index{}; index < value.size() - 1; ++index)
    {
        weight += distances_(value.at(index), value.at(index + 1));
    }

    return weight;
}

double SA::GenerateRandomInteger(int32_t min, int32_t max)
{
    static std::mt19937 rng{ std::random_device{}() };
    std::uniform_real_distribution<double> dist(min, max);
    return dist(rng);
}
} // namespace tsp::algorithm