// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <array>
#include <set>
#include <vector>

#include "Args.h"
#include "pre_types.h"

using namespace matfp;

namespace pre_matfp {

void find_feature_edges(const Args& args,
                        const std::vector<Vector3>& input_vertices,
                        const std::vector<Vector3i>& input_faces,
                        std::set<std::array<int, 2>>& s_edges,
                        std::set<std::array<int, 2>>& cc_edges,
                        std::set<int>& corners);

}  // namespace pre_matfp
