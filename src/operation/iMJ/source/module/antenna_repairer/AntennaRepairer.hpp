// ***************************************************************************************
// Copyright (c) 2023-2025 Peng Cheng Laboratory
// Copyright (c) 2023-2025 Institute of Computing Technology, Chinese Academy of Sciences
// Copyright (c) 2023-2025 Beijing Institute of Open Source Chip
//
// iEDA is licensed under Mulan PSL v2.
// You can use this software according to the terms and conditions of the Mulan PSL v2.
// You may obtain a copy of Mulan PSL v2 at:
// http://license.coscl.org.cn/MulanPSL2
//
// THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND,
// EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT,
// MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
//
// See the Mulan PSL v2 for more details.
// ***************************************************************************************
#pragma once

#include <any>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "ACViolation.hpp"
#include "IdbDesign.h"
#include "IdbGeometry.h"

namespace imj {

#define MJAR (imj::AntennaRepairer::getInst())

struct ARPlacedDiode
{
  std::string instance_name;
  std::string master_name;
  std::string net_name;
  std::string target_pin;
  int32_t x_dbu = 0;
  int32_t y_dbu = 0;
  int wire_segments_added = 0;
  int vias_added = 0;
};

struct ARSummary
{
  int violations_before = 0;
  int violations_after = 0;
  int diodes_inserted = 0;
  int nets_repaired = 0;
  int nets_unrepaired = 0;
  int iterations = 0;
  std::vector<std::string> diode_masters_used;
  std::vector<ARPlacedDiode> placed_diodes;
};

class AntennaRepairer
{
 public:
  static void initInst();
  static AntennaRepairer& getInst();
  static void destroyInst();

  ARSummary repair(std::map<std::string, std::any> config_map);

 private:
  static AntennaRepairer* _ar_instance;

  AntennaRepairer() = default;
  AntennaRepairer(const AntennaRepairer& other) = delete;
  AntennaRepairer(AntennaRepairer&& other) = delete;
  ~AntennaRepairer() = default;

  AntennaRepairer& operator=(const AntennaRepairer& other) = delete;
  AntennaRepairer& operator=(AntennaRepairer&& other) = delete;
};

}  // namespace imj
