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
#include <unistd.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "AntennaChecker.hpp"
#include "IdbDesign.h"
#include "IdbGeometry.h"
#include "IdbInstance.h"
#include "IdbNet.h"
#include "IdbRegularWire.h"
#include "MJInterface.hpp"
#include "idm.h"
#include "json.hpp"

namespace {

void require(bool condition, const std::string& message)
{
  if (!condition) {
    throw std::runtime_error(message);
  }
}

std::string readText(const std::filesystem::path& path)
{
  std::ifstream in(path);
  require(in.is_open(), "Failed to read file: " + path.string());
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

void writeText(const std::filesystem::path& path, const std::string& text)
{
  std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path);
  require(out.is_open(), "Failed to write file: " + path.string());
  out << text;
}

bool rectsViolateSpacing(const idb::IdbRect& a, const idb::IdbRect& b, int32_t spacing)
{
  const int64_t sp = spacing;
  return static_cast<int64_t>(a.get_low_x()) < static_cast<int64_t>(b.get_high_x()) + sp &&
         static_cast<int64_t>(a.get_high_x()) + sp > static_cast<int64_t>(b.get_low_x()) &&
         static_cast<int64_t>(a.get_low_y()) < static_cast<int64_t>(b.get_high_y()) + sp &&
         static_cast<int64_t>(a.get_high_y()) + sp > static_cast<int64_t>(b.get_low_y());
}

const char* kTechLef = R"(VERSION 5.8 ;
DIVIDERCHAR "/" ;
BUSBITCHARS "[]" ;
UNITS
  DATABASE MICRONS 1000 ;
END UNITS
MANUFACTURINGGRID 0.001 ;
LAYER M1
  TYPE ROUTING ;
  DIRECTION HORIZONTAL ;
  PITCH 2 ;
  WIDTH 1 ;
  SPACING 1 ;
  THICKNESS 0.1 ;
  ANTENNAAREARATIO 2 ;
  ANTENNACUMAREARATIO 3 ;
  ANTENNADIFFAREARATIO 20 ;
  ANTENNACUMDIFFAREARATIO 30 ;
END M1
LAYER V1
  TYPE CUT ;
  SPACING 0.1 ;
  ANTENNAAREARATIO 0.05 ;
  ANTENNADIFFAREARATIO 2 ;
END V1
LAYER M2
  TYPE ROUTING ;
  DIRECTION VERTICAL ;
  PITCH 2 ;
  WIDTH 1 ;
  SPACING 1 ;
  THICKNESS 0.1 ;
  ANTENNAAREARATIO 2 ;
  ANTENNACUMAREARATIO 7 ;
  ANTENNADIFFAREARATIO 20 ;
  ANTENNACUMDIFFAREARATIO 30 ;
END M2
VIA V12 DEFAULT
  LAYER M1 ;
    RECT -0.5 -0.5 0.5 0.5 ;
  LAYER V1 ;
    RECT -0.1 -0.1 0.1 0.1 ;
  LAYER M2 ;
    RECT -0.5 -0.5 0.5 0.5 ;
END V12
SITE CORE
  CLASS CORE ;
  SIZE 2 BY 6 ;
  SYMMETRY X Y ;
END CORE
END LIBRARY
)";

const char* kCellLef = R"(VERSION 5.8 ;
DIVIDERCHAR "/" ;
BUSBITCHARS "[]" ;
MACRO GATE_M1
  CLASS CORE ;
  ORIGIN 0 0 ;
  SIZE 2 BY 6 ;
  SITE CORE ;
  PIN A
    DIRECTION INPUT ;
    USE SIGNAL ;
    ANTENNAGATEAREA 1 ;
    PORT
      LAYER M1 ;
        RECT 0.5 2.5 1.5 3.5 ;
    END
  END A
END GATE_M1
MACRO GATE_M2
  CLASS CORE ;
  ORIGIN 0 0 ;
  SIZE 2 BY 6 ;
  SITE CORE ;
  PIN A
    DIRECTION INPUT ;
    USE SIGNAL ;
    ANTENNAGATEAREA 1 ;
    PORT
      LAYER M2 ;
        RECT 0.5 2.5 1.5 3.5 ;
    END
  END A
END GATE_M2
MACRO FILL_CELL
  CLASS CORE SPACER ;
  ORIGIN 0 0 ;
  SIZE 2 BY 6 ;
  SITE CORE ;
END FILL_CELL
MACRO ANTENNA_DIODE
  CLASS CORE ;
  ORIGIN 0 0 ;
  SIZE 2 BY 6 ;
  SITE CORE ;
  PIN DIODE
    DIRECTION INPUT ;
    USE SIGNAL ;
    ANTENNADIFFAREA 5 ;
    PORT
      LAYER M1 ;
        RECT 0.5 2.5 1.5 3.5 ;
    END
  END DIODE
END ANTENNA_DIODE
END LIBRARY
)";

class AntennaRepairTest
{
 public:
  void run()
  {
    setenv("ECC_LOGGER_THROW_ON_ERROR", "1", 1);
    std::filesystem::create_directories(_root);
    writeText(_tech_lef, kTechLef);
    writeText(_cell_lef, kCellLef);

    testM1ObstacleDetourRepair();
    testViaBridgeOverBlockedM1();
    testLayerQualifiedDiodeDiscovery();
    testCongestedNoSiteRollback();

    dmInst->reset();
    std::filesystem::remove_all(_root);
  }

 private:
  std::filesystem::path _root = std::filesystem::temp_directory_path() / ("test_imj_antenna_repair_" + std::to_string(getpid()));
  std::filesystem::path _tech_lef = _root / "tech.lef";
  std::filesystem::path _cell_lef = _root / "cells.lef";

  void loadCase(const std::string& name, const std::string& def_text)
  {
    std::filesystem::path case_dir = _root / name;
    std::filesystem::path def_path = case_dir / "input.def";
    writeText(def_path, def_text);

    dmInst->reset();
    require(dmInst->readLef(std::vector<std::string>{_tech_lef.string()}, true), "read tech LEF");
    require(dmInst->readLef(std::vector<std::string>{_cell_lef.string()}), "read cell LEF");
    require(dmInst->readDef(def_path.string()), "read DEF for " + name);
  }

  void verifyRoundTripClean(const std::filesystem::path& def_out, const std::string& mj_tmp)
  {
    require(dmInst->saveDef(def_out.string()), "save repaired DEF");
    dmInst->reset();
    require(dmInst->readLef(std::vector<std::string>{_tech_lef.string()}, true), "reload tech LEF");
    require(dmInst->readLef(std::vector<std::string>{_cell_lef.string()}), "reload cell LEF");
    require(dmInst->readDef(def_out.string()), "reload repaired DEF");

    MJI.initMJ({{"-temp_directory_path", mj_tmp}});
    std::filesystem::path check_json = def_out.parent_path() / "recheck.json";
    MJI.checkAntenna({{"-step_json", check_json.string()}});
    MJI.destroyMJ();

    nlohmann::json parsed = nlohmann::json::parse(readText(check_json));
    require(parsed.at("antenna_check").at("violation_count").get<int>() == 0,
            "expected 0 violations after DEF save/reload round-trip, got: " + parsed.dump());
  }

  void testM1ObstacleDetourRepair()
  {
    // Sites between x=12000 and x=18000: blocker occupies x=14000..16000 + foreign M1 obstacle at x=16000, y=[2500..3500].
    // Nearest free site left of blocker is (12000, 0), forcing an M1 detour around the foreign M1 wire.
    const std::string def = R"(VERSION 5.8 ;
DIVIDERCHAR "/" ;
BUSBITCHARS "[]" ;
DESIGN repair_m1_detour ;
UNITS DISTANCE MICRONS 1000 ;
DIEAREA ( 0 0 ) ( 24000 12000 ) ;
ROW ROW_0 CORE 0 0 N DO 12 BY 1 STEP 2000 0 ;
COMPONENTS 7 ;
- u_block_0 FILL_CELL + PLACED ( 0 0 ) N ;
- u_block_1 FILL_CELL + PLACED ( 2000 0 ) N ;
- u_block_2 FILL_CELL + PLACED ( 4000 0 ) N ;
- u_block_3 FILL_CELL + PLACED ( 6000 0 ) N ;
- u_block_4 FILL_CELL + PLACED ( 8000 0 ) N ;
- u_block_6 FILL_CELL + PLACED ( 14000 0 ) N ;
- u_gate GATE_M1 + PLACED ( 18000 0 ) N ;
END COMPONENTS
NETS 2 ;
- viol_net ( u_gate A )
  + USE SIGNAL
  + ROUTED M1 ( 19000 3000 ) ( 19000 10000 ) ;
- obs_net
  + USE SIGNAL
  + ROUTED M1 ( 16000 2500 ) ( 16000 3500 ) ;
END NETS
END DESIGN
)";

    loadCase("m1_detour", def);
    std::filesystem::path case_dir = _root / "m1_detour";
    std::filesystem::path step_json = case_dir / "repair.step.json";
    std::filesystem::path rpt_path = case_dir / "antenna_repair.rpt";

    MJI.initMJ({{"-temp_directory_path", (case_dir / "mj").string()}});
    MJI.repairAntenna({{"-max_iter", 2},
                       {"-search_radius_um", 30.0},
                       {"-report_path", rpt_path.string()},
                       {"-step_json", step_json.string()}});
    MJI.destroyMJ();

    nlohmann::json payload = nlohmann::json::parse(readText(step_json));
    const auto& ar = payload.at("antenna_repair");
    require(ar.at("violations_before").get<int>() > 0, "expected violations before M1 detour repair");
    require(ar.at("violations_after").get<int>() == 0, "expected 0 violations after M1 detour repair: " + payload.dump());
    require(ar.at("diodes_inserted").get<int>() >= 1, "expected diode insertion");

    // Verify no spacing violation or short with foreign M1 obstacle net obs_net.
    idb::IdbDesign* design = dmInst->get_idb_design();
    idb::IdbNet* viol_net = design->get_net_list()->find_net("viol_net");
    idb::IdbNet* obs_net = design->get_net_list()->find_net("obs_net");
    require(viol_net && obs_net, "expected both nets");
    idb::IdbRect obs_rect = obs_net->get_wire_list()->get_wire_list()[0]->get_segment_list()[0]->get_segment_rect();
    for (idb::IdbRegularWire* w : viol_net->get_wire_list()->get_wire_list()) {
      for (idb::IdbRegularWireSegment* s : w->get_segment_list()) {
        if (s->is_wire() && s->get_layer() && s->get_layer()->get_name() == "M1") {
          idb::IdbRect wr = s->get_segment_rect();
          require(!rectsViolateSpacing(wr, obs_rect, 1000), "M1 stub violated spacing against foreign M1 obstacle");
        }
      }
    }

    verifyRoundTripClean(case_dir / "repaired.def", (case_dir / "mj_recheck").string());
  }

  void testViaBridgeOverBlockedM1()
  {
    // GATE_M2 has a violating M2 antenna, and a full-height foreign M1 wall spans x=16000, y=[500..11500].
    // The diode pin is on M1 at x<=14000, so the stub must insert V12 vias and bridge across on M2.
    const std::string def = R"(VERSION 5.8 ;
DIVIDERCHAR "/" ;
BUSBITCHARS "[]" ;
DESIGN repair_via_bridge ;
UNITS DISTANCE MICRONS 1000 ;
DIEAREA ( 0 0 ) ( 24000 12000 ) ;
ROW ROW_0 CORE 0 0 N DO 12 BY 1 STEP 2000 0 ;
COMPONENTS 8 ;
- u_block_0 FILL_CELL + PLACED ( 0 0 ) N ;
- u_block_1 FILL_CELL + PLACED ( 2000 0 ) N ;
- u_block_2 FILL_CELL + PLACED ( 4000 0 ) N ;
- u_block_3 FILL_CELL + PLACED ( 6000 0 ) N ;
- u_block_4 FILL_CELL + PLACED ( 8000 0 ) N ;
- u_block_6 FILL_CELL + PLACED ( 14000 0 ) N ;
- u_block_7 FILL_CELL + PLACED ( 16000 0 ) N ;
- u_gate GATE_M2 + PLACED ( 18000 0 ) N ;
END COMPONENTS
NETS 2 ;
- viol_m2 ( u_gate A )
  + USE SIGNAL
  + ROUTED M2 ( 19000 3000 ) ( 19000 11000 ) ;
- m1_wall
  + USE SIGNAL
  + ROUTED M1 ( 16000 500 ) ( 16000 11500 ) ;
END NETS
END DESIGN
)";

    loadCase("via_bridge", def);
    std::filesystem::path case_dir = _root / "via_bridge";
    std::filesystem::path step_json = case_dir / "repair.step.json";

    MJI.initMJ({{"-temp_directory_path", (case_dir / "mj").string()}});
    MJI.repairAntenna({{"-diode", "ANTENNA_DIODE"},
                       {"-max_iter", 2},
                       {"-search_radius_um", 30.0},
                       {"-step_json", step_json.string()}});
    MJI.destroyMJ();

    nlohmann::json payload = nlohmann::json::parse(readText(step_json));
    const auto& ar = payload.at("antenna_repair");
    require(ar.at("violations_before").get<int>() > 0, "expected violations before via-bridge repair");
    require(ar.at("violations_after").get<int>() == 0, "expected 0 violations after via-bridge repair: " + payload.dump());
    require(ar.at("diodes_inserted").get<int>() == 1, "expected 1 diode inserted");
    require(payload.at("placed_diodes").at(0).at("vias_added").get<int>() >= 1, "expected via bridge across blocked M1");

    idb::IdbDesign* design = dmInst->get_idb_design();
    idb::IdbNet* viol_net = design->get_net_list()->find_net("viol_m2");
    idb::IdbNet* wall_net = design->get_net_list()->find_net("m1_wall");
    idb::IdbRect wall_rect = wall_net->get_wire_list()->get_wire_list()[0]->get_segment_list()[0]->get_segment_rect();
    for (idb::IdbRegularWire* w : viol_net->get_wire_list()->get_wire_list()) {
      for (idb::IdbRegularWireSegment* s : w->get_segment_list()) {
        if (s->is_wire() && s->get_layer() && s->get_layer()->get_name() == "M1") {
          require(!rectsViolateSpacing(s->get_segment_rect(), wall_rect, 1000),
                  "M1 segment shorted or violated spacing against full-height M1 wall");
        }
      }
    }

    verifyRoundTripClean(case_dir / "repaired.def", (case_dir / "mj_recheck").string());
  }

  void testLayerQualifiedDiodeDiscovery()
  {
    const std::string layer_cell_lef = R"(VERSION 5.8 ;
DIVIDERCHAR "/" ;
BUSBITCHARS "[]" ;
MACRO GATE_OX1
  CLASS CORE ;
  ORIGIN 0 0 ;
  SIZE 2 BY 6 ;
  SITE CORE ;
  PIN A
    DIRECTION INPUT ;
    USE SIGNAL ;
    ANTENNAMODEL OXIDE1 ;
      ANTENNAGATEAREA 1 LAYER M1 ;
    PORT
      LAYER M1 ;
        RECT 0.5 2.5 1.5 3.5 ;
    END
  END A
END GATE_OX1
MACRO ANTENNA_M1_DIO
  CLASS CORE ;
  ORIGIN 0 0 ;
  SIZE 2 BY 6 ;
  SITE CORE ;
  PIN A
    DIRECTION INPUT ;
    USE SIGNAL ;
    ANTENNADIFFAREA 5 LAYER M1 ;
    PORT
      LAYER M1 ;
        RECT 0.5 2.5 1.5 3.5 ;
    END
  END A
END ANTENNA_M1_DIO
END LIBRARY
)";
    const std::string def = R"(VERSION 5.8 ;
DIVIDERCHAR "/" ;
BUSBITCHARS "[]" ;
DESIGN repair_layer_diode ;
UNITS DISTANCE MICRONS 1000 ;
DIEAREA ( 0 0 ) ( 16000 6000 ) ;
ROW ROW_0 CORE 0 0 N DO 8 BY 1 STEP 2000 0 ;
COMPONENTS 1 ;
- u_gate GATE_OX1 + PLACED ( 8000 0 ) N ;
END COMPONENTS
NETS 1 ;
- viol_net ( u_gate A )
  + USE SIGNAL
  + ROUTED M1 ( 9000 3000 ) ( 9000 6000 ) ;
END NETS
END DESIGN
)";
    std::filesystem::path case_dir = _root / "layer_diode";
    std::filesystem::path lcell = case_dir / "layer_cells.lef";
    std::filesystem::path ldef = case_dir / "input.def";
    std::filesystem::path step_json = case_dir / "repair.step.json";
    writeText(lcell, layer_cell_lef);
    writeText(ldef, def);

    dmInst->reset();
    require(dmInst->readLef(std::vector<std::string>{_tech_lef.string()}, true), "read tech LEF");
    require(dmInst->readLef(std::vector<std::string>{lcell.string()}), "read layer-qualified cell LEF");
    require(dmInst->readDef(ldef.string()), "read DEF");

    MJI.initMJ({{"-temp_directory_path", (case_dir / "mj").string()}});
    MJI.repairAntenna({{"-max_iter", 2},
                       {"-search_radius_um", 20.0},
                       {"-step_json", step_json.string()}});
    MJI.destroyMJ();

    nlohmann::json payload = nlohmann::json::parse(readText(step_json));
    const auto& ar = payload.at("antenna_repair");
    require(ar.at("violations_before").get<int>() > 0, "expected violation before layer-qualified diode repair");
    require(ar.at("violations_after").get<int>() == 0, "expected 0 violations after layer-qualified diode repair: " + payload.dump());
    require(ar.at("diodes_inserted").get<int>() == 1, "expected 1 auto-discovered ANTENNA_M1_DIO diode");
  }

  void testCongestedNoSiteRollback()
  {
    const std::string def = R"(VERSION 5.8 ;
DIVIDERCHAR "/" ;
BUSBITCHARS "[]" ;
DESIGN repair_full_row ;
UNITS DISTANCE MICRONS 1000 ;
DIEAREA ( 0 0 ) ( 8000 6000 ) ;
ROW ROW_0 CORE 0 0 N DO 4 BY 1 STEP 2000 0 ;
COMPONENTS 4 ;
- u_fill0 FILL_CELL + PLACED ( 0 0 ) N ;
- u_fill1 FILL_CELL + PLACED ( 2000 0 ) N ;
- u_gate GATE_M1 + PLACED ( 4000 0 ) N ;
- u_fill3 FILL_CELL + PLACED ( 6000 0 ) N ;
END COMPONENTS
NETS 1 ;
- viol_net ( u_gate A )
  + USE SIGNAL
  + ROUTED M1 ( 5000 3000 ) ( 5000 6000 ) ;
END NETS
END DESIGN
)";

    loadCase("congested", def);
    std::filesystem::path case_dir = _root / "congested";
    std::filesystem::path before_def = case_dir / "before.def";
    std::filesystem::path after_def = case_dir / "after.def";
    std::filesystem::path step_json = case_dir / "repair.step.json";

    require(dmInst->saveDef(before_def.string()), "save baseline DEF");

    MJI.initMJ({{"-temp_directory_path", (case_dir / "mj").string()}});
    MJI.repairAntenna({{"-step_json", step_json.string()}});
    MJI.destroyMJ();

    require(dmInst->saveDef(after_def.string()), "save post-repair DEF");
    require(readText(before_def) == readText(after_def), "unrepaired congested DEF must remain unchanged");

    nlohmann::json payload = nlohmann::json::parse(readText(step_json));
    const auto& ar = payload.at("antenna_repair");
    require(ar.at("diodes_inserted").get<int>() == 0, "expected 0 diodes when row is full");
    require(ar.at("nets_unrepaired").get<int>() == 1, "expected 1 unrepaired net when row is full");
    require(ar.at("violations_before").get<int>() == ar.at("violations_after").get<int>(),
            "expected identical violation count");
  }
};

}  // namespace

int main()
{
  try {
    AntennaRepairTest test;
    test.run();
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "test_antenna_repair failed: " << ex.what() << std::endl;
    return 1;
  }
}
