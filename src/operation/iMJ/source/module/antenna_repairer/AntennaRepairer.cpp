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
#include "AntennaRepairer.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include "AntennaChecker.hpp"
#include "DataManager.hpp"
#include "IdbBlockages.h"
#include "IdbCellMaster.h"
#include "IdbDesign.h"
#include "IdbInstance.h"
#include "IdbLayer.h"
#include "IdbNet.h"
#include "IdbObs.h"
#include "IdbOrientTransform.h"
#include "IdbPins.h"
#include "IdbRegularWire.h"
#include "IdbRow.h"
#include "IdbSite.h"
#include "IdbSpecialNet.h"
#include "IdbSpecialWire.h"
#include "IdbTerm.h"
#include "IdbTrackGrid.h"
#include "IdbVias.h"
#include "Logger.hpp"
#include "Monitor.hpp"
#include "Utility.hpp"
#include "idm.h"
#include "json.hpp"

namespace imj {

namespace {

constexpr double kEps = 1e-9;

int64_t floorDiv(int64_t a, int64_t b)
{
  int64_t q = a / b;
  if ((a ^ b) < 0 && a % b != 0) {
    --q;
  }
  return q;
}

int64_t ceilDiv(int64_t a, int64_t b)
{
  return -floorDiv(-a, b);
}

std::string toUpperAscii(std::string s)
{
  for (char& ch : s) {
    ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
  }
  return s;
}

std::vector<std::string> splitTokens(const std::string& input)
{
  std::vector<std::string> tokens;
  std::string cur;
  for (char ch : input) {
    if (ch == ',' || std::isspace(static_cast<unsigned char>(ch))) {
      if (!cur.empty()) {
        tokens.push_back(cur);
        cur.clear();
      }
    } else {
      cur.push_back(ch);
    }
  }
  if (!cur.empty()) {
    tokens.push_back(cur);
  }
  return tokens;
}

std::string extractStringConfig(const std::map<std::string, std::any>& config_map, const std::string& key,
                                const std::string& default_value)
{
  auto it = config_map.find(key);
  if (it == config_map.end() || !it->second.has_value()) {
    return default_value;
  }
  if (const auto* s = std::any_cast<std::string>(&it->second)) {
    return s->empty() ? default_value : *s;
  }
  if (const auto* cs = std::any_cast<const char*>(&it->second)) {
    return (*cs != nullptr && **cs != '\0') ? std::string(*cs) : default_value;
  }
  return default_value;
}

int extractIntConfig(const std::map<std::string, std::any>& config_map, const std::string& key, int default_value)
{
  auto it = config_map.find(key);
  if (it == config_map.end() || !it->second.has_value()) {
    return default_value;
  }
  if (const auto* v = std::any_cast<int>(&it->second)) {
    return *v;
  }
  if (const auto* v64 = std::any_cast<int64_t>(&it->second)) {
    return static_cast<int>(*v64);
  }
  if (const auto* vd = std::any_cast<double>(&it->second)) {
    return static_cast<int>(*vd);
  }
  if (const auto* s = std::any_cast<std::string>(&it->second)) {
    try {
      return std::stoi(*s);
    } catch (...) {
      return default_value;
    }
  }
  return default_value;
}

double extractDoubleConfig(const std::map<std::string, std::any>& config_map, const std::string& key, double default_value)
{
  auto it = config_map.find(key);
  if (it == config_map.end() || !it->second.has_value()) {
    return default_value;
  }
  if (const auto* vd = std::any_cast<double>(&it->second)) {
    return *vd;
  }
  if (const auto* vf = std::any_cast<float>(&it->second)) {
    return static_cast<double>(*vf);
  }
  if (const auto* vi = std::any_cast<int>(&it->second)) {
    return static_cast<double>(*vi);
  }
  if (const auto* v64 = std::any_cast<int64_t>(&it->second)) {
    return static_cast<double>(*v64);
  }
  if (const auto* s = std::any_cast<std::string>(&it->second)) {
    try {
      return std::stod(*s);
    } catch (...) {
      return default_value;
    }
  }
  return default_value;
}

struct DiodeCandidateMaster
{
  idb::IdbCellMaster* master = nullptr;
  idb::IdbTerm* signal_term = nullptr;
  double diff_area = 0.0;
  bool name_hint = false;
  uint32_t width = 0;
};

struct LayerShapeBox
{
  int layer_idx = -1;
  idb::IdbRect rect;
};

struct RowOccupancy
{
  idb::IdbRow* row = nullptr;
  int32_t x0 = 0;
  int32_t y0 = 0;
  int32_t step_x = 0;
  int32_t height = 0;
  int32_t num_sites = 0;
  idb::IdbOrient orient = idb::IdbOrient::kN_R0;
  std::vector<bool> occupied;

  void blockRect(const idb::IdbRect& r)
  {
    if (step_x <= 0 || height <= 0 || num_sites <= 0) {
      return;
    }
    if (r.get_low_y() >= y0 + height || r.get_high_y() <= y0) {
      return;
    }
    if (r.get_high_x() <= r.get_low_x()) {
      return;
    }
    int64_t first_s = std::max<int64_t>(0, floorDiv(static_cast<int64_t>(r.get_low_x()) - x0, step_x));
    int64_t last_s = std::min<int64_t>(num_sites - 1, floorDiv(static_cast<int64_t>(r.get_high_x()) - x0 - 1, step_x));
    for (int64_t s = first_s; s <= last_s; ++s) {
      occupied[static_cast<size_t>(s)] = true;
    }
  }
};

struct RoutingLayerInfo
{
  idb::IdbLayerRouting* layer = nullptr;
  int order = -1;
  int32_t width = 1;
  int32_t half_width = 1;
  int32_t spacing = 1;
  int32_t pitch_x = 1;
  int32_t pitch_y = 1;
  bool is_horizontal = true;
  bool is_vertical = false;
};

struct ViaOption
{
  idb::IdbVia* tmpl = nullptr;
  idb::IdbLayer* bottom_layer = nullptr;
  idb::IdbLayer* top_layer = nullptr;
  std::vector<idb::IdbRect> bottom_rects;
  std::vector<idb::IdbRect> top_rects;
};

struct StubMove
{
  bool is_via = false;
  int layer_idx = 0;      // For wire: routing layer idx. For via: lower routing layer idx.
  int next_layer_idx = 0; // For via: upper or lower target layer idx.
  int32_t x1 = 0;
  int32_t y1 = 0;
  int32_t x2 = 0;
  int32_t y2 = 0;
};

struct RoutePlan
{
  bool found = false;
  std::vector<StubMove> moves;
};

idb::IdbRect transformMasterRect(const idb::IdbRect& master_rect, int32_t inst_x, int32_t inst_y, idb::IdbOrient orient,
                                 uint32_t master_w, uint32_t master_h)
{
  idb::IdbRect rect(master_rect);
  rect.moveByStep(inst_x, inst_y);
  idb::IdbCoordinate<int32_t> origin(inst_x, inst_y);
  idb::IdbOrientTransform tf(orient, &origin, static_cast<int32_t>(master_w), static_cast<int32_t>(master_h));
  tf.transformRect(&rect);
  return rect;
}

bool rectClearOfObstacle(const idb::IdbRect& a, const idb::IdbRect& b, int32_t spacing)
{
  // Violation occurs when horizontal gap < spacing AND vertical gap < spacing.
  const int64_t sp = std::max<int32_t>(0, spacing);
  const bool x_too_close = static_cast<int64_t>(a.get_low_x()) < static_cast<int64_t>(b.get_high_x()) + sp &&
                           static_cast<int64_t>(a.get_high_x()) + sp > static_cast<int64_t>(b.get_low_x());
  const bool y_too_close = static_cast<int64_t>(a.get_low_y()) < static_cast<int64_t>(b.get_high_y()) + sp &&
                           static_cast<int64_t>(a.get_high_y()) + sp > static_cast<int64_t>(b.get_low_y());
  return !(x_too_close && y_too_close);
}

bool pointInsideRect(int32_t x, int32_t y, const idb::IdbRect& r)
{
  return x >= r.get_low_x() && x <= r.get_high_x() && y >= r.get_low_y() && y <= r.get_high_y();
}

idb::IdbRect shiftedRect(const idb::IdbRect& local_rect, int32_t cx, int32_t cy)
{
  return idb::IdbRect(local_rect.get_low_x() + cx, local_rect.get_low_y() + cy, local_rect.get_high_x() + cx,
                      local_rect.get_high_y() + cy);
}

std::vector<DiodeCandidateMaster> discoverDiodes(idb::IdbLayout* layout, const std::string& explicit_spec)
{
  std::vector<DiodeCandidateMaster> result;
  if (!layout || !layout->get_cell_master_list()) {
    return result;
  }
  auto* master_list = layout->get_cell_master_list();

  auto inspectMaster = [&](idb::IdbCellMaster* master, bool require_core_nonfiller) -> std::optional<DiodeCandidateMaster> {
    if (!master) {
      return std::nullopt;
    }
    if (require_core_nonfiller && (!master->is_core() || master->is_core_filler() || master->is_endcap())) {
      return std::nullopt;
    }
    idb::IdbTerm* best_term = nullptr;
    double best_diff = 0.0;
    for (idb::IdbTerm* term : master->get_term_list()) {
      if (!term || term->is_special_net() || term->is_power() || term->is_ground()) {
        continue;
      }
      if (term->get_type() == idb::IdbConnectType::kPower || term->get_type() == idb::IdbConnectType::kGround ||
          term->get_type() == idb::IdbConnectType::kClock) {
        continue;
      }
      if (term->get_port_list().empty()) {
        continue;
      }
      double gate = term->has_antenna_gate_area() ? term->get_antenna_gate_area() : 0.0;
      for (const auto& [lname, area] : term->get_antenna_gate_area_by_layer()) {
        (void) lname;
        gate += std::max(0.0, area);
      }
      if (gate > kEps) {
        return std::nullopt;  // Has gate oxide area; not a pure protection diode.
      }
      double diff = term->has_antenna_diff_area() ? std::max(0.0, term->get_antenna_diff_area()) : 0.0;
      for (const auto& [lname, area] : term->get_antenna_diff_area_by_layer()) {
        (void) lname;
        diff += std::max(0.0, area);
      }
      bool has_routing_shape = false;
      for (idb::IdbPort* port : term->get_port_list()) {
        if (!port) {
          continue;
        }
        for (idb::IdbLayerShape* sh : port->get_layer_shape()) {
          if (sh && sh->get_layer() && sh->get_layer()->is_routing() && !sh->get_rect_list().empty()) {
            has_routing_shape = true;
            break;
          }
        }
        if (has_routing_shape) {
          break;
        }
      }
      if (has_routing_shape && diff > best_diff + kEps) {
        best_diff = diff;
        best_term = term;
      }
    }
    if (!best_term) {
      return std::nullopt;
    }
    std::string upper = toUpperAscii(master->get_name());
    bool hint = (upper.find("ANT") != std::string::npos) || (upper.find("DIODE") != std::string::npos);
    return DiodeCandidateMaster{master, best_term, best_diff, hint, master->get_width()};
  };

  if (!explicit_spec.empty()) {
    for (const std::string& name : splitTokens(explicit_spec)) {
      if (auto cand = inspectMaster(master_list->find_cell_master(name), false)) {
        result.push_back(*cand);
      } else {
        MJLOG.warn(Loc::current(), "Requested diode master is missing or invalid: ", name);
      }
    }
    return result;
  }

  for (idb::IdbCellMaster* master : master_list->get_cell_master()) {
    if (auto cand = inspectMaster(master, true)) {
      result.push_back(*cand);
    }
  }

  std::sort(result.begin(), result.end(), [](const DiodeCandidateMaster& a, const DiodeCandidateMaster& b) {
    if (a.name_hint != b.name_hint) {
      return a.name_hint > b.name_hint;
    }
    if (a.width != b.width) {
      return a.width < b.width;
    }
    if (std::abs(a.diff_area - b.diff_area) > kEps) {
      return a.diff_area > b.diff_area;
    }
    return a.master->get_name() < b.master->get_name();
  });
  return result;
}

class LocalRouterContext
{
 public:
  bool init(idb::IdbDesign* design, int max_layers_to_use)
  {
    _design = design;
    if (!_design || !_design->get_layout() || !_design->get_layout()->get_layers()) {
      return false;
    }
    auto* idb_layers = _design->get_layout()->get_layers();
    std::vector<idb::IdbLayerRouting*> rlayers;
    for (idb::IdbLayer* l : idb_layers->get_routing_layers()) {
      if (auto* r = dynamic_cast<idb::IdbLayerRouting*>(l)) {
        rlayers.push_back(r);
      }
    }
    std::sort(rlayers.begin(), rlayers.end(),
              [](idb::IdbLayerRouting* a, idb::IdbLayerRouting* b) { return a->get_order() < b->get_order(); });
    if (rlayers.empty()) {
      return false;
    }

    int use_count = std::min<int>(static_cast<int>(rlayers.size()), std::max(2, max_layers_to_use));
    _layers.clear();
    _order_to_idx.clear();
    for (int i = 0; i < use_count; ++i) {
      idb::IdbLayerRouting* r = rlayers[i];
      RoutingLayerInfo info;
      info.layer = r;
      info.order = static_cast<int>(r->get_order());
      int32_t w = r->get_width() > 0 ? r->get_width() : r->get_min_width();
      info.width = std::max<int32_t>(1, w);
      info.half_width = info.width / 2;
      int32_t sp = r->get_spacing(info.width);
      if (sp <= 0 && r->get_spacing_list()) {
        sp = r->get_spacing_list()->get_spacing(info.width);
      }
      info.spacing = sp > 0 ? sp : info.width;
      int32_t px = r->get_pitch_x();
      int32_t py = r->get_pitch_y();
      int32_t def_pitch = info.width + info.spacing;
      info.pitch_x = px > 0 ? px : (py > 0 ? py : def_pitch);
      info.pitch_y = py > 0 ? py : (px > 0 ? px : def_pitch);
      info.is_horizontal = r->is_horizontal();
      info.is_vertical = r->is_vertical();
      if (!info.is_horizontal && !info.is_vertical) {
        info.is_horizontal = (i % 2 == 0);
        info.is_vertical = !info.is_horizontal;
      }
      _order_to_idx[info.order] = i;
      _layers.push_back(info);
    }

    buildVias();
    rebuildStaticObstacles();
    return true;
  }

  int layerIdxByOrder(int order) const
  {
    auto it = _order_to_idx.find(order);
    return it == _order_to_idx.end() ? -1 : it->second;
  }

  const std::vector<RoutingLayerInfo>& layers() const { return _layers; }

  const ViaOption* viaBetween(int layer_a, int layer_b) const
  {
    int lo = std::min(layer_a, layer_b);
    if (lo < 0 || lo + 1 >= static_cast<int>(_layers.size())) {
      return nullptr;
    }
    auto it = _vias_by_low_idx.find(lo);
    if (it == _vias_by_low_idx.end() || !it->second.tmpl) {
      return nullptr;
    }
    return &it->second;
  }

  void appendRoutedRect(int layer_idx, idb::IdbNet* owner, const idb::IdbRect& r)
  {
    if (layer_idx < 0 || layer_idx >= static_cast<int>(_layers.size())) {
      return;
    }
    _net_shapes[layer_idx].push_back({owner, r});
  }

  void appendDiodeObstacle(int layer_idx, const idb::IdbRect& r)
  {
    if (layer_idx < 0 || layer_idx >= static_cast<int>(_layers.size())) {
      return;
    }
    _fixed_obstacles[layer_idx].push_back(r);
  }

  RoutePlan routeStub(idb::IdbNet* target_net, const std::vector<LayerShapeBox>& start_shapes,
                      const std::vector<LayerShapeBox>& goal_shapes,
                      const std::vector<LayerShapeBox>& candidate_temp_obstacles, int32_t search_radius_dbu) const
  {
    RoutePlan failed;
    if (start_shapes.empty() || goal_shapes.empty()) {
      return failed;
    }

    // Center of primary start shape.
    const idb::IdbRect& srect0 = start_shapes[0].rect;
    const int32_t sx0 = (srect0.get_low_x() + srect0.get_high_x()) / 2;
    const int32_t sy0 = (srect0.get_low_y() + srect0.get_high_y()) / 2;

    // Sort goal shapes by Manhattan distance from candidate diode pin center; keep closest subset.
    std::vector<std::pair<int64_t, LayerShapeBox>> ranked_goals;
    ranked_goals.reserve(goal_shapes.size());
    for (const LayerShapeBox& gs : goal_shapes) {
      if (gs.layer_idx < 0 || gs.layer_idx >= static_cast<int>(_layers.size())) {
        continue;
      }
      int32_t cx = std::clamp(sx0, gs.rect.get_low_x(), gs.rect.get_high_x());
      int32_t cy = std::clamp(sy0, gs.rect.get_low_y(), gs.rect.get_high_y());
      int64_t dist = static_cast<int64_t>(std::abs(cx - sx0)) + std::abs(cy - sy0);
      ranked_goals.push_back({dist, gs});
    }
    if (ranked_goals.empty()) {
      return failed;
    }
    std::sort(ranked_goals.begin(), ranked_goals.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    if (ranked_goals.size() > 6) {
      ranked_goals.resize(6);
    }

    int32_t min_px = _layers[0].pitch_x;
    int32_t min_py = _layers[0].pitch_y;
    for (const auto& info : _layers) {
      min_px = std::min(min_px, info.pitch_x);
      min_py = std::min(min_py, info.pitch_y);
    }
    const int32_t pad = std::max<int32_t>(std::max(min_px, min_py) * 4, search_radius_dbu / 4);

    int32_t win_lx = sx0;
    int32_t win_hx = sx0;
    int32_t win_ly = sy0;
    int32_t win_hy = sy0;
    for (const LayerShapeBox& ss : start_shapes) {
      win_lx = std::min(win_lx, ss.rect.get_low_x());
      win_hx = std::max(win_hx, ss.rect.get_high_x());
      win_ly = std::min(win_ly, ss.rect.get_low_y());
      win_hy = std::max(win_hy, ss.rect.get_high_y());
    }
    for (const auto& [d, gs] : ranked_goals) {
      (void) d;
      int32_t cx = std::clamp(sx0, gs.rect.get_low_x(), gs.rect.get_high_x());
      int32_t cy = std::clamp(sy0, gs.rect.get_low_y(), gs.rect.get_high_y());
      win_lx = std::min(win_lx, cx);
      win_hx = std::max(win_hx, cx);
      win_ly = std::min(win_ly, cy);
      win_hy = std::max(win_hy, cy);
    }
    win_lx -= pad;
    win_hx += pad;
    win_ly -= pad;
    win_hy += pad;

    const idb::IdbRect win_rect(win_lx, win_ly, win_hx, win_hy);
    const int num_layers = static_cast<int>(_layers.size());

    // Build per-layer filtered obstacle lists inside window (+ clearance margin).
    std::vector<std::vector<idb::IdbRect>> obs_by_layer(num_layers);
    auto collectObstacle = [&](int lidx, const idb::IdbRect& r) {
      int32_t margin = _layers[lidx].half_width + _layers[lidx].spacing + 2 * _layers[lidx].width;
      if (r.get_high_x() + margin < win_lx || r.get_low_x() - margin > win_hx ||
          r.get_high_y() + margin < win_ly || r.get_low_y() - margin > win_hy) {
        return;
      }
      obs_by_layer[lidx].push_back(r);
    };
    for (int l = 0; l < num_layers; ++l) {
      for (const idb::IdbRect& r : _fixed_obstacles[l]) {
        collectObstacle(l, r);
      }
      for (const auto& [owner, r] : _net_shapes[l]) {
        if (owner != target_net) {
          collectObstacle(l, r);
        }
      }
    }
    for (const LayerShapeBox& ts : candidate_temp_obstacles) {
      if (ts.layer_idx >= 0 && ts.layer_idx < num_layers) {
        collectObstacle(ts.layer_idx, ts.rect);
      }
    }

    // Build sparse non-uniform X and Y coordinate sets.
    std::vector<int32_t> must_xs;
    std::vector<int32_t> must_ys;
    std::vector<int32_t> extra_xs;
    std::vector<int32_t> extra_ys;

    for (const LayerShapeBox& ss : start_shapes) {
      must_xs.push_back((ss.rect.get_low_x() + ss.rect.get_high_x()) / 2);
      must_ys.push_back((ss.rect.get_low_y() + ss.rect.get_high_y()) / 2);
    }
    for (const auto& [d, gs] : ranked_goals) {
      (void) d;
      must_xs.push_back(std::clamp(sx0, gs.rect.get_low_x(), gs.rect.get_high_x()));
      must_ys.push_back(std::clamp(sy0, gs.rect.get_low_y(), gs.rect.get_high_y()));
      must_xs.push_back((gs.rect.get_low_x() + gs.rect.get_high_x()) / 2);
      must_ys.push_back((gs.rect.get_low_y() + gs.rect.get_high_y()) / 2);
    }

    // Add track grid lines + obstacle bypass lines around obstacles in window.
    for (int l = 0; l < num_layers; ++l) {
      const auto& info = _layers[l];
      for (idb::IdbTrackGrid* tg : info.layer->get_track_grid_list()) {
        if (!tg || !tg->get_track() || tg->get_track()->get_pitch() == 0) {
          continue;
        }
        int32_t start = static_cast<int32_t>(tg->get_track()->get_start());
        int32_t pitch = static_cast<int32_t>(tg->get_track()->get_pitch());
        int32_t ntracks = static_cast<int32_t>(tg->get_track_num());
        int32_t max_coord = start + std::max<int32_t>(0, ntracks - 1) * pitch;
        if (tg->get_track()->is_track_direction_x()) {
          int64_t k0 = std::max<int64_t>(0, ceilDiv(static_cast<int64_t>(win_lx) - start, pitch));
          int64_t k1 = std::min<int64_t>(ntracks > 0 ? ntracks - 1 : floorDiv(static_cast<int64_t>(win_hx) - start, pitch),
                                         floorDiv(static_cast<int64_t>(win_hx) - start, pitch));
          int64_t step = std::max<int64_t>(1, (k1 - k0 + 1) / 24);
          for (int64_t k = k0; k <= k1; k += step) {
            int32_t c = static_cast<int32_t>(start + k * pitch);
            if (ntracks == 0 || c <= max_coord) {
              extra_xs.push_back(c);
            }
          }
        } else if (tg->get_track()->is_track_direction_y()) {
          int64_t k0 = std::max<int64_t>(0, ceilDiv(static_cast<int64_t>(win_ly) - start, pitch));
          int64_t k1 = std::min<int64_t>(ntracks > 0 ? ntracks - 1 : floorDiv(static_cast<int64_t>(win_hy) - start, pitch),
                                         floorDiv(static_cast<int64_t>(win_hy) - start, pitch));
          int64_t step = std::max<int64_t>(1, (k1 - k0 + 1) / 24);
          for (int64_t k = k0; k <= k1; k += step) {
            int32_t c = static_cast<int32_t>(start + k * pitch);
            if (ntracks == 0 || c <= max_coord) {
              extra_ys.push_back(c);
            }
          }
        }
      }

      // Local pitch rays around start/goal points for designs without track grids or fine snapping.
      int32_t gx0 = must_xs.size() > 1 ? must_xs[1] : sx0;
      int32_t gy0 = must_ys.size() > 1 ? must_ys[1] : sy0;
      for (int k = -4; k <= 4; ++k) {
        for (int32_t base_x : {sx0, gx0}) {
          int32_t cx = base_x + k * info.pitch_x;
          if (cx >= win_lx && cx <= win_hx) {
            extra_xs.push_back(cx);
          }
        }
        for (int32_t base_y : {sy0, gy0}) {
          int32_t cy = base_y + k * info.pitch_y;
          if (cy >= win_ly && cy <= win_hy) {
            extra_ys.push_back(cy);
          }
        }
      }

      // Obstacle escape coordinates (exact minimum clearance offset around local obstacles).
      int32_t clear = info.half_width + info.spacing + 1;
      int obs_added = 0;
      for (const idb::IdbRect& ob : obs_by_layer[l]) {
        if (rectClearOfObstacle(ob, win_rect, 1)) {
          continue;
        }
        if (obs_added >= 24) {
          break;
        }
        ++obs_added;
        if (ob.get_low_x() - clear >= win_lx && ob.get_low_x() - clear <= win_hx) {
          extra_xs.push_back(ob.get_low_x() - clear);
        }
        if (ob.get_high_x() + clear >= win_lx && ob.get_high_x() + clear <= win_hx) {
          extra_xs.push_back(ob.get_high_x() + clear);
        }
        if (ob.get_low_y() - clear >= win_ly && ob.get_low_y() - clear <= win_hy) {
          extra_ys.push_back(ob.get_low_y() - clear);
        }
        if (ob.get_high_y() + clear >= win_ly && ob.get_high_y() + clear <= win_hy) {
          extra_ys.push_back(ob.get_high_y() + clear);
        }
      }
    }

    auto compactAxis = [](std::vector<int32_t> must, std::vector<int32_t> extra, int32_t lo, int32_t hi,
                          size_t max_count) -> std::vector<int32_t> {
      auto dedup = [&](std::vector<int32_t>& v) {
        v.erase(std::remove_if(v.begin(), v.end(), [&](int32_t c) { return c < lo || c > hi; }), v.end());
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
      };
      dedup(must);
      dedup(extra);
      if (must.size() >= max_count) {
        return must;
      }
      std::vector<int32_t> merged = must;
      size_t budget = max_count - must.size();
      if (extra.size() <= budget) {
        merged.insert(merged.end(), extra.begin(), extra.end());
      } else {
        for (size_t i = 0; i < budget; ++i) {
          size_t idx = (i * extra.size()) / budget;
          merged.push_back(extra[idx]);
        }
      }
      dedup(merged);
      return merged;
    };

    const std::vector<int32_t> xs = compactAxis(must_xs, extra_xs, win_lx, win_hx, 64);
    const std::vector<int32_t> ys = compactAxis(must_ys, extra_ys, win_ly, win_hy, 64);
    const int nx = static_cast<int>(xs.size());
    const int ny = static_cast<int>(ys.size());
    if (nx == 0 || ny == 0) {
      return failed;
    }

    auto segValid = [&](int lidx, int32_t x1, int32_t y1, int32_t x2, int32_t y2) -> bool {
      int32_t hw = _layers[lidx].half_width;
      int32_t sp = _layers[lidx].spacing;
      idb::IdbRect wire_box(std::min(x1, x2) - hw, std::min(y1, y2) - hw, std::max(x1, x2) + hw,
                            std::max(y1, y2) + hw);
      for (const idb::IdbRect& ob : obs_by_layer[lidx]) {
        if (!rectClearOfObstacle(wire_box, ob, sp)) {
          return false;
        }
      }
      return true;
    };

    auto viaValid = [&](int low_idx, int32_t vx, int32_t vy) -> bool {
      const ViaOption* vo = viaBetween(low_idx, low_idx + 1);
      if (!vo) {
        return false;
      }
      for (const idb::IdbRect& br : vo->bottom_rects) {
        idb::IdbRect placed = shiftedRect(br, vx, vy);
        for (const idb::IdbRect& ob : obs_by_layer[low_idx]) {
          if (!rectClearOfObstacle(placed, ob, _layers[low_idx].spacing)) {
            return false;
          }
        }
      }
      for (const idb::IdbRect& tr : vo->top_rects) {
        idb::IdbRect placed = shiftedRect(tr, vx, vy);
        for (const idb::IdbRect& ob : obs_by_layer[low_idx + 1]) {
          if (!rectClearOfObstacle(placed, ob, _layers[low_idx + 1].spacing)) {
            return false;
          }
        }
      }
      return true;
    };

    auto isGoalCell = [&](int ix, int iy, int l) -> bool {
      int32_t px = xs[ix];
      int32_t py = ys[iy];
      for (const auto& [d, gs] : ranked_goals) {
        (void) d;
        if (gs.layer_idx == l && pointInsideRect(px, py, gs.rect)) {
          return true;
        }
      }
      return false;
    };

    auto heuristic = [&](int ix, int iy) -> int64_t {
      int64_t best = std::numeric_limits<int64_t>::max() / 4;
      int32_t px = xs[ix];
      int32_t py = ys[iy];
      for (const auto& [d, gs] : ranked_goals) {
        (void) d;
        int32_t cx = std::clamp(px, gs.rect.get_low_x(), gs.rect.get_high_x());
        int32_t cy = std::clamp(py, gs.rect.get_low_y(), gs.rect.get_high_y());
        best = std::min<int64_t>(best, static_cast<int64_t>(std::abs(px - cx)) + std::abs(py - cy));
      }
      return best;
    };

    const int total = num_layers * nx * ny;
    auto nodeId = [&](int ix, int iy, int l) -> int { return (l * ny + iy) * nx + ix; };

    std::vector<int64_t> g_cost(total, std::numeric_limits<int64_t>::max() / 4);
    std::vector<int> parent(total, -1);
    std::vector<bool> closed(total, false);

    struct PQItem
    {
      int64_t f;
      int64_t g;
      int id;
      bool operator>(const PQItem& o) const { return f > o.f; }
    };
    std::priority_queue<PQItem, std::vector<PQItem>, std::greater<PQItem>> pq;

    for (const LayerShapeBox& ss : start_shapes) {
      if (ss.layer_idx < 0 || ss.layer_idx >= num_layers) {
        continue;
      }
      for (int ix = 0; ix < nx; ++ix) {
        if (xs[ix] < ss.rect.get_low_x() || xs[ix] > ss.rect.get_high_x()) {
          continue;
        }
        for (int iy = 0; iy < ny; ++iy) {
          if (ys[iy] < ss.rect.get_low_y() || ys[iy] > ss.rect.get_high_y()) {
            continue;
          }
          int sid = nodeId(ix, iy, ss.layer_idx);
          if (g_cost[sid] > 0) {
            g_cost[sid] = 0;
            pq.push({heuristic(ix, iy), 0, sid});
          }
        }
      }
    }

    const int64_t via_penalty = std::max<int64_t>(6 * std::max(min_px, min_py), 4000);
    int goal_id = -1;

    while (!pq.empty()) {
      PQItem top = pq.top();
      pq.pop();
      if (closed[top.id] || top.g > g_cost[top.id]) {
        continue;
      }
      closed[top.id] = true;

      int rem = top.id;
      int ix = rem % nx;
      rem /= nx;
      int iy = rem % ny;
      int l = rem / ny;

      if (isGoalCell(ix, iy, l)) {
        goal_id = top.id;
        break;
      }

      const int32_t x = xs[ix];
      const int32_t y = ys[iy];
      const bool pref_h = _layers[l].is_horizontal;
      const bool pref_v = _layers[l].is_vertical;

      auto relax = [&](int nx_i, int ny_i, int nl, int64_t step_cost) {
        int nid = nodeId(nx_i, ny_i, nl);
        if (closed[nid]) {
          return;
        }
        int64_t ng = top.g + step_cost;
        if (ng < g_cost[nid]) {
          g_cost[nid] = ng;
          parent[nid] = top.id;
          pq.push({ng + heuristic(nx_i, ny_i), ng, nid});
        }
      };

      if (ix > 0 && segValid(l, x, y, xs[ix - 1], y)) {
        int64_t dx = x - xs[ix - 1];
        relax(ix - 1, iy, l, pref_h ? dx : (dx * 6) / 5);
      }
      if (ix + 1 < nx && segValid(l, x, y, xs[ix + 1], y)) {
        int64_t dx = xs[ix + 1] - x;
        relax(ix + 1, iy, l, pref_h ? dx : (dx * 6) / 5);
      }
      if (iy > 0 && segValid(l, x, y, x, ys[iy - 1])) {
        int64_t dy = y - ys[iy - 1];
        relax(ix, iy - 1, l, pref_v ? dy : (dy * 6) / 5);
      }
      if (iy + 1 < ny && segValid(l, x, y, x, ys[iy + 1])) {
        int64_t dy = ys[iy + 1] - y;
        relax(ix, iy + 1, l, pref_v ? dy : (dy * 6) / 5);
      }
      if (l + 1 < num_layers && viaValid(l, x, y)) {
        relax(ix, iy, l + 1, via_penalty);
      }
      if (l > 0 && viaValid(l - 1, x, y)) {
        relax(ix, iy, l - 1, via_penalty);
      }
    }

    if (goal_id < 0) {
      return failed;
    }

    std::vector<int> chain;
    for (int cur = goal_id; cur >= 0; cur = parent[cur]) {
      chain.push_back(cur);
    }
    std::reverse(chain.begin(), chain.end());

    struct GridPt
    {
      int32_t x;
      int32_t y;
      int l;
    };
    std::vector<GridPt> pts;
    pts.reserve(chain.size());
    for (int id : chain) {
      int rem = id;
      int ix = rem % nx;
      rem /= nx;
      int iy = rem % ny;
      int l = rem / ny;
      pts.push_back({xs[ix], ys[iy], l});
    }

    RoutePlan plan;
    plan.found = true;
    for (size_t i = 0; i + 1 < pts.size(); ++i) {
      const GridPt& a = pts[i];
      const GridPt& b = pts[i + 1];
      if (a.l != b.l) {
        StubMove vm;
        vm.is_via = true;
        vm.layer_idx = std::min(a.l, b.l);
        vm.next_layer_idx = std::max(a.l, b.l);
        vm.x1 = a.x;
        vm.y1 = a.y;
        vm.x2 = a.x;
        vm.y2 = a.y;
        plan.moves.push_back(vm);
        continue;
      }
      if (a.x == b.x && a.y == b.y) {
        continue;
      }
      if (!plan.moves.empty() && !plan.moves.back().is_via && plan.moves.back().layer_idx == a.l) {
        StubMove& prev = plan.moves.back();
        if (prev.x1 == prev.x2 && prev.x2 == a.x && a.x == b.x) {
          prev.y2 = b.y;
          continue;
        }
        if (prev.y1 == prev.y2 && prev.y2 == a.y && a.y == b.y) {
          prev.x2 = b.x;
          continue;
        }
      }
      StubMove wm;
      wm.is_via = false;
      wm.layer_idx = a.l;
      wm.next_layer_idx = a.l;
      wm.x1 = a.x;
      wm.y1 = a.y;
      wm.x2 = b.x;
      wm.y2 = b.y;
      plan.moves.push_back(wm);
    }
    return plan;
  }

 private:
  idb::IdbDesign* _design = nullptr;
  std::vector<RoutingLayerInfo> _layers;
  std::map<int, int> _order_to_idx;
  std::map<int, ViaOption> _vias_by_low_idx;
  std::vector<std::vector<idb::IdbRect>> _fixed_obstacles;
  std::vector<std::vector<std::pair<idb::IdbNet*, idb::IdbRect>>> _net_shapes;

  void buildVias()
  {
    _vias_by_low_idx.clear();
    std::map<int, std::tuple<bool, bool, int64_t>> best_score;

    auto scanViaList = [&](idb::IdbVias* vias) {
      if (!vias) {
        return;
      }
      for (idb::IdbVia* via : vias->get_via_list()) {
        if (!via || !via->get_instance()) {
          continue;
        }
        idb::IdbViaMaster* vm = via->get_instance();
        idb::IdbLayerShape* bs = vm->get_bottom_layer_shape();
        idb::IdbLayerShape* ts = vm->get_top_layer_shape();
        if (!bs || !ts || !bs->get_layer() || !ts->get_layer()) {
          continue;
        }
        int bi = layerIdxByOrder(static_cast<int>(bs->get_layer()->get_order()));
        int ti = layerIdxByOrder(static_cast<int>(ts->get_layer()->get_order()));
        if (bi < 0 || ti < 0 || std::abs(bi - ti) != 1) {
          continue;
        }
        int lo = std::min(bi, ti);
        int64_t area = 0;
        for (idb::IdbRect* r : bs->get_rect_list()) {
          if (r) {
            area += r->get_area();
          }
        }
        auto score = std::make_tuple(vm->is_default(), vm->isOneCut(), -area);
        if (!best_score.count(lo) || score > best_score[lo]) {
          best_score[lo] = score;
          ViaOption opt;
          opt.tmpl = via;
          opt.bottom_layer = (bi < ti) ? bs->get_layer() : ts->get_layer();
          opt.top_layer = (bi < ti) ? ts->get_layer() : bs->get_layer();
          for (idb::IdbRect* r : (bi < ti ? bs : ts)->get_rect_list()) {
            if (r) {
              opt.bottom_rects.push_back(*r);
            }
          }
          for (idb::IdbRect* r : (bi < ti ? ts : bs)->get_rect_list()) {
            if (r) {
              opt.top_rects.push_back(*r);
            }
          }
          _vias_by_low_idx[lo] = std::move(opt);
        }
      }
    };

    if (_design && _design->get_layout()) {
      scanViaList(_design->get_layout()->get_via_list());
    }
    if (_design) {
      scanViaList(_design->get_via_list());
    }
  }

  void rebuildStaticObstacles()
  {
    const int num_layers = static_cast<int>(_layers.size());
    _fixed_obstacles.assign(num_layers, {});
    _net_shapes.assign(num_layers, {});

    if (!_design) {
      return;
    }

    // 1. Placed instances' OBS shapes & IO/signal/unconnected pins.
    if (_design->get_instance_list()) {
      for (idb::IdbInstance* inst : _design->get_instance_list()->get_instance_list()) {
        if (!inst || inst->is_unplaced()) {
          continue;
        }
        for (idb::IdbLayerShape* obs : inst->get_obs_box_list()) {
          if (!obs || !obs->get_layer()) {
            continue;
          }
          int li = layerIdxByOrder(static_cast<int>(obs->get_layer()->get_order()));
          if (li < 0) {
            continue;
          }
          for (idb::IdbRect* r : obs->get_rect_list()) {
            if (r) {
              _fixed_obstacles[li].push_back(*r);
            }
          }
        }
        if (inst->get_pin_list()) {
          for (idb::IdbPin* pin : inst->get_pin_list()->get_pin_list()) {
            if (!pin) {
              continue;
            }
            idb::IdbNet* net = pin->get_net();
            for (idb::IdbLayerShape* sh : pin->get_port_box_list()) {
              if (!sh || !sh->get_layer()) {
                continue;
              }
              int li = layerIdxByOrder(static_cast<int>(sh->get_layer()->get_order()));
              if (li < 0) {
                continue;
              }
              for (idb::IdbRect* r : sh->get_rect_list()) {
                if (!r) {
                  continue;
                }
                if (net) {
                  _net_shapes[li].push_back({net, *r});
                } else {
                  _fixed_obstacles[li].push_back(*r);
                }
              }
            }
          }
        }
      }
    }

    // 2. IO pins.
    if (_design->get_io_pin_list()) {
      for (idb::IdbPin* pin : _design->get_io_pin_list()->get_pin_list()) {
        if (!pin) {
          continue;
        }
        idb::IdbNet* net = pin->get_net();
        for (idb::IdbLayerShape* sh : pin->get_port_box_list()) {
          if (!sh || !sh->get_layer()) {
            continue;
          }
          int li = layerIdxByOrder(static_cast<int>(sh->get_layer()->get_order()));
          if (li < 0) {
            continue;
          }
          for (idb::IdbRect* r : sh->get_rect_list()) {
            if (!r) {
              continue;
            }
            if (net) {
              _net_shapes[li].push_back({net, *r});
            } else {
              _fixed_obstacles[li].push_back(*r);
            }
          }
        }
      }
    }

    // 3. Regular nets' wires & vias.
    if (_design->get_net_list()) {
      for (idb::IdbNet* net : _design->get_net_list()->get_net_list()) {
        if (!net || !net->get_wire_list()) {
          continue;
        }
        for (idb::IdbRegularWire* wire : net->get_wire_list()->get_wire_list()) {
          if (!wire) {
            continue;
          }
          for (idb::IdbRegularWireSegment* seg : wire->get_segment_list()) {
            if (!seg) {
              continue;
            }
            if (seg->is_wire() || seg->is_rect()) {
              if (seg->get_layer()) {
                int li = layerIdxByOrder(static_cast<int>(seg->get_layer()->get_order()));
                if (li >= 0) {
                  _net_shapes[li].push_back({net, seg->get_segment_rect()});
                }
              }
            }
            if (seg->is_via()) {
              for (idb::IdbVia* via : seg->get_via_list()) {
                if (!via) {
                  continue;
                }
                for (idb::IdbLayerShape sh : {via->get_bottom_layer_shape(), via->get_top_layer_shape()}) {
                  if (!sh.get_layer()) {
                    continue;
                  }
                  int li = layerIdxByOrder(static_cast<int>(sh.get_layer()->get_order()));
                  if (li < 0) {
                    continue;
                  }
                  for (idb::IdbRect* r : sh.get_rect_list()) {
                    if (r) {
                      _net_shapes[li].push_back({net, *r});
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // 4. Special nets (VDD/VSS PG grids/rails).
    if (_design->get_special_net_list()) {
      for (idb::IdbSpecialNet* snet : _design->get_special_net_list()->get_net_list()) {
        if (!snet || !snet->get_wire_list()) {
          continue;
        }
        for (idb::IdbSpecialWire* wire : snet->get_wire_list()->get_wire_list()) {
          if (!wire) {
            continue;
          }
          for (idb::IdbSpecialWireSegment* seg : wire->get_segment_list()) {
            if (!seg) {
              continue;
            }
            if (seg->is_via() && seg->get_via()) {
              idb::IdbVia* via = seg->get_via();
              for (idb::IdbLayerShape sh : {via->get_bottom_layer_shape(), via->get_top_layer_shape()}) {
                if (!sh.get_layer()) {
                  continue;
                }
                int li = layerIdxByOrder(static_cast<int>(sh.get_layer()->get_order()));
                if (li < 0) {
                  continue;
                }
                for (idb::IdbRect* r : sh.get_rect_list()) {
                  if (r) {
                    _fixed_obstacles[li].push_back(*r);
                  }
                }
              }
            } else if (seg->is_rect() && seg->get_layer() && seg->get_delta_rect()) {
              int li = layerIdxByOrder(static_cast<int>(seg->get_layer()->get_order()));
              if (li >= 0) {
                idb::IdbRect rect(seg->get_delta_rect());
                if (seg->get_point_start()) {
                  rect.moveByStep(seg->get_point_start()->get_x(), seg->get_point_start()->get_y());
                }
                _fixed_obstacles[li].push_back(rect);
              }
            } else if (seg->get_layer() && seg->get_point_num() >= 2) {
              int li = layerIdxByOrder(static_cast<int>(seg->get_layer()->get_order()));
              if (li >= 0) {
                idb::IdbCoordinate<int32_t>* p1 = seg->get_point_start();
                idb::IdbCoordinate<int32_t>* p2 = seg->get_point_second();
                int32_t w = seg->get_route_width() > 0 ? seg->get_route_width() : _layers[li].width;
                int32_t hw = (w + 1) / 2;
                if (p1 && p2) {
                  _fixed_obstacles[li].emplace_back(std::min(p1->get_x(), p2->get_x()) - hw,
                                                    std::min(p1->get_y(), p2->get_y()) - hw,
                                                    std::max(p1->get_x(), p2->get_x()) + hw,
                                                    std::max(p1->get_y(), p2->get_y()) + hw);
                }
              }
            }
          }
        }
      }
    }

    // 5. Routing blockages.
    if (_design->get_blockage_list()) {
      for (idb::IdbBlockage* blk : _design->get_blockage_list()->get_blockage_list()) {
        auto* rblk = dynamic_cast<idb::IdbRoutingBlockage*>(blk);
        if (!rblk || !rblk->get_layer()) {
          continue;
        }
        int li = layerIdxByOrder(static_cast<int>(rblk->get_layer()->get_order()));
        if (li < 0) {
          continue;
        }
        for (idb::IdbRect* r : rblk->get_rect_list()) {
          if (r) {
            _fixed_obstacles[li].push_back(*r);
          }
        }
      }
    }
  }
};

std::vector<RowOccupancy> buildRowOccupancy(idb::IdbDesign* design)
{
  std::vector<RowOccupancy> rows;
  if (!design || !design->get_layout() || !design->get_layout()->get_rows()) {
    return rows;
  }
  for (idb::IdbRow* row : design->get_layout()->get_rows()->get_row_list()) {
    if (!row || !row->is_horizontal() || !row->get_site() || !row->get_original_coordinate()) {
      continue;
    }
    int32_t step_x = row->get_step_x() > 0 ? row->get_step_x() : row->get_site()->get_width();
    int32_t height = row->get_site()->get_height();
    int32_t num_sites = row->get_row_num_x();
    if (step_x <= 0 || height <= 0 || num_sites <= 0) {
      continue;
    }
    RowOccupancy ro;
    ro.row = row;
    ro.x0 = row->get_original_coordinate()->get_x();
    ro.y0 = row->get_original_coordinate()->get_y();
    ro.step_x = step_x;
    ro.height = height;
    ro.num_sites = num_sites;
    ro.orient = row->get_orient();
    ro.occupied.assign(static_cast<size_t>(num_sites), false);
    rows.push_back(std::move(ro));
  }

  if (design->get_instance_list()) {
    for (idb::IdbInstance* inst : design->get_instance_list()->get_instance_list()) {
      if (!inst || inst->is_unplaced() || !inst->has_placed() || !inst->get_bounding_box()) {
        continue;
      }
      const idb::IdbRect& bb = *inst->get_bounding_box();
      for (RowOccupancy& ro : rows) {
        ro.blockRect(bb);
      }
    }
  }

  if (design->get_blockage_list()) {
    for (idb::IdbBlockage* blk : design->get_blockage_list()->get_blockage_list()) {
      auto* pblk = dynamic_cast<idb::IdbPlacementBlockage*>(blk);
      if (!pblk) {
        continue;
      }
      for (idb::IdbRect* r : pblk->get_rect_list()) {
        if (!r) {
          continue;
        }
        for (RowOccupancy& ro : rows) {
          ro.blockRect(*r);
        }
      }
    }
  }

  return rows;
}

void connectDiodePGPins(idb::IdbDesign* design, idb::IdbInstance* inst, idb::IdbTerm* sig_term)
{
  if (!design || !inst || !inst->get_pin_list() || !design->get_special_net_list()) {
    return;
  }
  auto* snets = design->get_special_net_list();
  idb::IdbSpecialNet* default_pwr = nullptr;
  idb::IdbSpecialNet* default_gnd = nullptr;
  for (idb::IdbSpecialNet* sn : snets->get_net_list()) {
    if (!sn) {
      continue;
    }
    if (!default_pwr && sn->get_connect_type() == idb::IdbConnectType::kPower) {
      default_pwr = sn;
    }
    if (!default_gnd && sn->get_connect_type() == idb::IdbConnectType::kGround) {
      default_gnd = sn;
    }
  }
  for (idb::IdbPin* pin : inst->get_pin_list()->get_pin_list()) {
    if (!pin || !pin->get_term() || pin->get_term() == sig_term) {
      continue;
    }
    idb::IdbTerm* term = pin->get_term();
    idb::IdbSpecialNet* sn = snets->find_net(term->get_name());
    if (!sn && term->is_power()) {
      sn = default_pwr;
    } else if (!sn && term->is_ground()) {
      sn = default_gnd;
    }
    if (sn) {
      design->connectPinToSpecialNet(pin, sn);
    }
  }
}

void writeRepairOutputs(const ARSummary& summary, const std::string& report_path, const std::string& step_json_path)
{
  if (!report_path.empty()) {
    std::filesystem::path rpath(report_path);
    if (!rpath.parent_path().empty()) {
      std::filesystem::create_directories(rpath.parent_path());
    }
    std::ofstream stream(rpath);
    if (stream.is_open()) {
      stream << "========================================================================\n"
             << "                      Antenna Repair Report                           \n"
             << "========================================================================\n"
             << "Violations Before : " << summary.violations_before << "\n"
             << "Violations After  : " << summary.violations_after << "\n"
             << "Diodes Inserted   : " << summary.diodes_inserted << "\n"
             << "Nets Repaired     : " << summary.nets_repaired << "\n"
             << "Nets Unrepaired   : " << summary.nets_unrepaired << "\n"
             << "Iterations        : " << summary.iterations << "\n\n";
      if (!summary.placed_diodes.empty()) {
        stream << "Placed Diodes:\n"
               << "------------------------------------------------------------------------\n";
        for (const ARPlacedDiode& d : summary.placed_diodes) {
          stream << "  " << d.instance_name << " (" << d.master_name << ") -> net " << d.net_name
                 << " near " << d.target_pin << " at (" << d.x_dbu << ", " << d.y_dbu << ")"
                 << " [wires=" << d.wire_segments_added << ", vias=" << d.vias_added << "]\n";
        }
      }
    } else {
      MJLOG.error(Loc::current(), "Failed to write antenna repair report: ", report_path);
    }
  }

  if (!step_json_path.empty()) {
    std::filesystem::path jpath(step_json_path);
    if (!jpath.parent_path().empty()) {
      std::filesystem::create_directories(jpath.parent_path());
    }
    nlohmann::json payload = {
        {"schema_version", 1},
        {"step", "antennaRepair"},
        {"antenna_repair",
         {{"violations_before", summary.violations_before},
          {"violations_after", summary.violations_after},
          {"diodes_inserted", summary.diodes_inserted},
          {"nets_repaired", summary.nets_repaired},
          {"nets_unrepaired", summary.nets_unrepaired},
          {"iterations", summary.iterations},
          {"diode_masters_used", summary.diode_masters_used}}},
        {"placed_diodes", nlohmann::json::array()}};
    for (const ARPlacedDiode& d : summary.placed_diodes) {
      payload["placed_diodes"].push_back({
          {"instance", d.instance_name},
          {"master", d.master_name},
          {"net", d.net_name},
          {"target_pin", d.target_pin},
          {"x_dbu", d.x_dbu},
          {"y_dbu", d.y_dbu},
          {"wire_segments_added", d.wire_segments_added},
          {"vias_added", d.vias_added},
      });
    }
    std::ofstream jstream(jpath);
    if (!jstream.is_open()) {
      MJLOG.error(Loc::current(), "Failed to write antenna repair step JSON: ", step_json_path);
      return;
    }
    jstream << payload.dump(2) << '\n';
  }
}

}  // namespace

void AntennaRepairer::initInst()
{
  if (_ar_instance == nullptr) {
    _ar_instance = new AntennaRepairer();
  }
}

AntennaRepairer& AntennaRepairer::getInst()
{
  if (_ar_instance == nullptr) {
    MJLOG.error(Loc::current(), "The instance is not initialized!");
  }
  return *_ar_instance;
}

void AntennaRepairer::destroyInst()
{
  if (_ar_instance != nullptr) {
    delete _ar_instance;
    _ar_instance = nullptr;
  }
}

ARSummary AntennaRepairer::repair(std::map<std::string, std::any> config_map)
{
  Monitor monitor;
  MJLOG.info(Loc::current(), "Starting AntennaRepairer...");

  ARSummary summary;
  const std::string explicit_diode = extractStringConfig(config_map, "-diode", "");
  const int max_iter = std::max(1, extractIntConfig(config_map, "-max_iter", 3));
  const double search_radius_um = std::max(1.0, extractDoubleConfig(config_map, "-search_radius_um", 30.0));
  const std::string report_dir = extractStringConfig(config_map, "-report_dir", MJDM.getConfig().ar_temp_directory_path);
  const std::string report_path =
      extractStringConfig(config_map, "-report_path", report_dir.empty() ? "" : (report_dir + "/antenna_repair.rpt"));
  const std::string step_json_path = extractStringConfig(config_map, "-step_json", "");

  idb::IdbDesign* design = (dmInst != nullptr) ? dmInst->get_idb_design() : nullptr;
  if (!design || !design->get_layout()) {
    MJLOG.warn(Loc::current(), "No design loaded; skipping antenna repair.");
    writeRepairOutputs(summary, report_path, step_json_path);
    return summary;
  }

  idb::IdbLayout* layout = design->get_layout();
  int32_t dbu = (design->get_units() && design->get_units()->get_micron_dbu() > 0)
                    ? design->get_units()->get_micron_dbu()
                    : (layout->get_units() && layout->get_units()->get_micron_dbu() > 0 ? layout->get_units()->get_micron_dbu()
                                                                                       : 1000);
  const int32_t radius_dbu = static_cast<int32_t>(std::llround(search_radius_um * dbu));

  AntennaChecker::initInst();
  std::vector<ACViolation> violations = MJAC.checkAndCollect({}, false);
  summary.violations_before = static_cast<int>(violations.size());
  summary.violations_after = summary.violations_before;

  if (violations.empty()) {
    MJLOG.info(Loc::current(), "No antenna violations found; nothing to repair.");
    writeRepairOutputs(summary, report_path, step_json_path);
    AntennaChecker::destroyInst();
    MJLOG.info(Loc::current(), "Completed", monitor.getStatsInfo());
    return summary;
  }

  std::vector<DiodeCandidateMaster> diodes = discoverDiodes(layout, explicit_diode);
  if (diodes.empty()) {
    MJLOG.warn(Loc::current(), "No usable antenna diode master found in LEF library.");
    std::set<std::string> uniq_nets;
    for (const ACViolation& v : violations) {
      uniq_nets.insert(v.net_name);
    }
    summary.nets_unrepaired = static_cast<int>(uniq_nets.size());
    writeRepairOutputs(summary, report_path, step_json_path);
    AntennaChecker::destroyInst();
    return summary;
  }

  std::vector<RowOccupancy> rows = buildRowOccupancy(design);
  std::set<std::string> used_masters;
  std::set<std::string> repaired_nets;
  std::set<std::string> initial_viol_nets;
  for (const ACViolation& v : violations) {
    initial_viol_nets.insert(v.net_name);
  }

  for (int iter = 0; iter < max_iter && !violations.empty(); ++iter) {
    summary.iterations = iter + 1;

    // Group current violations by net and determine highest goal layer order needed.
    std::map<std::string, std::vector<const ACViolation*>> by_net;
    for (const ACViolation& v : violations) {
      by_net[v.net_name].push_back(&v);
    }

    int max_needed_idx = 2;  // Use at least M1..M3 for obstacle-avoiding bridge detours.
    for (const auto& [net_name, nvs] : by_net) {
      (void) nvs;
      idb::IdbNet* net = design->get_net_list() ? design->get_net_list()->find_net(net_name) : nullptr;
      if (!net || !net->get_instance_pin_list()) {
        continue;
      }
      for (idb::IdbPin* pin : net->get_instance_pin_list()->get_pin_list()) {
        if (!pin) {
          continue;
        }
        for (idb::IdbLayerShape* sh : pin->get_port_box_list()) {
          if (sh && sh->get_layer() && layout->get_layers()) {
            int k = 0;
            for (idb::IdbLayer* rl : layout->get_layers()->get_routing_layers()) {
              if (rl && rl->get_order() == sh->get_layer()->get_order()) {
                max_needed_idx = std::max(max_needed_idx, k + 1);
              }
              ++k;
            }
          }
        }
      }
    }

    LocalRouterContext router;
    if (!router.init(design, std::max(3, max_needed_idx))) {
      MJLOG.warn(Loc::current(), "Failed to initialize routing layer context for antenna repair.");
      break;
    }

    int inserted_this_iter = 0;

    for (const auto& [net_name, net_viols] : by_net) {
      idb::IdbNet* net = design->get_net_list() ? design->get_net_list()->find_net(net_name) : nullptr;
      if (!net) {
        continue;
      }

      // Identify anchor gate pins matching reported violations (fallback to all gate pins on net).
      std::vector<idb::IdbPin*> target_pins;
      std::unordered_set<idb::IdbPin*> seen_pins;
      for (const ACViolation* v : net_viols) {
        if (!net->get_instance_pin_list()) {
          continue;
        }
        for (idb::IdbPin* pin : net->get_instance_pin_list()->get_pin_list()) {
          if (!pin || seen_pins.count(pin)) {
            continue;
          }
          std::string inst_name = pin->get_instance() ? pin->get_instance()->get_name() : "";
          if (pin->get_pin_name() == v->pin_name && inst_name == v->instance_name && !pin->get_port_box_list().empty()) {
            seen_pins.insert(pin);
            target_pins.push_back(pin);
          }
        }
      }
      if (target_pins.empty() && net->get_instance_pin_list()) {
        for (idb::IdbPin* pin : net->get_instance_pin_list()->get_pin_list()) {
          if (!pin || !pin->get_term() || seen_pins.count(pin) || pin->get_port_box_list().empty()) {
            continue;
          }
          if (pin->get_term()->has_antenna_gate_area() || !pin->get_term()->get_antenna_gate_area_by_layer().empty()) {
            seen_pins.insert(pin);
            target_pins.push_back(pin);
          }
        }
      }
      if (target_pins.empty()) {
        continue;
      }

      // Collect goal shapes belonging to target_net (anchor pin port shapes + existing net wire/via shapes).
      std::vector<LayerShapeBox> goal_shapes;
      for (idb::IdbPin* gpin : target_pins) {
        for (idb::IdbLayerShape* sh : gpin->get_port_box_list()) {
          if (!sh || !sh->get_layer()) {
            continue;
          }
          int li = router.layerIdxByOrder(static_cast<int>(sh->get_layer()->get_order()));
          if (li < 0) {
            continue;
          }
          for (idb::IdbRect* r : sh->get_rect_list()) {
            if (r) {
              goal_shapes.push_back({li, *r});
            }
          }
        }
      }
      if (net->get_wire_list()) {
        for (idb::IdbRegularWire* wire : net->get_wire_list()->get_wire_list()) {
          if (!wire) {
            continue;
          }
          for (idb::IdbRegularWireSegment* seg : wire->get_segment_list()) {
            if (!seg) {
              continue;
            }
            if ((seg->is_wire() || seg->is_rect()) && seg->get_layer()) {
              int li = router.layerIdxByOrder(static_cast<int>(seg->get_layer()->get_order()));
              if (li >= 0) {
                goal_shapes.push_back({li, seg->get_segment_rect()});
              }
            }
            if (seg->is_via()) {
              for (idb::IdbVia* via : seg->get_via_list()) {
                if (!via) {
                  continue;
                }
                for (idb::IdbLayerShape sh : {via->get_bottom_layer_shape(), via->get_top_layer_shape()}) {
                  if (!sh.get_layer()) {
                    continue;
                  }
                  int li = router.layerIdxByOrder(static_cast<int>(sh.get_layer()->get_order()));
                  if (li < 0) {
                    continue;
                  }
                  for (idb::IdbRect* r : sh.get_rect_list()) {
                    if (r) {
                      goal_shapes.push_back({li, *r});
                    }
                  }
                }
              }
            }
          }
        }
      }
      if (goal_shapes.empty()) {
        continue;
      }

      // Compute anchor point around first target pin.
      idb::IdbPin* primary_pin = target_pins[0];
      idb::IdbCoordinate<int32_t>* avg = primary_pin->get_average_coordinate();
      const int32_t anchor_x = avg ? avg->get_x() : (goal_shapes[0].rect.get_low_x() + goal_shapes[0].rect.get_high_x()) / 2;
      const int32_t anchor_y = avg ? avg->get_y() : (goal_shapes[0].rect.get_low_y() + goal_shapes[0].rect.get_high_y()) / 2;

      bool placed_for_net = false;
      for (const DiodeCandidateMaster& dm : diodes) {
        if (placed_for_net) {
          break;
        }
        struct SiteCandidate
        {
          size_t row_idx;
          int32_t site_idx;
          int32_t sites_needed;
          int32_t inst_x;
          int32_t inst_y;
          int64_t dist;
          std::vector<LayerShapeBox> start_shapes;
          std::vector<LayerShapeBox> internal_obs;
        };
        std::vector<SiteCandidate> candidates;

        for (size_t ri = 0; ri < rows.size(); ++ri) {
          const RowOccupancy& ro = rows[ri];
          if (std::abs(ro.y0 + ro.height / 2 - anchor_y) > radius_dbu + ro.height) {
            continue;
          }
          if (dm.master->get_height() > static_cast<uint32_t>(ro.height) + 1) {
            continue;
          }
          int32_t sites_needed = static_cast<int32_t>(std::max<int64_t>(1, ceilDiv(dm.master->get_width(), ro.step_x)));
          if (sites_needed > ro.num_sites) {
            continue;
          }
          int32_t min_s = static_cast<int32_t>(
              std::max<int64_t>(0, floorDiv(static_cast<int64_t>(anchor_x - radius_dbu) - ro.x0, ro.step_x)));
          int32_t max_s = static_cast<int32_t>(
              std::min<int64_t>(ro.num_sites - sites_needed,
                                floorDiv(static_cast<int64_t>(anchor_x + radius_dbu) - ro.x0, ro.step_x)));

          for (int32_t s = min_s; s <= max_s; ++s) {
            bool free_run = true;
            for (int32_t k = 0; k < sites_needed; ++k) {
              if (ro.occupied[static_cast<size_t>(s + k)]) {
                free_run = false;
                break;
              }
            }
            if (!free_run) {
              continue;
            }
            int32_t inst_x = ro.x0 + s * ro.step_x;
            int32_t inst_y = ro.y0;
            int64_t approx_dist = static_cast<int64_t>(std::abs(inst_x - anchor_x)) + std::abs(inst_y - anchor_y);
            candidates.push_back({ri, s, sites_needed, inst_x, inst_y, approx_dist, {}, {}});
          }
        }

        std::sort(candidates.begin(), candidates.end(),
                  [](const SiteCandidate& a, const SiteCandidate& b) { return a.dist < b.dist; });
        if (candidates.size() > 24) {
          candidates.resize(24);
        }

        for (SiteCandidate& cand : candidates) {
          const RowOccupancy& ro = rows[cand.row_idx];
          for (idb::IdbPort* port : dm.signal_term->get_port_list()) {
            if (!port) {
              continue;
            }
            for (idb::IdbLayerShape* sh : port->get_layer_shape()) {
              if (!sh || !sh->get_layer()) {
                continue;
              }
              int li = router.layerIdxByOrder(static_cast<int>(sh->get_layer()->get_order()));
              if (li < 0) {
                continue;
              }
              for (idb::IdbRect* r : sh->get_rect_list()) {
                if (r) {
                  cand.start_shapes.push_back(
                      {li, transformMasterRect(*r, cand.inst_x, cand.inst_y, ro.orient, dm.master->get_width(),
                                               dm.master->get_height())});
                }
              }
            }
          }
          if (cand.start_shapes.empty()) {
            continue;
          }

          // Diode's non-signal pins (PG) & macro OBS become temporary routing obstacles.
          for (idb::IdbTerm* term : dm.master->get_term_list()) {
            if (!term || term == dm.signal_term) {
              continue;
            }
            for (idb::IdbPort* port : term->get_port_list()) {
              if (!port) {
                continue;
              }
              for (idb::IdbLayerShape* sh : port->get_layer_shape()) {
                if (!sh || !sh->get_layer()) {
                  continue;
                }
                int li = router.layerIdxByOrder(static_cast<int>(sh->get_layer()->get_order()));
                if (li < 0) {
                  continue;
                }
                for (idb::IdbRect* r : sh->get_rect_list()) {
                  if (r) {
                    cand.internal_obs.push_back(
                        {li, transformMasterRect(*r, cand.inst_x, cand.inst_y, ro.orient, dm.master->get_width(),
                                                 dm.master->get_height())});
                  }
                }
              }
            }
          }
          for (idb::IdbObs* obs : dm.master->get_obs_list()) {
            if (!obs) {
              continue;
            }
            for (idb::IdbObsLayer* ol : obs->get_obs_layer_list()) {
              if (!ol || !ol->get_shape() || !ol->get_shape()->get_layer()) {
                continue;
              }
              int li = router.layerIdxByOrder(static_cast<int>(ol->get_shape()->get_layer()->get_order()));
              if (li < 0) {
                continue;
              }
              for (idb::IdbRect* r : ol->get_shape()->get_rect_list()) {
                if (r) {
                  cand.internal_obs.push_back(
                      {li, transformMasterRect(*r, cand.inst_x, cand.inst_y, ro.orient, dm.master->get_width(),
                                               dm.master->get_height())});
                }
              }
            }
          }

          RoutePlan plan = router.routeStub(net, cand.start_shapes, goal_shapes, cand.internal_obs, radius_dbu);
          if (!plan.found) {
            continue;
          }

          // Commit placement & routing atomically.
          std::string inst_name = design->makeUniqueInstanceName("ant_diode");
          idb::IdbInstance* inst = design->createInstance(inst_name, dm.master->get_name(), idb::IdbInstanceType::kDist,
                                                          idb::IdbPlacementStatus::kPlaced, ro.orient, cand.inst_x, cand.inst_y,
                                                          idb::IdbCreatePolicy::kErrorIfExists);
          if (!inst) {
            continue;
          }
          idb::IdbPin* sig_pin = inst->get_pin(dm.signal_term->get_name());
          if (!sig_pin || !design->connectPinToNet(sig_pin, net)) {
            design->removeInstanceSafe(inst_name);
            continue;
          }
          connectDiodePGPins(design, inst, dm.signal_term);

          idb::IdbRegularWire* target_wire = nullptr;
          if (net->get_wire_list() && !net->get_wire_list()->get_wire_list().empty()) {
            target_wire = net->get_wire_list()->get_wire_list().front();
          } else if (net->get_wire_list()) {
            target_wire = net->get_wire_list()->add_wire();
            if (target_wire) {
              target_wire->set_wire_state(idb::IdbWiringStatement::kRouted);
            }
          }

          int wire_added = 0;
          int via_added = 0;
          if (target_wire) {
            bool first_seg = (target_wire->get_num() == 0);
            for (const StubMove& mv : plan.moves) {
              if (mv.is_via) {
                const ViaOption* vo = router.viaBetween(mv.layer_idx, mv.next_layer_idx);
                if (!vo || !vo->tmpl) {
                  continue;
                }
                idb::IdbRegularWireSegment* seg = target_wire->add_segment();
                if (!first_seg) {
                  seg->set_layer_as_new();
                }
                first_seg = false;
                seg->set_layer(vo->bottom_layer);
                seg->set_layer_name(vo->bottom_layer->get_name());
                seg->set_is_via(true);
                seg->add_point(mv.x1, mv.y1);
                idb::IdbVia* v = seg->copy_via(vo->tmpl);
                if (v) {
                  v->set_coordinate(mv.x1, mv.y1);
                  for (const idb::IdbRect& br : vo->bottom_rects) {
                    router.appendRoutedRect(mv.layer_idx, net, shiftedRect(br, mv.x1, mv.y1));
                  }
                  for (const idb::IdbRect& tr : vo->top_rects) {
                    router.appendRoutedRect(mv.next_layer_idx, net, shiftedRect(tr, mv.x1, mv.y1));
                  }
                }
                ++via_added;
              } else {
                const auto& linfo = router.layers()[mv.layer_idx];
                idb::IdbRegularWireSegment* seg = target_wire->add_segment();
                if (!first_seg) {
                  seg->set_layer_as_new();
                }
                first_seg = false;
                seg->set_layer(linfo.layer);
                seg->set_layer_name(linfo.layer->get_name());
                seg->add_point(mv.x1, mv.y1);
                seg->add_point(mv.x2, mv.y2);
                idb::IdbRect rbox(std::min(mv.x1, mv.x2) - linfo.half_width, std::min(mv.y1, mv.y2) - linfo.half_width,
                                  std::max(mv.x1, mv.x2) + linfo.half_width, std::max(mv.y1, mv.y2) + linfo.half_width);
                router.appendRoutedRect(mv.layer_idx, net, rbox);
                ++wire_added;
              }
            }
          }

          for (int32_t k = 0; k < cand.sites_needed; ++k) {
            rows[cand.row_idx].occupied[static_cast<size_t>(cand.site_idx + k)] = true;
          }
          for (const LayerShapeBox& sh : cand.start_shapes) {
            router.appendRoutedRect(sh.layer_idx, net, sh.rect);
          }
          for (const LayerShapeBox& obs : cand.internal_obs) {
            router.appendDiodeObstacle(obs.layer_idx, obs.rect);
          }

          ARPlacedDiode pd;
          pd.instance_name = inst_name;
          pd.master_name = dm.master->get_name();
          pd.net_name = net_name;
          pd.target_pin = (primary_pin->get_instance() ? primary_pin->get_instance()->get_name() + "/" : "") +
                          primary_pin->get_pin_name();
          pd.x_dbu = cand.inst_x;
          pd.y_dbu = cand.inst_y;
          pd.wire_segments_added = wire_added;
          pd.vias_added = via_added;
          summary.placed_diodes.push_back(std::move(pd));

          used_masters.insert(dm.master->get_name());
          repaired_nets.insert(net_name);
          ++summary.diodes_inserted;
          ++inserted_this_iter;
          placed_for_net = true;
          break;
        }
      }
    }

    if (inserted_this_iter == 0) {
      break;
    }
    violations = MJAC.checkAndCollect({}, false);
    summary.violations_after = static_cast<int>(violations.size());
  }

  std::set<std::string> remaining_nets;
  for (const ACViolation& v : violations) {
    remaining_nets.insert(v.net_name);
  }
  summary.nets_unrepaired = static_cast<int>(remaining_nets.size());
  int resolved_count = 0;
  for (const std::string& n : repaired_nets) {
    if (!remaining_nets.count(n)) {
      ++resolved_count;
    }
  }
  summary.nets_repaired = resolved_count;
  summary.diode_masters_used.assign(used_masters.begin(), used_masters.end());

  writeRepairOutputs(summary, report_path, step_json_path);
  AntennaChecker::destroyInst();
  MJLOG.info(Loc::current(), "Completed", monitor.getStatsInfo());
  return summary;
}

AntennaRepairer* AntennaRepairer::_ar_instance = nullptr;

}  // namespace imj
