#include "descartes_planner/ladder_graph_dag_search_lazy_collision.h"

namespace descartes_planner
{

DAGSearchLazyCollision::DAGSearchLazyCollision(const LadderGraph &graph)
  : graph_(graph)
{
  // On creating an object, let's allocate everything we need
  solution_.resize(graph.size());

  for (size_t i = 0; i < graph.size(); ++i)
  {
    const auto n_vertices = graph.rungSize(i);
    solution_[i].distance.resize(n_vertices);
    solution_[i].predecessor.resize(n_vertices);
    solution_[i].valid.resize(n_vertices, Checked::UNKNOWN); // all points start as valid
  }
}

double DAGSearchLazyCollision::run(const std::vector<planning_scene::PlanningScenePtr> &scenes,
                                   const std::vector<DAGSearchLazyCollision::size_type> &scene_max_indices)
{

  const auto* jmg = scenes.front()->getRobotModel()->getJointModelGroup("manipulator");

  size_t counter = 0;
  while (true)
  {
    ROS_INFO("Search iter %lu", counter++);
    // Run the search
    const bool has_possible_path = runOne();
    if (!has_possible_path)
    {
      return std::numeric_limits<double>::max();
    }

    // get the path's indices
    auto possible_path = shortestPath();

    moveit::core::RobotState state (scenes.front()->getRobotModel());
    state.setToDefaultValues();
    state.update();

    bool all_valid = true;
    size_t current_index = 0;
    for (std::size_t i = 0; i < scene_max_indices.size(); ++i)
    {
      const size_t start_rung = current_index;
      const size_t end_rung = start_rung + scene_max_indices[i];

      for (std::size_t rung = start_rung; rung < end_rung; ++rung)
      {
//        ROS_INFO_STREAM("Rung " << rung);
        const auto idx = possible_path[rung];

        if (valid(rung, idx) == Checked::UNKNOWN)
        {
          const auto* data = graph_.vertex(rung, idx);

          state.setJointGroupPositions(jmg, data);
          state.update();

          const planning_scene::PlanningScene& checker = *scenes[i];
          bool is_valid = !checker.isStateColliding(state, "manipulator", false);
          all_valid = all_valid && is_valid;
          valid(rung, idx) = is_valid ? Checked::FREE : Checked::COLLIDING;
          if (!is_valid)
          {
            ROS_WARN("Collision on index %lu", rung);
          }
        }
        else if (valid(rung, idx) == Checked::COLLIDING)
        {
          all_valid = false;
        } // otherwise we're good!
      }
      current_index = end_rung;
    }

    if (all_valid)
    {
      break;
    }
  }

  return *std::min_element(solution_.back().distance.begin(), solution_.back().distance.end());
}

bool DAGSearchLazyCollision::runOne()
{
  // Cost to the first rung should be set to zero
  std::fill(solution_.front().distance.begin(), solution_.front().distance.end(), 0.0);
  // Other rows initialize to zero
  for (size_type i = 1; i < solution_.size(); ++i)
  {
    std::fill(solution_[i].distance.begin(), solution_[i].distance.end(), std::numeric_limits<double>::max());
  }

  // Now we iterate over the graph in 'topological' order
  for (size_type rung = 0; rung < solution_.size() - 1; ++rung)
  {
    const auto n_vertices = graph_.rungSize(rung);
    const auto next_rung = rung + 1;
    // For each vertex in the out edge list
    for (size_t index = 0; index < n_vertices; ++index)
    {
      if (valid(rung, index) == Checked::COLLIDING) continue;

      const auto u_cost = distance(rung, index);
      const auto& edges = graph_.getEdges(rung)[index];
      // for each out edge
      for (const auto& edge : edges)
      {
        auto dv = u_cost + edge.cost; // new cost
        if (dv < distance(next_rung, edge.idx) && valid(next_rung, edge.idx) != Checked::COLLIDING)
        {
          distance(next_rung, edge.idx) = dv;
          predecessor(next_rung, edge.idx) = index; // the predecessor's rung is implied to be the current rung
        }
      }
    } // vertex for loop
  } // rung for loop

  return *std::min_element(solution_.back().distance.begin(), solution_.back().distance.end()) != std::numeric_limits<double>::max();
}



std::vector<DAGSearchLazyCollision::predecessor_t> DAGSearchLazyCollision::shortestPath() const
{
  auto min_it = std::min_element(solution_.back().distance.begin(), solution_.back().distance.end());
  auto min_idx = std::distance(solution_.back().distance.begin(), min_it);
  assert(min_idx >= 0);

  std::vector<predecessor_t> path (solution_.size());

  size_type current_rung = path.size() - 1;
  size_type current_index = min_idx;

  for (unsigned i = 0; i < path.size(); ++i)
  {
    auto count = path.size() - 1 - i;
    assert(current_rung == count);
    path[count] = current_index;
    current_index = predecessor(current_rung, current_index);
    current_rung -= 1;
  }

  return path;
}

} // namespace descartes_planner
