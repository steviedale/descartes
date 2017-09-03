#include "descartes_planner/graph_builder.h"
namespace descartes_planner
{
typedef boost::function<double(const double*, const double*)> CostFunction;
}
#include "descartes_planner/planning_graph_edge_policy.h"

namespace // anon namespace to hide utility functions
{

using PositionVector = std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>;

// Generate evenly sampled point at some discretization 'ds' between start and stop.
// ds must be > 0
PositionVector discretizePositions(const Eigen::Vector3d& start, const Eigen::Vector3d& stop, const double ds)
{
  auto dist = (stop - start).norm();

  size_t n_intermediate_points = 0;
  if (dist > ds)
    n_intermediate_points = static_cast<size_t>(std::lround(dist / ds));

  const auto total_points = 2 + n_intermediate_points;

  PositionVector result;
  result.reserve(total_points);

  for (std::size_t i = 0; i < total_points; ++i)
  {
    const double r = i / static_cast<double>(total_points - 1);
    Eigen::Vector3d point = start + (stop - start) * r;
    result.push_back(point);
  }
  return result;
}

int getDiscretizeNum(const double dist, const double ds)
{
  int n_intermediate_points = 0;
  if (dist > ds)
  {
    n_intermediate_points = std::lround(dist / ds);
  }

  int total_points_num = 2 + n_intermediate_points;
  return total_points_num;
}

Eigen::Affine3d makePose(const Eigen::Vector3d& position, const Eigen::Matrix3d& orientation,
                         const double z_axis_angle)
{
  Eigen::Affine3d m = Eigen::Affine3d::Identity();
  m.matrix().block<3,3>(0,0) = orientation;
  m.matrix().col(3).head<3>() = position;

  Eigen::AngleAxisd z_rot (z_axis_angle, Eigen::Vector3d::UnitZ());

  return m * z_rot;
}

descartes_planner::LadderGraph sampleSingleConfig(const descartes_core::RobotModel& model, const PositionVector& ps,
                                                  const double dt, const Eigen::Matrix3d& orientation,
                                                  const double z_axis_angle)
{
  descartes_planner::LadderGraph graph {model.getDOF()};
  graph.resize(ps.size());

  const descartes_core::TimingConstraint timing (dt);

  // Solve IK for each point
  for (std::size_t i = 0; i < ps.size(); ++i) //const auto& point : ps)
  {
    const auto& point = ps[i];
    const Eigen::Affine3d pose = makePose(point, orientation, z_axis_angle);
    std::vector<std::vector<double>> joint_poses;
    model.getAllIK(pose, joint_poses);
    graph.assignRung(i, descartes_core::TrajectoryID::make_nil(), timing, joint_poses);
  }

  bool has_edges_t = true;

  // now we have a graph with data in the 'rungs' and we need to compute the edges
  for (std::size_t i = 0; i < graph.size() - 1; ++i)
  {
    const auto start_idx = i;
    const auto end_idx = i + 1;
    const auto& joints1 = graph.getRung(start_idx).data;
    const auto& joints2 = graph.getRung(end_idx).data;
    const auto& tm = graph.getRung(end_idx).timing;
    const auto dof = model.getDOF();

    const auto start_size = joints1.size() / dof;
    const auto end_size = joints2.size() / dof;

    descartes_planner::DefaultEdgesWithTime builder (start_size, end_size, dof, tm.upper,
                                                     model.getJointVelocityLimits());
    for (size_t k = 0; k < start_size; k++) // from rung
    {
      const auto start_index = k * dof;

      for (size_t j = 0; j < end_size; j++) // to rung
      {
        const auto end_index = j * dof;

        builder.consider(&joints1[start_index], &joints2[end_index], j);
      }
      builder.next(k);
    }

    std::vector<descartes_planner::LadderGraph::EdgeList> edges = builder.result();
    if (!builder.hasEdges())
    {
//      ROS_WARN("No edges");
    }
    has_edges_t = has_edges_t && builder.hasEdges();
    graph.assignEdges(i, std::move(edges));
  } // end edge loop


//  if (has_edges_t) ROS_ERROR("Lots of edges");
//  else ROS_ERROR("No edges...");
  return graph;
}

void concatenate(descartes_planner::LadderGraph& dest, const descartes_planner::LadderGraph& src)
{
  assert(dest.size() > 0);
  assert(src.size() > 0);
  assert(dest.size() == src.size()); // same number of rungs
  // TODO
  // Combines the joints and edges from one graph, src, into another, dest
  // The joints copy straight over, but the index of the points needs to be transformed
  for (std::size_t i = 0; i < src.size(); ++i)
  {
    // Copy the joints
    auto& dest_joints = dest.getRung(i).data;
    const auto& src_joints = src.getRung(i).data;
    dest_joints.insert(dest_joints.end(), src_joints.begin(), src_joints.end());

    if (i != src.size() - 1)
    {
      // Copy the edges and transform them
      const auto next_rung_size = dest.rungSize(i + 1);
      auto& dest_edges = dest.getEdges(i);
      const auto& src_edges = src.getEdges(i);
      for (const auto& edge : src_edges)
      {
        auto edge_copy = edge;
        for (auto& e : edge_copy)
          e.idx += next_rung_size;
        dest_edges.push_back(edge_copy);
      }
    }
  }
}

} // end anon utility function ns

descartes_planner::LadderGraph descartes_planner::sampleConstrainedPaths(const descartes_core::RobotModel& model,
                                                                         ConstrainedSegment& segment)
{
  // Determine the linear points
  auto points = discretizePositions(segment.start, segment.end, segment.linear_disc);
  segment.process_pt_num = points.size();
  segment.retract_start_pt_num = 0;
  segment.retract_end_pt_num = 0;

  // Compute the number of angle steps
  static const auto min_angle = -M_PI;
  static const auto max_angle = M_PI;
  const auto n_angle_disc = std::lround( (max_angle - min_angle) / segment.z_axis_disc);
  const auto angle_step = (max_angle - min_angle) / n_angle_disc;

  // Compute the expected time step for each linear point
  double traverse_length = (segment.end - segment.start).norm();
  const auto dt =  traverse_length / segment.linear_vel;

  // fill in retract pts number once for later path re-identification
//  segment.retract_start_pt_num = getDiscretizeNum(segment.retract_dist, segment.linear_disc);
//  segment.retract_end_pt_num = getDiscretizeNum(segment.retract_dist, segment.linear_disc);

  LadderGraph graph {model.getDOF()};
  // there will be a ladder rung for each point that we must solve
  graph.resize(points.size());

  ROS_INFO_STREAM("Point has " << segment.orientations.size() << " orientations");
  // We will build up our graph one configuration at a time: a configuration is a single orientation and z angle disc
  std::vector<LadderGraph> graphs (segment.orientations.size() * n_angle_disc, LadderGraph(model.getDOF()));

  #pragma omp parallel for
  for (std::size_t j = 0; j < segment.orientations.size(); ++j) //const auto& orientation : segment.orientations)
  {
    const auto& orientation = segment.orientations[j];

//    ROS_INFO_STREAM("Orientation:\n" << orientation);

    // add retract pts according to orientation
    const PositionVector& process_pts = points;

    for (long i = 0; i < n_angle_disc; ++i)
    {
      const auto angle = angle_step * i;
//      LadderGraph single_config_graph = sampleSingleConfig(model, process_pts, dt, orientation, angle);
      graphs[j * n_angle_disc + i] =  sampleSingleConfig(model, process_pts, dt, orientation, angle);
    }
  }

  // reduction
  for (const auto& g : graphs)
  {
    concatenate(graph, g);
  }

//  concatenate(graph, single_config_graph);

  return graph;
}

void descartes_planner::appendInTime(LadderGraph &current, const LadderGraph &next)
{
  const auto ref_size = current.size();
  const auto new_total_size = ref_size + next.size();

  // So step 1 is to literally add the two sets of joints and edges together to make
  // one longer graph
  current.resize(new_total_size);
  for (std::size_t i = 0; i < next.size(); ++i)
  {
    current.getRung(ref_size + i) = next.getRung(i);
  }

  // The second problem is that we now need to 'connect' the two graphs.
  if (ref_size > 0 && next.size() > 0)
  {
    const auto dof = current.dof();
    auto& a_rung = current.getRung(ref_size - 1);
    auto& b_rung = current.getRung(ref_size);

    const auto n_start = a_rung.data.size() / dof;
    const auto n_end = b_rung.data.size() / dof;

    descartes_planner::DefaultEdgesWithoutTime builder (n_start, n_end, dof);

    for (size_t k = 0; k < n_start; k++) // from rung
    {
      const auto start_index = k * dof;

      for (size_t j = 0; j < n_end; j++) // to rung
      {
        const auto end_index = j * dof;

        builder.consider(&a_rung.data[start_index], &b_rung.data[end_index], j);
      }
      builder.next(k);
    }

    std::vector<descartes_planner::LadderGraph::EdgeList> edges = builder.result();
    current.assignEdges(ref_size - 1, std::move(edges));
  }

}
