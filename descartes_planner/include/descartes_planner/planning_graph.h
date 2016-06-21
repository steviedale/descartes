/*
 * Software License Agreement (Apache License)
 *
 * Copyright (c) 2014, Southwest Research Institute
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/*
 * Planning_graph.h
 *
 *  Created on: Jun 5, 2014
 *      Author: Dan Solomon
 */

#ifndef PLANNING_GRAPH_H_
#define PLANNING_GRAPH_H_

#include <boost/function.hpp>
#include "descartes_core/trajectory_pt.h"
#include "descartes_trajectory/cart_trajectory_pt.h"
#include "descartes_trajectory/joint_trajectory_pt.h"

#include "descartes_planner/ladder_graph.h"

namespace descartes_planner
{

typedef boost::function<double(const std::vector<double> &, const std::vector<double> &)> CostFunction;


class PlanningGraph
{
public:
  PlanningGraph(descartes_core::RobotModelConstPtr model);
  PlanningGraph(descartes_core::RobotModelConstPtr model, CostFunction cost_function_callback);

  virtual ~PlanningGraph();

  /** \brief Clear all previous graph data */
  void clear();

  /** @brief initial population of graph trajectory elements
   * @param points list of trajectory points to be used to construct the graph
   * @return True if the graph was successfully created
   */
  bool insertGraph(const std::vector<descartes_core::TrajectoryPtPtr> *points);

  /** @brief adds a single trajectory point to the graph
   * @param point The new point to add to the graph
   * @return True if the point was successfully added
   */
  bool addTrajectory(descartes_core::TrajectoryPtPtr point, descartes_core::TrajectoryPt::ID previous_id,
                     descartes_core::TrajectoryPt::ID next_id);

  bool modifyTrajectory(descartes_core::TrajectoryPtPtr point);

  bool removeTrajectory(descartes_core::TrajectoryPtPtr point);

  /** @brief Calculate and return the shortest path from the given joint solution indices
   * @param startIndex The index of the joint solution at which to start
   * @param endIndex The index of the joint solution at which to end
   * @param cost The cost of the returned path
   * @param path The sequence of points (joint solutions) for the path (TODO: change to JointTrajectoryPt?)
   * @return True if a valid path is found
   */
  bool getShortestPath(double &cost, std::list<descartes_trajectory::JointTrajectoryPt> &path);

  descartes_core::RobotModelConstPtr getRobotModel();

protected:
  descartes_planner::LadderGraph graph_;
  descartes_core::RobotModelConstPtr robot_model_;
  CostFunction custom_cost_function_;
  /**
   * @brief A pair indicating the validity of the edge, and if valid, the cost associated
   *        with that edge
   */
  typedef std::pair<bool, double> EdgeWeightResult;

  /** @brief function for computing edge weight based on specified cost function */
  EdgeWeightResult edgeWeight(const descartes_trajectory::JointTrajectoryPt &start,
                              const descartes_trajectory::JointTrajectoryPt &end) const;

  bool calculateJointSolutions(const descartes_core::TrajectoryPtPtr* points, const std::size_t count,
                               std::vector<std::vector<std::vector<double>>>& poses);

  /** @brief (Re)create the actual graph nodes(vertices) from the list of joint solutions (vertices) */
  bool populateGraphVertices(const std::vector<descartes_core::TrajectoryPtPtr> &points,
                             std::vector<std::vector<descartes_trajectory::JointTrajectoryPt>> &poses);

  /** @brief calculate weights fro each start point to each end point */
  std::vector<LadderGraph::EdgeList> calculateEdgeWeights(const std::vector<double> &start_joints,
                                                          const std::vector<double> &end_joints,
                                                          size_t dof);
};

} /* namespace descartes_planner */

#endif /* PLANNING_GRAPH_H_ */
