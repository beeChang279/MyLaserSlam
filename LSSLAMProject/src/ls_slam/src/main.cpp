#include <gaussian_newton.h>
#include <readfile.h>

#include <ros/ros.h>
#include <visualization_msgs/MarkerArray.h>

#include <chrono>
#include "ceres_opt.h"
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/core/sparse_optimizer_terminate_action.h>
#include <g2o/types/slam2d/vertex_se2.h>
#include <g2o/types/slam2d/edge_se2.h>

//#define GUASSIAN_NEWTON
//#define CERES
#define G2O

//for visual
void PublishGraphForVisulization(ros::Publisher* pub,
                                 std::vector<Eigen::Vector3d>& Vertexs,
                                 std::vector<Edge>& Edges,
                                 int color = 0)
{
    visualization_msgs::MarkerArray marray;

    //point--red
    visualization_msgs::Marker m;
    m.header.frame_id = "map";
    m.header.stamp = ros::Time::now();
    m.id = 0;
    m.ns = "ls-slam";
    m.type = visualization_msgs::Marker::SPHERE;
    m.pose.position.x = 0.0;
    m.pose.position.y = 0.0;
    m.pose.position.z = 0.0;
    m.scale.x = 0.1;
    m.scale.y = 0.1;
    m.scale.z = 0.1;

    if(color == 0)
    {
        m.color.r = 1.0;
        m.color.g = 0.0;
        m.color.b = 0.0;
    }
    else
    {
        m.color.r = 0.0;
        m.color.g = 1.0;
        m.color.b = 0.0;
    }

    m.color.a = 1.0;
    m.lifetime = ros::Duration(0);

    //linear--blue
    visualization_msgs::Marker edge;
    edge.header.frame_id = "map";
    edge.header.stamp = ros::Time::now();
    edge.action = visualization_msgs::Marker::ADD;
    edge.ns = "karto";
    edge.id = 0;
    edge.type = visualization_msgs::Marker::LINE_STRIP;
    edge.scale.x = 0.1;
    edge.scale.y = 0.1;
    edge.scale.z = 0.1;

    if(color == 0)
    {
        edge.color.r = 0.0;
        edge.color.g = 0.0;
        edge.color.b = 1.0;
    }
    else
    {
        edge.color.r = 1.0;
        edge.color.g = 0.0;
        edge.color.b = 1.0;
    }
    edge.color.a = 1.0;

    m.action = visualization_msgs::Marker::ADD;
    uint id = 0;

    //加入节点
    for (uint i=0; i<Vertexs.size(); i++)
    {
        m.id = id;
        m.pose.position.x = Vertexs[i](0);
        m.pose.position.y = Vertexs[i](1);
        marray.markers.push_back(visualization_msgs::Marker(m));
        id++;
    }

    //加入边
    for(int i = 0; i < Edges.size();i++)
    {
        Edge tmpEdge = Edges[i];
        edge.points.clear();

        geometry_msgs::Point p;
        p.x = Vertexs[tmpEdge.xi](0);
        p.y = Vertexs[tmpEdge.xi](1);
        edge.points.push_back(p);

        p.x = Vertexs[tmpEdge.xj](0);
        p.y = Vertexs[tmpEdge.xj](1);
        edge.points.push_back(p);
        edge.id = id;

        marray.markers.push_back(visualization_msgs::Marker(edge));
        id++;
    }

    pub->publish(marray);
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "ls_slam");

    ros::NodeHandle nodeHandle;

    // beforeGraph
    ros::Publisher beforeGraphPub,afterGraphPub;
    beforeGraphPub = nodeHandle.advertise<visualization_msgs::MarkerArray>("beforePoseGraph",1,true);
    afterGraphPub  = nodeHandle.advertise<visualization_msgs::MarkerArray>("afterPoseGraph",1,true);


    std::string VertexPath = "/home/beechang/shenlanHW/HW6/LSSLAMProject/src/ls_slam/data/test_quadrat-v.dat";
    std::string EdgePath = "/home/beechang/shenlanHW/HW6/LSSLAMProject/src/ls_slam/data/test_quadrat-e.dat";

//    std::string VertexPath = "/home/eventec/LSSLAMProject/src/ls_slam/data/intel-v.dat";
//    std::string EdgePath = "/home/eventec/LSSLAMProject/src/ls_slam/data/intel-e.dat";

    std::vector<Eigen::Vector3d> Vertexs;
    std::vector<Edge> Edges;

    ReadVertexInformation(VertexPath,Vertexs);
    ReadEdgesInformation(EdgePath,Edges);

    PublishGraphForVisulization(&beforeGraphPub,
                                Vertexs,
                                Edges);

    double initError = ComputeError(Vertexs,Edges);
    std::cout << "initError:" << initError << std::endl;


#if defined CERES

    //ceres优化
    // 构建最小二乘问题
    ceres::Problem problem;
    ceres::LossFunction* loss_function = NULL;
    ceres::LocalParameterization* angle_local_parameterization =
        AngleLocalParameterization::Create();
    for (int i=0; i < Edges.size(); ++i)
    {
        Edge tmpEdge = Edges[i];

        const Eigen::Matrix3d sqrt_information =
            tmpEdge.infoMatrix.llt().matrixL();
        // Ceres will take ownership of the pointer.
        ceres::CostFunction* cost_function = PoseGraph2dErrorTerm::Create(
            tmpEdge.measurement(0), tmpEdge.measurement(1), tmpEdge.measurement(2), sqrt_information);
        problem.AddResidualBlock(cost_function,
                                loss_function,
                                &Vertexs[tmpEdge.xi](0),
                                &Vertexs[tmpEdge.xi](1),
                                &Vertexs[tmpEdge.xi](2),
                                &Vertexs[tmpEdge.xj](0),
                                &Vertexs[tmpEdge.xj](1),
                                &Vertexs[tmpEdge.xj](2));

        problem.SetParameterization(&Vertexs[tmpEdge.xi](2),
                                    angle_local_parameterization);
        problem.SetParameterization(&Vertexs[tmpEdge.xj](2),
                                    angle_local_parameterization);
    }
    problem.SetParameterBlockConstant(&Vertexs[0](0));
    problem.SetParameterBlockConstant(&Vertexs[0](1));
    problem.SetParameterBlockConstant(&Vertexs[0](2));

    ceres::Solver::Options options;
    options.max_num_iterations = 100;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    options.gradient_tolerance = 10e-4;
    // options.function_tolerance = 10e-4;
    // options.parameter_tolerance = 10e-4;
    options.trust_region_strategy_type = ceres::DOGLEG;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::cout << summary.FullReport() << '\n';

#endif

#if defined G2O

    typedef g2o::BlockSolver<g2o::BlockSolverTraits<3, 3>> Block; // 每个误差项优化变量维度为3，误差值维度为3
    //1、创建一个线性求解器LinearSolver，这里使用CSparse法，继承自LinearSolverCCS
    //Block::LinearSolverType* linearSolver = new g2o::LinearSolverCSparse<Block::PoseMatrixType>();
    std::unique_ptr<Block::LinearSolverType> linearSolver (new g2o::LinearSolverCSparse<Block::PoseMatrixType>());
    //2、创建BlockSolver，并用上面定义的线性求解器初始化
    //Block* solver_ptr = new Block(linearSolver);
    std::unique_ptr<Block> solver_ptr(new Block(std::move(linearSolver))); 
    //3、创建总求解器solver，并从GN, LM, DogLeg中选一个，再用上述块求解器BlockSolver初始化
    g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton(std::move(solver_ptr));
    //4、创建终极大boss，稀疏优化器（SparseOptimizer）
    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm(solver);

    //5、定义图的顶点和边，并添加到SparseOptimizer中
    for (size_t i = 0; i < Vertexs.size(); i++) {
        g2o::VertexSE2* v = new g2o::VertexSE2();
        v->setEstimate(Vertexs[i]);
        v->setId(i);
        if (i == 0) {
            v->setFixed(true);
        }
        optimizer.addVertex(v);
    }

    for (size_t i = 0; i < Edges.size(); i++) {
        g2o::EdgeSE2* edge = new g2o::EdgeSE2();

        Edge tmpEdge = Edges[i];

        edge->setId(i);
        edge->setVertex(0, optimizer.vertices()[tmpEdge.xi]);
        edge->setVertex(1, optimizer.vertices()[tmpEdge.xj]);

        edge->setMeasurement(tmpEdge.measurement);
        edge->setInformation(tmpEdge.infoMatrix);
        optimizer.addEdge(edge);
    }

    //6、设置优化参数，开始执行优化
    optimizer.setVerbose(true);
    optimizer.initializeOptimization();
    g2o::SparseOptimizerTerminateAction* terminateAction = new g2o::SparseOptimizerTerminateAction;
    terminateAction->setGainThreshold(1e-4);
    optimizer.addPostIterationAction(terminateAction);
    optimizer.optimize(100);

    for (size_t i = 0; i < Vertexs.size(); i++) {
        g2o::VertexSE2* v = static_cast<g2o::VertexSE2*>(optimizer.vertices()[i]);
        Vertexs[i] = v->estimate().toVector();
    }

#endif

#if defined GUASSIAN_NEWTON
    int maxIteration = 100;
    double epsilon = 1e-4;
    // std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
    for(int i = 0; i < maxIteration;i++)
    {
        std::cout << "Iterations:" << i << std::endl;
        
        Eigen::VectorXd dx = LinearizeAndSolve(Vertexs,Edges);
        
        //进行更新
        //TODO--Start
        for(int i = 0; i < Vertexs.size(); i++)
        {
            Vertexs[i] += dx.block(i * 3, 0, 3, 1);//更新回环中的每个点的位姿
            if (Vertexs[i](2) > M_PI)
                Vertexs[i](2) -= 2 * M_PI;
            else if (Vertexs[i](2) < -M_PI)
                Vertexs[i](2) += 2 * M_PI;
        }
        //TODO--End

        double maxError = -1;
        for(int k = 0; k < 3 * Vertexs.size();k++)
        {
            if(maxError < std::fabs(dx(k)))
            {
                maxError = std::fabs(dx(k));
            }
        }

        if(maxError < epsilon)
            break;
    }
    // std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
    // std::chrono::duration<double> time_used = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    // std::cout << "总用时: " << time_used.count() << std::endl;
#endif
    double finalError  = ComputeError(Vertexs,Edges);

    std::cout <<"FinalError:"<<finalError<<std::endl;

    PublishGraphForVisulization(&afterGraphPub,
                                Vertexs,
                                Edges,1);

    ros::spin();

    return 0;
}




