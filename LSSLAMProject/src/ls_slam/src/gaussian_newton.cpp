#include "gaussian_newton.h"
#include <eigen3/Eigen/Jacobi>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Householder>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>
#include <chrono>

#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseCholesky>

#include <iostream>

//位姿-->转换矩阵
Eigen::Matrix3d PoseToTrans(Eigen::Vector3d x)
{
    Eigen::Matrix3d trans;
    trans << cos(x(2)),-sin(x(2)),x(0),
             sin(x(2)), cos(x(2)),x(1),
                     0,         0,    1;

    return trans;
}


//转换矩阵－－＞位姿
Eigen::Vector3d TransToPose(Eigen::Matrix3d trans)
{
    Eigen::Vector3d pose;
    pose(0) = trans(0,2);
    pose(1) = trans(1,2);
    pose(2) = atan2(trans(1,0),trans(0,0));

    return pose;
}

//计算整个pose-graph的误差
double ComputeError(std::vector<Eigen::Vector3d>& Vertexs,
                    std::vector<Edge>& Edges)
{
    double sumError = 0;
    for(int i = 0; i < Edges.size();i++)
    {
        Edge tmpEdge = Edges[i];
        Eigen::Vector3d xi = Vertexs[tmpEdge.xi];
        Eigen::Vector3d xj = Vertexs[tmpEdge.xj];
        Eigen::Vector3d z = tmpEdge.measurement;
        Eigen::Matrix3d infoMatrix = tmpEdge.infoMatrix;

        Eigen::Matrix3d Xi = PoseToTrans(xi);
        Eigen::Matrix3d Xj = PoseToTrans(xj);
        Eigen::Matrix3d Z  = PoseToTrans(z);

        Eigen::Matrix3d Ei = Z.inverse() *  Xi.inverse() * Xj;

        Eigen::Vector3d ei = TransToPose(Ei);
        
        sumError += ei.transpose() * infoMatrix * ei;
    }
    return sumError;
}


/**
 * @brief CalcJacobianAndError
 *         计算jacobian矩阵和error
 * @param xi    fromIdx
 * @param xj    toIdx
 * @param z     观测值:xj相对于xi的坐标
 * @param ei    计算的误差
 * @param Ai    相对于xi的Jacobian矩阵
 * @param Bi    相对于xj的Jacobian矩阵
 */
void CalcJacobianAndError(Eigen::Vector3d xi,Eigen::Vector3d xj,Eigen::Vector3d z,
                          Eigen::Vector3d& ei,Eigen::Matrix3d& Ai,Eigen::Matrix3d& Bi)
{
    //TODO--Start
    // std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();    
    //定义变量
    double angle = 0;
    ei.setZero();
    Ai.setZero();
    Bi.setZero();
    Eigen::Vector2d ti(0, 0), tj(0, 0), tij(0, 0), ei_t(0, 0);
    Eigen::Matrix2d Ri;
    Eigen::Matrix2d Rij;
    Eigen::Matrix2d dRi_T_theta;
  
    //初始化值
    Ri.setZero();//xi到世界坐标系的旋转矩阵
    Rij.setZero();//xj到xi坐标系的旋转矩阵
    dRi_T_theta.setZero();//Ri的转置对theta_i的偏导数
    ti(0) = xi(0);//xi的平移向量x分量
    ti(1) = xi(1);//xi的平移向量y分量
    tj(0) = xj(0);//xj的平移向量x分量
    tj(1) = xj(1);//xj的平移向量y分量
    tij(0) = z(0);//xj相对于xi的平移向量x分量
    tij(1) = z(1);//xj相对于xi的平移向量y分量
    Ri << cos(xi(2)), -sin(xi(2)),
            sin(xi(2)), cos(xi(2));
       
    Rij << cos(z(2)), -sin(z(2)),
            sin(z(2)), cos(z(2)); 
    ei_t = Rij.transpose() * ((Ri.transpose() * (tj - ti)) - tij);
    ei(0) = ei_t(0);
    ei(1) = ei_t(1);
    angle = xj(2) - xi(2) - z(2);
    if(angle > M_PI)//角度归一化，角度未归一化图形变形严重
    {
        angle -= 2 * M_PI;
    }else if(angle < -M_PI)
    {
        angle += 2 * M_PI;
    }      
    ei(2) = angle;

    dRi_T_theta << -sin(xi(2)), cos(xi(2)),
                   -cos(xi(2)), -sin(xi(2));
    Ai.block(0, 0, 2, 2) = -Rij.transpose() * Ri.transpose();
    Ai.block(0, 2, 2, 1) = Rij.transpose() * dRi_T_theta * (tj - ti);
    Ai(2, 0) = 0;
    Ai(2, 1) = 0;
    Ai(2, 2) = -1;

    Bi.block(0, 0, 2, 2) = Rij.transpose() * Ri.transpose();
    Bi(0, 2) = 0;
    Bi(1, 2) = 0;
    Bi(2, 0) = 0;
    Bi(2, 1) = 0;
    Bi(2, 2) = 1;
    // std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
    // std::chrono::duration<double> time_used = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    // std::cout << "雅可比和误差求解用时: " << time_used.count() << std::endl;
    //求解一次雅可比和误差用时约0.000038
    //TODO--end
}

/**
 * @brief LinearizeAndSolve
 *        高斯牛顿方法的一次迭代．
 * @param Vertexs   图中的所有节点
 * @param Edges     图中的所有边
 * @return          位姿的增量
 */
Eigen::VectorXd LinearizeAndSolve(std::vector<Eigen::Vector3d>& Vertexs,
                                   std::vector<Edge>& Edges)
{
    //申请内存
    Eigen::MatrixXd H(Vertexs.size() * 3,Vertexs.size() * 3);
    Eigen::VectorXd b(Vertexs.size() * 3);

    H.setZero();
    b.setZero();

    //固定第一帧
    Eigen::Matrix3d I;
    I.setIdentity();
    H.block(0,0,3,3) += I;

    //构造H矩阵　＆ b向量
    for(int i = 0; i < Edges.size();i++)
    {
        //提取信息
        Edge tmpEdge = Edges[i];
        Eigen::Vector3d xi = Vertexs[tmpEdge.xi];
        Eigen::Vector3d xj = Vertexs[tmpEdge.xj];
        Eigen::Vector3d z = tmpEdge.measurement;
        Eigen::Matrix3d infoMatrix = tmpEdge.infoMatrix;

        //计算误差和对应的Jacobian
        Eigen::Vector3d ei; //eij
        Eigen::Matrix3d Ai; //Aij
        Eigen::Matrix3d Bi; //Bij
        CalcJacobianAndError(xi,xj,z,ei,Ai,Bi);

        //TODO--Start
        // std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
#if 0
        //对当前回环中的位姿求雅可比矩阵Jij，Jij是一个3行，当前回环位姿点数*3的矩阵
        //该写法低数据量有效，killian数据集H和b计算一次花费近2.3s，intel数据集H和b计算一次花费近0.3s，可能是J和H太大了
        Eigen::MatrixXd Jij(3, Vertexs.size() * 3);
        Jij.setZero();

        Jij.block(0, tmpEdge.xi * 3, 3, 3) = Ai;
        Jij.block(0, tmpEdge.xj * 3, 3, 3) = Bi;

        //对当前回环位姿求H和b，Hij = Jij_T * 信息矩阵 * Jij; bij = Jij_T * 信息矩阵 * eij
        H += (Jij.transpose() * infoMatrix * Jij);
        b += (Jij.transpose() * infoMatrix * ei);
#endif

#if 1  
        //该写法对于大数据量killian数据集仍能保持H和b计算一次约0.000071s
       H.block(3*tmpEdge.xi,3*tmpEdge.xi,3,3) += Ai.transpose()*infoMatrix*Ai;
       H.block(3*tmpEdge.xi,3*tmpEdge.xj,3,3) += Ai.transpose()*infoMatrix*Bi;
       H.block(3*tmpEdge.xj,3*tmpEdge.xi,3,3) += Bi.transpose()*infoMatrix*Ai;
       H.block(3*tmpEdge.xj,3*tmpEdge.xj,3,3) += Bi.transpose()*infoMatrix*Bi;
       b.block(3*tmpEdge.xi,0,3,1) += Ai.transpose()*infoMatrix*ei;
       b.block(3*tmpEdge.xj,0,3,1) += Bi.transpose()*infoMatrix*ei;
#endif
        // std::cout << "Edges.size(): " << Edges.size() << std::endl;
        // std::cout << "i: " << i << std::endl;
        // std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
        // std::chrono::duration<double> time_used = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
        // std::cout << "H和b求解用时: " << time_used.count() << std::endl;
        //求解一次H和b用时约0.000071
        //TODO--End
    }

    //求解
    Eigen::VectorXd dx;
 
    //TODO--Start
    // std::cout << "线性方程组求解中..." << std::endl;
    // std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
    //稠密矩阵解法，对于killian和intel数据集会一直卡在此处计算，耗时严重，killian数据集求解一次约2042s(半小时)
    //test_quadrat直接求解,用时约0.0002，3次迭代收敛总用时0.0020682,要求可逆
    //dx = H.inverse() * (-b);
    
    //QR分解,test_quadrat数据集用时约0.000242818，3次迭代收敛总用时0.00239064，无
    //dx = H.colPivHouseholderQr().solve(-b); 
   
    //choleskey分解,test_quadrat数据集用时约0.000139042，3次迭代收敛总用时0.00190253，要求正定
    //dx = (H.transpose() * H).llt().solve(H.transpose() * (-b));

    //改进的choleskey分解,test_quadrat数据集用时约0.000157566，4次迭代完成总用时0.002577，要求正定或负半定
    //dx = (H.transpose() * H).ldlt().solve(H.transpose() * (-b));

    //eigen稀疏矩阵求解方法，test_quadrat数据集求解一次约0.00015407s，3次迭代收敛，总耗时0.00194703；
    //killian数据集求解一次约0.878693s，6次迭代收敛，总耗时10.3719
    //intel数据集求解一次约0.176618s，5次迭代收敛，总耗时2.54404
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;
    Eigen::SparseMatrix<double> A = H.sparseView();
    dx = solver.compute(A).solve(-b);
    
    // std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
    // std::chrono::duration<double> time_used = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    // std::cout << "线性方程组求解用时: " << time_used.count() << std::endl;
    //TODO-End

    return dx;
}











