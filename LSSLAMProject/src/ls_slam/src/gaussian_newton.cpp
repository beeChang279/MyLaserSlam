#include "gaussian_newton.h"
#include <eigen3/Eigen/Jacobi>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Householder>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>
#include <chrono>

#include <iostream>
const double GN_PI = 3.1415926;

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
    
    //定义变量
    double angle = 0;
    ei.setZero();
    Ai.setZero();
    Bi.setZero();
    Eigen::Vector2d ti(0, 0), tj(0, 0), tij(0, 0), ei_t(0, 0);
    Eigen::Matrix2d Ri;
    Eigen::Matrix2d Rij;
    Eigen::Matrix2d dRi_T_theta;
    Eigen::Matrix3d Ti;
    Eigen::Matrix3d Tj;
    
    //初始化值
    Ri.setZero();//xi到世界坐标系的旋转矩阵
    Rij.setZero();//xj到xi坐标系的旋转矩阵
    dRi_T_theta.setZero();//Ri的转置对theta_i的偏导数
    Ti.setZero();//xi到世界坐标系的变换矩阵
    Tj.setZero();//xj到世界坐标系的变换矩阵
    ti(0) = xi(0);//xi的平移向量x分量
    ti(1) = xi(1);//xi的平移向量y分量
    tj(0) = xj(0);//xj的平移向量x分量
    tj(1) = xj(1);//xj的平移向量y分量
    tij(0) = z(0);//xj相对于xi的平移向量x分量
    tij(1) = z(1);//xj相对于xi的平移向量y分量
    Ti = PoseToTrans(xi);
    Tj = PoseToTrans(xj);
    Ri << cos(xi(2)), -sin(xi(2)),
            sin(xi(2)), cos(xi(2));//Ri = Ti.block(0, 0, 2, 2);
   
    // Rij << cos(z(2)), -sin(z(2)),
    //         sin(z(2)), cos(z(2));//用z算出来的不对
    Rij = (Ti.inverse() * Tj).block(0, 0, 2, 2); //xj到世界坐标系再到xi坐标系  
    ei_t = Rij.transpose() * ((Ri.transpose() * (tj - ti)) - tij);
    ei(0) = ei_t(0);
    ei(1) = ei_t(1);
    angle = xj(2) - xi(2) - z(2);
    if(angle > GN_PI)//角度归一化，角度未归一化图形变形严重
    {
        angle -= 2 * GN_PI;
    }else if(angle < -GN_PI)
    {
        angle += 2 * GN_PI;
    }      
    ei(2) = angle;

    // dRi_T_theta << -(sin(xi(2))), cos(xi(2)),
    //                -(cos(xi(2))), -(sin(xi(2)));//Ri的转置不需要对theta_i求导？求完导出错
    dRi_T_theta << cos(xi(2)), sin(xi(2)),
                    -sin(xi(2)), cos(xi(2));
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
        
        //对当前回环中的位姿求雅可比矩阵Jij，Jij是一个3行，当前回环位姿点数*3的矩阵
        Eigen::MatrixXd Jij(3, Vertexs.size() * 3);
        Jij.setZero();

        Jij.block(0, tmpEdge.xi * 3, 3, 3) = Ai;
        Jij.block(0, tmpEdge.xj * 3, 3, 3) = Bi;

        //对当前回环位姿求H和b，Hij = Jij_T * 信息矩阵 * Jij; bij = Jij_T * 信息矩阵 * eij
        H += (Jij.transpose() * infoMatrix * Jij);
        b += (Jij.transpose() * infoMatrix * ei);

        //TODO--End
    }

    //求解
    Eigen::VectorXd dx;

    //TODO--Start
    //std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
    //dx = H.inverse() * (-b);//直接求解,用时约0.00012，58次迭代完成总用时0.0403664,要求可逆
    //dx = H.colPivHouseholderQr().solve(-b); //QR分解,用时约0.00018，58次迭代完成总用时0.0737936，无
    dx = (H.transpose() * H).llt().solve(H.transpose() * (-b));//choleskey分解,用时约0.00008，58次迭代完成总用时0.0357698，要求正定
    //dx = (H.transpose() * H).ldlt().solve(H.transpose() * (-b));//改进的choleskey分解,用时约0.00013，58次迭代完成总用时0.06514，要求正定或负半定
    // std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
    // std::chrono::duration<double> time_used = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    // std::cout << "线性方程组求解用时: " << time_used.count() << std::endl;
    //TODO-End

    return dx;
}











