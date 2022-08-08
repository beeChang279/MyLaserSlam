#include "gaussian_newton.h"
#include <eigen3/Eigen/Jacobi>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Householder>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>

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
        std::cout << "sumError: " << sumError << std::endl;
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
{std::cout << "xi: " << xi << std::endl;
std::cout << "xj: " << xj << std::endl;
    //TODO--Start
#if 1
    ei.setZero();
    Ai.setZero();
    Bi.setZero();
    Eigen::Vector2d ti(0, 0), tj(0, 0), tij(0, 0), ei_t(0, 0);
    Eigen::Matrix2d Ri;
    Eigen::Matrix2d Rij;
    Eigen::Matrix2d dRi_T_theta;
    Ri.setZero();
    Rij.setZero();
    dRi_T_theta.setZero();
    ti(0) = xi(0);
    ti(1) = xi(1);
    tj(0) = xj(0);
    tj(1) = xj(1);
    tij(0) = z(0);
    tij(1) = z(1);
    Ri << cos(xi(2)), -sin(xi(2)),
            sin(xi(2)), cos(xi(2));
   
    Rij << cos(z(2)), -sin(z(2)),
            sin(z(2)), cos(z(2));
    
std::cout << "Rij: " << Rij << std::endl;    

    ei_t = Rij.transpose() * ((Ri.transpose() * (tj - ti)) - tij);
    ei(0) = ei_t(0);
    ei(1) = ei_t(1);
    ei(2) = xj(2) - xi(2) - z(2);
std::cout << "ei: " << ei << std::endl;
    dRi_T_theta << -(sin(xi(2))), cos(xi(2)),
                   -(cos(xi(2))), -(sin(xi(2)));
    Ai.block(0, 0, 2, 2) = -Rij.inverse() * Ri.inverse();
    Ai.block(0, 2, 2, 1) = Rij.inverse() * dRi_T_theta * (tj - ti);
    Ai(2, 0) = 0;
    Ai(2, 1) = 0;
    Ai(2, 2) = -1;

    Bi.block(0, 0, 2, 2) = Rij.inverse() * Ri.inverse();
    Bi(0, 2) = 0;
    Bi(1, 2) = 0;
    Bi(2, 0) = 0;
    Bi(2, 1) = 0;
    Bi(2, 2) = -1;
#endif
#if 0
    Eigen::Matrix2d Ri,Ri_t,Rj,Rj_t,Rij ,Rij_t ,d_Ri_T,e1;
    Eigen::Vector2d ti,tj,tij,e2;
    Eigen::Matrix3d Ti = PoseToTrans(xi);
    Eigen::Matrix3d Tj = PoseToTrans(xj);
    Eigen::Matrix3d Tij = Ti.transpose()*Tj;
    Ri  = Ti.block(0,0,2,2);
    Rj = Tj.block(0,0,2,2);
    Rij = Tij.block(0,0,2,2);
    std::cout << "xi: " << xi << std::endl;
    std::cout << "xj: " << xj << std::endl;
    std::cout << "Ti: " << Ti << std::endl;
    std::cout << "Tj: " << Tj << std::endl;
    std::cout << "Tij: " << Tij << std::endl;
    std::cout << "Ri: " << Ri << std::endl;
    std::cout << "Rij: " << Rij << std::endl;
    Ri_t = Ti.transpose().block(0,0,2,2);
    Rj_t = Tj.transpose().block(0,0,2,2);
    Rij_t = Tij.transpose().block(0,0,2,2);      
    ti =Eigen::Vector2d(xi(0),xi(1));
    tj =Eigen::Vector2d(xj(0),xj(1));
    tij=Eigen::Vector2d(z(0),z(1));
    d_Ri_T = Ti.transpose().block(0,0,2,2);
    Eigen::Vector2d error = Rij_t*(Ri_t*(tj- ti)-tij);
    double col = xj(2)-xi(2)-z(2);
    ei  = Eigen::Vector3d(error(0),error(1),col);
std::cout << "ei: " << ei << std::endl;
    e1 =  Rij_t*Ri_t;
    e2 =  Rij_t* d_Ri_T*( tj- ti);

    Ai.block(0,0,2,2) = -e1;
    Ai.block(0,2,2,1) = e2;
    Ai.block(2,0,1,3) = Eigen::Vector3d(0.0,0.0,-1.0).transpose();
    Bi.setIdentity();
    Bi.block(0,0,2,2) = e1;
#endif
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
        Eigen::MatrixXd Jij(3, Vertexs.size() * 3);
        Jij.setZero();
// std::cout << "Ai: " << Ai << std::endl;
// std::cout << "Bi: " << Bi << std::endl;
        Jij.block(0, tmpEdge.xi * 3, 3, 3) = Ai;
        Jij.block(0, tmpEdge.xj * 3, 3, 3) = Bi;
//std::cout << "Jij: " << Jij << std::endl;
        H += (Jij.transpose() * infoMatrix * Jij);
        b += (Jij.transpose() * infoMatrix * ei);
        //TODO--End
    }

    //求解
    Eigen::VectorXd dx;
std::cout << "H: " << H << std::endl;
std::cout << "b: " << b << std::endl;
    //TODO--Start
    dx = -H.inverse() * b;
std::cout << "dx: " << dx << std::endl;   
    //TODO-End

    return dx;
}











