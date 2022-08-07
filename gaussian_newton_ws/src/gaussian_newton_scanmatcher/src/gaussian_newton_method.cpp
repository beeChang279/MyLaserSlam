#include <map.h>
#include "gaussian_newton_method.h"

const double GN_PI = 3.1415926;
static double E_T = 0;

//进行角度正则化．
double GN_NormalizationAngle(double angle)
{
    if(angle > GN_PI)
        angle -= 2*GN_PI;
    else if(angle < -GN_PI)
        angle += 2*GN_PI;

    return angle;
}

Eigen::Matrix3d GN_V2T(Eigen::Vector3d vec)
{
    Eigen::Matrix3d T;
    T  << cos(vec(2)),-sin(vec(2)),vec(0),
            sin(vec(2)), cos(vec(2)),vec(1),
            0,           0,     1;

    return T;
}

//对某一个点进行转换．
Eigen::Vector2d GN_TransPoint(Eigen::Vector2d pt,Eigen::Matrix3d T)
{
    Eigen::Vector3d tmp_pt(pt(0),pt(1),1);
    tmp_pt = T * tmp_pt;
    return Eigen::Vector2d(tmp_pt(0),tmp_pt(1));
}



//用激光雷达数据创建势场．
map_t* CreateMapFromLaserPoints(Eigen::Vector3d map_origin_pt,
                                std::vector<Eigen::Vector2d> laser_pts,
                                double resolution)
{
    map_t* map = map_alloc();

    map->origin_x = map_origin_pt(0);
    map->origin_y = map_origin_pt(1);
    map->resolution = resolution;

    //固定大小的地图，必要时可以扩大．
    map->size_x = 10000;
    map->size_y = 10000;

    map->cells = (map_cell_t*)malloc(sizeof(map_cell_t)*map->size_x*map->size_y);

    //高斯平滑的sigma－－固定死
    map->likelihood_sigma = 0.5;

    Eigen::Matrix3d Trans = GN_V2T(map_origin_pt);

    //设置障碍物
    for(int i = 0; i < laser_pts.size();i++)
    {
        Eigen::Vector2d tmp_pt = GN_TransPoint(laser_pts[i],Trans);

        int cell_x,cell_y;
        cell_x = MAP_GXWX(map,tmp_pt(0));
        cell_y = MAP_GYWY(map,tmp_pt(1));

        map->cells[MAP_INDEX(map,cell_x,cell_y)].occ_state = CELL_STATUS_OCC;
    }

    //进行障碍物的膨胀--最大距离固定死．
    map_update_cspace(map,0.5);

    return map;
}


/**
 * @brief InterpMapValueWithDerivatives
 * 在地图上的进行插值，得到coords处的势场值和对应的关于位置的梯度．
 * 返回值为Eigen::Vector3d ans
 * ans(0)表示势场值
 * ans(1:2)表示梯度
 * @param map
 * @param coords
 * @return
 */
Eigen::Vector3d InterpMapValueWithDerivatives(map_t* map,Eigen::Vector2d& coords)
{
    Eigen::Vector3d ans;
    //TODO
    //坐标变换后的激光点在世界坐标系下的坐标转换到地图坐标系下，找到x，y的最近整数坐标(对应的栅格)
    double x = coords[0];
    double y = coords[1];

    int map_x0 = std::floor(MAP_GXWX_DOUBLE(map, x));//x下界
    int map_y0 = std::floor(MAP_GYWY_DOUBLE(map, y));//y下界
    int map_x1 = std::ceil(MAP_GXWX_DOUBLE(map, x));//x上界
    int map_y1 = std::ceil(MAP_GYWY_DOUBLE(map, y));//y上界

    //x，y附近四个栅格似然场值
    double z1 = map->cells[MAP_INDEX(map, map_x0, map_y0)].score;
    double z2 = map->cells[MAP_INDEX(map, map_x1, map_y0)].score;
    double z3 = map->cells[MAP_INDEX(map, map_x1, map_y1)].score;
    double z4 = map->cells[MAP_INDEX(map, map_x0, map_y1)].score;

    //四个栅格坐标转换到世界坐标下
    double x0 = MAP_WXGX(map, map_x0);
    double y0 = MAP_WYGY(map, map_y0);
    double x1 = MAP_WXGX(map, map_x1);
    double y1 = MAP_WYGY(map, map_y1);

    double dx0 = (x - x0) / (x1 - x0);
    double dy0 = (y - y0) / (y1 - y0);
    double dx1 = (x1 - x) / (x1 - x0);
    double dy1 = (y1 - y) / (y1 - y0);
  
    //通过插值函数求解该点似然场值和x方向的梯度和y方向的梯度
    ans(0) = dy0 * (dx0 * z3 + dx1 * z4) + dy1 * (dx0 * z2 + dx1 * z1);
    ans(1) = dy0 * ((z3 - z4) / (x1 - x0)) + dy1 * ((z2 - z1) / (x1 - x0));
    ans(2) = (dx0 * z3 + dx1 * z4) / (y1 - y0) - (dx0 * z2 - dx1 * z1) / (y1 - y0);

    //END OF TODO

    return ans;
}


/**
 * @brief ComputeCompleteHessianAndb
 * 计算H*dx = b中的H和b
 * @param map
 * @param now_pose
 * @param laser_pts
 * @param H
 * @param b
 */
void ComputeHessianAndb(map_t* map, Eigen::Vector3d now_pose,
                        std::vector<Eigen::Vector2d>& laser_pts,
                        Eigen::Matrix3d& H, Eigen::Vector3d& b)
{
    H = Eigen::Matrix3d::Zero();
    b = Eigen::Vector3d::Zero();

    //TODO
   
    //给定变换初值T，将里程计增量变换到雷达坐标系下，加上前一帧时刻雷达位姿得到
    Eigen::Matrix3d T;
    T.setZero();
    T << cos(now_pose(2)), -sin(now_pose(2)), now_pose(0),
        sin(now_pose(2)), cos(now_pose(2)), now_pose(1),
        0, 0, 1;

    Eigen::Vector3d LaserPose3d(0, 0, 0);
    Eigen::Vector3d Si_T(0, 0, 0);
    Eigen::Vector2d LaserPose2d(0, 0);
    Eigen::Matrix<double, 2, 3> dT_Si;
    Eigen::Vector3d Ans(0, 0, 0);
    Eigen::Vector2d dM_Si_T(0, 0);

    dT_Si.setZero();

    for(int i = 0; i < laser_pts.size(); i++){
        
        LaserPose3d << laser_pts[i](0), laser_pts[i](1), 1; //将当前帧激光雷达二维向量表示转为三维向量表示
        
        //求解Si(T)，第i个激光点位姿变换后的坐标
        Si_T = T * LaserPose3d;//第i个激光点初始位姿变换后的坐标
        LaserPose2d << Si_T(0), Si_T(1);//变换后的激光点转成二维向量，用于在似然场中进行拉格朗日插值
       
        //求解dSi(T)/dT
        dT_Si << 1, 0, (-sin(now_pose(2)) * LaserPose3d(0) - cos(now_pose(2)) * LaserPose3d(1)),
                 0, 1, (cos(now_pose(2)) * LaserPose3d(0) - sin(now_pose(2)) * LaserPose3d(1));//变换后的激光点对初始变换T的导数
        
        //求解M(Si(T))和dM(Si(T))
        //变换后的激光点坐标在似然场中做拉格朗日插值，返回Ans(0)-该点似然场值，Ans(1)-似然场在该点对x的偏导数(x方向上的梯度)，Ans(2)-似然场在该点对y的偏导数(y方向上的梯度)
        Ans = InterpMapValueWithDerivatives(map, LaserPose2d);
        dM_Si_T << Ans(1), Ans(2);//似然场对坐标位置的导数，由x，y两个维度的偏导数组成

        //累加计算H和b，Eigen::Vector为列向量，先转置
        b += (dM_Si_T.transpose() * dT_Si).transpose() * (1 - Ans(0));//计算b
        H += (dM_Si_T.transpose() * dT_Si).transpose() * (dM_Si_T.transpose() * dT_Si);//计算H，该H值用来近似海瑟矩阵H
        // std::cout << "laser_pts.size(): " << laser_pts.size() << std::endl; 
        // std::cout << "Ans(0): " << Ans(0) << std::endl; 
        // std::cout << "1-Ans(0): " << (1 - Ans(0)) << std::endl; 
        // std::cout << "pow((1 - Ans(0)), 2): " << pow((1 - Ans(0)), 2) << std::endl; 
        //计算目标函数值E(T):1-M(Si(T))的平方和，E_T为额外定义的静态全局变量static double E_T
        E_T += pow((1 - Ans(0)), 2);
        
    }
    
    //END OF TODO
}


/**
 * @brief GaussianNewtonOptimization
 * 进行高斯牛顿优化．
 * @param map
 * @param init_pose
 * @param laser_pts
 */
void GaussianNewtonOptimization(map_t*map,Eigen::Vector3d& init_pose,std::vector<Eigen::Vector2d>& laser_pts)
{
    //int Iteration_count = 0;
    int maxIteration = 20;
    Eigen::Vector3d now_pose = init_pose;
    double lastE_T = std::numeric_limits<double>::max();
    for(int i = 0; i < maxIteration;i++)
    {
        //Iteration_count = i;
        
        //TODO
        Eigen::Matrix3d H;
        Eigen::Vector3d b;
        Eigen::Vector3d delta_T;

        ComputeHessianAndb(map, now_pose, laser_pts, H, b);
        //std::cout << "E_T----------------: " << E_T << std::endl;
        if(lastE_T <= E_T || (E_T > 0 && E_T <= 0.001))//假定阈值为0.001，E_T为额外定义的静态全局变量static double E_T,表示目标函数值E(T):1-M(Si(T))的平方和
        {
            E_T = 0;
            break;
        }
        
        lastE_T = E_T;
        E_T = 0;
        delta_T = H.inverse() * b;//计算delta_T
        now_pose += delta_T;//T += delta_T,位姿变换增加一个梯度delta_T
        //END OF TODO

    }
    init_pose = now_pose;
    E_T = 0;
    //std::cout << "Iteration_count: " << Iteration_count << std::endl;
}
