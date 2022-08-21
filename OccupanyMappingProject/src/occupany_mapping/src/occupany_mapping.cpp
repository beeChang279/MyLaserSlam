#include "occupany_mapping.h"
#include "nav_msgs/GetMap.h"
#include "sensor_msgs/PointCloud.h"
#include "sensor_msgs/PointCloud2.h"
#include "geometry_msgs/Point32.h"

#define GRID
//#define COUNTMODEL
//#define TSDF

/**
 * Increments all the grid cells from (x0, y0) to (x1, y1);
 * //不包含(x1,y1)
 * 2D画线算法　来进行计算两个点之间的grid cell
 * @param x0
 * @param y0
 * @param x1
 * @param y1
 */
std::vector<GridIndex> TraceLine(int x0, int y0, int x1, int y1)
{
    GridIndex tmpIndex;
    std::vector<GridIndex> gridIndexVector;

    bool steep = abs(y1 - y0) > abs(x1 - x0);
    if (steep)
    {
        std::swap(x0, y0);
        std::swap(x1, y1);
    }
    if (x0 > x1)
    {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    int deltaX = x1 - x0;
    int deltaY = abs(y1 - y0);
    int error = 0;
    int ystep;
    int y = y0;

    if (y0 < y1)
    {
        ystep = 1;
    }
    else
    {
        ystep = -1;
    }

    int pointX;
    int pointY;
    for (int x = x0; x <= x1; x++)
    {
        if (steep)
        {
            pointX = y;
            pointY = x;
        }
        else
        {
            pointX = x;
            pointY = y;
        }

        error += deltaY;

        if (2 * error >= deltaX)
        {
            y += ystep;
            error -= deltaX;
        }

        //不包含最后一个点．
        if (pointX == x1 && pointY == y1)
            continue;

        //保存所有的点
        tmpIndex.SetIndex(pointX, pointY);

        gridIndexVector.push_back(tmpIndex);
    }

    return gridIndexVector;
}

void SetMapParams(void)
{
    mapParams.width = 1000; //栅格/像素个数
    mapParams.height = 1000;
    mapParams.resolution = 0.05;

    //每次被击中的log变化值，覆盖栅格建图算法需要的参数
    mapParams.log_free = -1;
    mapParams.log_occ = 2;

    //每个栅格的最大最小值．
    mapParams.log_max = 100.0;
    mapParams.log_min = 0.0;

    mapParams.origin_x = 0.0;
    mapParams.origin_y = 0.0;

    //地图的原点，在地图的正中间
    mapParams.offset_x = 500;
    mapParams.offset_y = 500;

    pMap = new unsigned char[mapParams.width * mapParams.height];//用于栅格建图

    //计数建图算法需要的参数
    //每个栅格被激光击中的次数
    pMapHits = new unsigned long[mapParams.width * mapParams.height];
    //每个栅格被激光通过的次数
    pMapMisses = new unsigned long[mapParams.width * mapParams.height];

    //TSDF建图算法需要的参数
    pMapW = new unsigned long[mapParams.width * mapParams.height];
    pMapTSDF = new double[mapParams.width * mapParams.height];

    //初始化
    for (int i = 0; i < mapParams.width * mapParams.height; i++)
    {
        pMap[i] = 50;
        pMapHits[i] = 0;
        pMapMisses[i] = 0;
        pMapW[i] = 0;
        pMapTSDF[i] = -1;
    }
}

//从世界坐标系转换到栅格坐标系
GridIndex ConvertWorld2GridIndex(double x, double y)
{
    GridIndex index;

    index.x = std::ceil((x - mapParams.origin_x) / mapParams.resolution) + mapParams.offset_x;
    index.y = std::ceil((y - mapParams.origin_y) / mapParams.resolution) + mapParams.offset_y;

    return index;
}

int GridIndexToLinearIndex(GridIndex index)
{
    int linear_index;
    linear_index = index.y + index.x * mapParams.width;
}

//判断index是否有效
bool isValidGridIndex(GridIndex index)
{
    if (index.x >= 0 && index.x < mapParams.width && index.y >= 0 && index.y < mapParams.height)
        return true;

    return false;
}

void DestoryMap()
{
    if (pMap != NULL)
        delete pMap;
}

//
void OccupanyMapping(std::vector<GeneralLaserScan> &scans, std::vector<Eigen::Vector3d> &robot_poses)
{
    std::cout << "开始建图，请稍后..." << std::endl;
    //枚举所有的激光雷达数据

    for (int i = 0; i < scans.size(); i++)
    {
        GeneralLaserScan scan = scans[i];
        Eigen::Vector3d robotPose = robot_poses[i];

        //机器人的下标
        GridIndex robotIndex = ConvertWorld2GridIndex(robotPose(0), robotPose(1));

        for (int id = 0; id < scan.range_readings.size(); id++)
        {
            double dist = scan.range_readings[id];
            double angle = -scan.angle_readings[id]; // 激光雷达逆时针转，角度取反

            if (std::isinf(dist) || std::isnan(dist))
                continue;

            //计算得到该激光点的世界坐标系的坐标
            double theta = -robotPose(2); // 激光雷达逆时针转，角度取反
            double laser_x = dist * cos(angle);
            double laser_y = dist * sin(angle);

            double world_x = cos(theta) * laser_x - sin(theta) * laser_y + robotPose(0);
            double world_y = sin(theta) * laser_x + cos(theta) * laser_y + robotPose(1);

            //start of TODO 对对应的map的cell信息进行更新．（1,2,3题内容）
#if defined GRID
            //找到该激光点在栅格地图中的坐标(击中坐标)
            GridIndex hitGridIndex = ConvertWorld2GridIndex(world_x, world_y);
            if(true == isValidGridIndex(hitGridIndex)) //如果是有效坐标，寻找两个栅格坐标连线上的栅格坐标
            {
                //寻找两个栅格坐标连线上的栅格坐标
                std::vector<GridIndex> freeGridIndex = TraceLine(robotIndex.x, robotIndex.y, hitGridIndex.x, hitGridIndex.y);
                //对穿过的栅格计算Misses值
                for(int j = 0; j < freeGridIndex.size(); j++)
                {
                    int misId = GridIndexToLinearIndex(freeGridIndex[j]);
                    // pMap[misId] += mapParams.log_free;
                    // if(pMap[misId] < mapParams.log_min)
                    // {
                    //     pMap[misId] = mapParams.log_min; //栅格的log值不能小于最小值
                    // }
                    int value = pMap[misId];
                    value += mapParams.log_free;
                    if(value < mapParams.log_min)
                    {
                        pMap[misId] = mapParams.log_min;
                    }else{
                        pMap[misId] = value;
                    }
                }
                
                //对击中的栅格计算Hits值
                int hitId = GridIndexToLinearIndex(hitGridIndex);
                // pMap[hitId] += mapParams.log_occ;
                // if(pMap[hitId] > mapParams.log_max)
                // {
                //     pMap[hitId] = mapParams.log_max; //栅格的log值不能大于最大值
                // }
                int value = pMap[hitId];
                value += mapParams.log_occ;
                if(value > mapParams.log_max)
                {
                    pMap[hitId] = mapParams.log_max;
                }else{
                    pMap[hitId] = value;
                }
            }
#endif

#if defined COUNTMODEL
            //找到该激光点在栅格地图中的坐标(击中坐标)
            GridIndex hitGridIndex = ConvertWorld2GridIndex(world_x, world_y);
            if(true == isValidGridIndex(hitGridIndex)) //如果是有效坐标，寻找两个栅格坐标连线上的栅格坐标
            {
                //寻找两个栅格坐标连线上的栅格坐标
                std::vector<GridIndex> freeGridIndex = TraceLine(robotIndex.x, robotIndex.y, hitGridIndex.x, hitGridIndex.y);
                //对穿过的栅格累计Misses值
                for(int j = 0; j < freeGridIndex.size(); j++)
                {
                    int misId = GridIndexToLinearIndex(freeGridIndex[j]);
                    pMapMisses[misId]++;                    
                }
                
                //对击中的栅格累计Hits值
                int hitId = GridIndexToLinearIndex(hitGridIndex);
                pMapHits[hitId]++;
            }
#endif

#if defined TSDF
            double t = 0.1;//设置截断距离t
            GridIndex hitGridIndex = ConvertWorld2GridIndex(world_x, world_y);
            if(true == isValidGridIndex(hitGridIndex)) //如果是有效坐标，寻找两个栅格坐标连线上的栅格坐标
            {
                //寻找两个栅格坐标连线上的栅格坐标
                std::vector<GridIndex> freeGridIndex = TraceLine(robotIndex.x, robotIndex.y, hitGridIndex.x, hitGridIndex.y);
                //计算该束激光穿过的栅格的TSDF值
                for(int j = 0; j < freeGridIndex.size(); j++)
                {
                    //计算每一个栅格到激光坐标系原点的距离
                    
                    //计算栅格在世界坐标下的表示
                    double w_x = (freeGridIndex[j].x - mapParams.offset_x) * mapParams.resolution + mapParams.origin_x;
                    double w_y = (freeGridIndex[j].y - mapParams.offset_y) * mapParams.resolution + mapParams.origin_y;
                    
                    double delta_x = w_x - robotPose(0);//机器人(激光雷达)的位姿即是机器人(激光雷达)原点在世界坐标系下的坐标
                    double delta_y = w_y - robotPose(1);

                    //计算sdf
                    double dGrid2Laser = sqrt((pow(delta_x, 2) + pow(delta_y, 2)));//栅格到激光雷达坐标系原点的距离
                    double dLaser = dist;//激光雷达测量距离
                    double sdf = dLaser - dGrid2Laser;

                    //计算tsdf
                    double tempMin = sdf / t;
                    if(1 < tempMin)
                    {
                        tempMin = 1;
                    }

                    double tempMax = tempMin;
                    if(-1 > tempMax)
                    {
                        tempMax = -1;
                    }
                    double tsdf = tempMax;

                    int tsdfId = GridIndexToLinearIndex(freeGridIndex[j]);//栅格坐标转一维线性坐标
                    //计算更新TSDFi(x),pMapW为该栅格累计观察次数，pMapTSDF为该栅格累计的平均值
                    pMapTSDF[tsdfId] = (pMapTSDF[tsdfId] * pMapW[tsdfId] + tsdf) / (pMapW[tsdfId] + 1);
                    pMapW[tsdfId]++;
                }
            }

#endif
            //end of TODO
        }
    }
    //start of TODO 通过计数建图算法或TSDF算法对栅格进行更新（2,3题内容）
#if defined COUNTMODEL
    //根据pMapHits和pMapMisses计算pMap
    for(int i = 0; i < mapParams.width * mapParams.height; i++)
    {
        if(0 != pMapHits[i] || 0 != pMapMisses[i])
        {
            //pMap[i] = ((double)pMapHits[i] / (pMapHits[i] + pMapMisses[i])) * 100;//没有几个边界点
            //计算结果刚开始没有转换为double类型，导致只有很少的边界点
            int temValue = ((double)pMapHits[i] / (pMapHits[i] + pMapMisses[i])) * 100;
            if(temValue > 20) //没有通过阈值分类也只有很少边界点
            {
                pMap[i] = 100;
            }else{
                pMap[i] = temValue;
            }
        }
    }

#endif

#if defined TSDF
    //根据TSDF计算pMap
    for(int i = 0; i < mapParams.height - 2; i++)  //设i为y轴
    {
        for(int j = 0; j < mapParams.width - 2; j++)  //设j为x轴  
        { 
            int gridId_x0_y0 = j + i * mapParams.height;
            int gridId_x1_y0 = gridId_x0_y0 + 1;//往x轴移动一格，移一格
            int gridId_x0_y1 = gridId_x0_y0 + mapParams.height;//往y轴方向移动一格，移一行

            //转世界坐标
            double x_0 = (j - mapParams.offset_x) * mapParams.resolution + mapParams.origin_x;
            double y_0 = (i - mapParams.offset_x) * mapParams.resolution + mapParams.origin_x;
            double x_1 = ((j + 1) - mapParams.offset_x) * mapParams.resolution + mapParams.origin_x;
            double y_1 = ((i + 1) - mapParams.offset_x) * mapParams.resolution + mapParams.origin_x;
            
            double tsdfx0_y0, tsdfx1_y0, tsdfx0_y1, x, y;
            
            tsdfx0_y0 = pMapTSDF[gridId_x0_y0]; tsdfx1_y0 = pMapTSDF[gridId_x1_y0]; tsdfx0_y1 = pMapTSDF[gridId_x0_y1];
            if(tsdfx0_y0 * tsdfx0_y1 < 0)
            { 
                //y方向上两个栅格tsdf值相反，插值应在这两个栅格之间 
                x = x_0;
                if(tsdfx0_y0 != tsdfx0_y1)
                {
                    y = (y_0 * tsdfx0_y1 - y_1 * tsdfx0_y0) / (tsdfx0_y1 - tsdfx0_y0);
                }else{
                    y = y_0;
                }
                pMap[GridIndexToLinearIndex(ConvertWorld2GridIndex(x, y))] = 100;
            }
            else if(tsdfx0_y0 * tsdfx1_y0 < 0)
            { 
                //x方向上两个栅格tsdf值相反，插值应在这两个栅格之间 
                if(tsdfx0_y0 != tsdfx1_y0)
                {
                    x = (x_0 * tsdfx1_y0 - x_1 * tsdfx0_y0) / (tsdfx1_y0 - tsdfx0_y0);
                }else{
                    x = x_0;
                }
                y = y_0;
                pMap[GridIndexToLinearIndex(ConvertWorld2GridIndex(x, y))] = 100;
            }
            // GridIndex id;
            // id.x = j;
            // id.y = i;
            // std::cout << (int)pMap[GridIndexToLinearIndex(id)] << " ";   
        }
    }

#endif

    //end of TODO
    std::cout << "建图完毕" << std::endl;
}

//发布地图．
void PublishMap(ros::Publisher &map_pub)
{
    nav_msgs::OccupancyGrid rosMap;

    rosMap.info.resolution = mapParams.resolution;
    rosMap.info.origin.position.x = 0.0;
    rosMap.info.origin.position.y = 0.0;
    rosMap.info.origin.position.z = 0.0;
    rosMap.info.origin.orientation.x = 0.0;
    rosMap.info.origin.orientation.y = 0.0;
    rosMap.info.origin.orientation.z = 0.0;
    rosMap.info.origin.orientation.w = 1.0;

    rosMap.info.origin.position.x = mapParams.origin_x;
    rosMap.info.origin.position.y = mapParams.origin_y;
    rosMap.info.width = mapParams.width;
    rosMap.info.height = mapParams.height;
    rosMap.data.resize(rosMap.info.width * rosMap.info.height);

    //0~100
    int cnt0, cnt1, cnt2;
    cnt0 = cnt1 = cnt2 = 100;
    for (int i = 0; i < mapParams.width * mapParams.height; i++)
    {
        if (pMap[i] == 50)
        {
            rosMap.data[i] = -1.0;
        }
        else
        {

            rosMap.data[i] = pMap[i];
        }
    }

    rosMap.header.stamp = ros::Time::now();
    rosMap.header.frame_id = "map";

    map_pub.publish(rosMap);
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "OccupanyMapping");

    ros::NodeHandle nodeHandler;

    ros::Publisher mapPub = nodeHandler.advertise<nav_msgs::OccupancyGrid>("laser_map", 1, true);

    std::vector<Eigen::Vector3d> robotPoses;
    std::vector<GeneralLaserScan> generalLaserScans;

    std::string basePath = "/home/beechang/shenlanHW/HW7/OccupanyMappingProject/src/data";

    std::string posePath = basePath + "/pose.txt";
    std::string anglePath = basePath + "/scanAngles.txt";
    std::string scanPath = basePath + "/ranges.txt";

    //读取数据
    ReadPoseInformation(posePath, robotPoses);

    ReadLaserScanInformation(anglePath,
                             scanPath,
                             generalLaserScans);

    //设置地图信息
    SetMapParams();

    OccupanyMapping(generalLaserScans, robotPoses);

    PublishMap(mapPub);

    ros::spin();

    DestoryMap();

    std::cout << "Release Memory!!" << std::endl;
}
