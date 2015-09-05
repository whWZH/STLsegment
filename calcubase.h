#pragma once
#include "solidMesh.h"
#include "octree.h"
#include "KD-tree/KDtree.h"
#include "KD-tree/Vec.h"
#include <queue>
#include <cmath>
#include "TNT/jama_lu.h" 
#include<memory>
using namespace JAMA;

class CPlane
{
public:
	CPlane(void);
	CPlane(POINT3D point,VECTOR3D normal);  //点法式方程
	~CPlane(void);

	POINT3D   m_PointIn;       //平面上的一个点
	VECTOR3D  m_Normal;        //平面法向量
	double    a,b,c,d;         //平面方程：ax+by+cz+d=0

	void GetParameter();       //由点法式获得方程参数
	CPlane BuildPlane(POINT3D a, POINT3D b, POINT3D c);
	CPlane BuildPlane(VECTOR3D n, POINT3D p);  // p为平面上一点
	CPlane BuildPlane(POINT3D p1, POINT3D p2, VECTOR3D NorF);
};
class Ccalcubase 
{
public:
	Ccalcubase(void);
	~Ccalcubase(void);
	//Octree* pOctree; 
	// 计算线段与平面的交点
	static BOOL IntersectLinePlane(PPOINT3D pStartP, PPOINT3D pEndP, CPlane Plane, POINT3D& InterP);
	static BOOL IntersectLinePlane(PHEDGE pHEdge, CPlane Plane, POINT3D& InterP);
	static BOOL IntersectLinePlane(POINT3D LinePt,VECTOR3D LineNor, PFACETTRI pFacTri, POINT3D& InterP);
	//三角面片与平面交点，2013.7.18
	static BOOL IntersectFacPlane(POINT3D PointPre,PHEDGE pHEdge,POINT3D StartP, CPlane Plane,PHEDGE& InterpHE, POINT3D& InterP);
	//求空间点是否在面片上,2014.1.4
	static BOOL IsPointInFac(POINT3D Pt,PFACETTRI pFacTri);
	//////////////////最短路径
	static void Dijkstra_Point_ST_EN_new(PVERT ST,PVERT EN,map<PVERT,PVERT>& Path,vector<PVERT>& BLOOP,vec_PVERT& m_vecPVert);//整个模型上,ST到EN的最短路径
	static BOOL IntersectTriangle(POINT3D orig,VECTOR3D dir,POINT3D v0,POINT3D v1,POINT3D v2,POINT3D &JPOT);//射线求交
	static BOOL IntersectTriangle(POINT3D orig,VECTOR3D dir,POINT3D v0,POINT3D v1,POINT3D v2);
	static POINT3D Center(PHEDGE trigle);//得到一个三角面片的中点
	static double dis(POINT3D A,POINT3D B);//计算点A到点B之间的距离
	static double calcuPOINTV(PVERT A);//计算点的凹凸度
	static double calcuNoV(PVERT A);//计算点的平均凹凸度
	static double calcuNoV_new(PVERT A);//计算一阶点的凹凸度
	static double calcuPOINTV(PFACETTRI A);//计算面片的凹凸度
	static BOOL FindOneRing(PVERT& pVer,vec_PVERT& vecpVer);//计算网格点的一阶领域点
	static BOOL FindOnePH(PVERT A,vec_PHEDGE& vecpH);//计算网格点的一阶领域半边
	static BOOL FindTwoRing(PVERT pVer,vec_PFACETTRI& vecpVerT);//计算网格点的二阶领域面
	static BOOL FindTwoRing(PVERT pVer,vec_PVERT& vecpVerT);//计算网格点的二阶领域点
	static BOOL FindTwoRing_new(PVERT pVer,vec_PVERT& vecpVerT);//计算网格点的二阶领域点
	static BOOL FindOneRFac(PVERT pVer,vec_PFACETTRI& vecpFac);////计算网格点的一阶领域面
	static BOOL findThreeRing(PVERT pVer,vec_PVERT& vecpVerR);///计算网格点的三阶领域面
	static BOOL FindPOneRFAC(PFACETTRI pFac,vec_PFACETTRI& vecpFacP);//面片的一阶邻域面片
	static BOOL FindPOneRFAC_NEW(PFACETTRI pFac,vec_PFACETTRI& vecpFacP);//面片的一阶邻域面片
	static void FindPtwoRFAC(PFACETTRI pFac, vec_PFACETTRI& vecpFacP);//面片的二阶领域面片
	static VECTOR3D CalcuVerNormal(PVERT pVer);//计算顶点的法矢,以点相邻面的面积和相邻夹角为权值（三角面片）
	static VECTOR3D CalcuPntInFacNormal(POINT3D Pt,PFACETTRI Fac);//计算三角面片上任意一点法矢
	static VECTOR3D CalcuPntInLineNormal(POINT3D Pt,  PHEDGE pHEdge);//计算边上交点的法矢
	static void CalcuVerPrinCurvature(PVERT pVer,double& K1,double& K2);//计算顶点的主曲率
	static double CalcuVerMeanCurvature(PVERT pVer);//计算顶点的平均曲率
	static double CalcuVerGausCurvature(PVERT pVer);//计算顶点的高斯曲率
	static double AreaTriMix(PPOINT3D pPoint1, PPOINT3D pPoint2, PPOINT3D pPoint3);//计算顶点一阶领域锐角三角形的混合面积（三角以外心为界分割，与顶点[pPoint1]相邻的部分的面积）
	static double AreaTriMixSum(PVERT pVer); //一阶领域混合面积总和
	static BOOL FindRHEd(PHEDGE A,vec_PFACETTRI& m_vecS);
	
	static void create_current_Frenet(MATRIX3D& currentM,VECTOR3D LineNor,POINT3D Center);//建立局部坐标系,返回一个局部坐标系
	static void create_BACK_Frenet(MATRIX3D& currentB,VECTOR3D LineNor,POINT3D Center);//局部坐标系还原成世界坐标
	static void create_current_matriX(MATRIX3D& matriX,double a);//绕X旋转a度
	static void create_current_matriZ(MATRIX3D& matriZ,double a);//建立绕Z轴旋转的坐标
	static void create_cone(vector<VECTOR3D>& m_vecVECTOR3D,VECTOR3D LineNor,POINT3D orig,double angA,double angB);//建立可视锥
	static void create_cone_new(vector<VECTOR3D>& m_vecVECTOR3D,VECTOR3D LineNorZ,VECTOR3D LineNor,POINT3D orig,double angA,int numB);//建立扇形可视锥
	static void cone_delel(vector<double>& vec_dis);//智能点的搜寻
	static void create_cone_SDF(vector<VECTOR3D>& m_vecVECTOR3D,POINT3D orig,VECTOR3D liner,double angA,double angB);//建立SDF锥
	
	////////////////////////////////////////////////////SDF聚类
	static int SDF_center(vector<double> &vec_center,double xi);//初步判断属于哪个质心
	static double SDF_new_center(vector<double> &vec_c);//每一类里求出新的质心
	static void SDF_delet(vector<double>& vec_SDF);//中位数平均差之间的
	static void SDF_nomal(vector<double>& vec_SDF,map<double,PFACETTRI>& map_SDF);//SDF归一化
	static void SDF_nomal_old(vector<double>& vec_SDF,map<double,PFACETTRI>& map_SDF);//SDF归一化
	static void SDF_nomal(vector<double>& vec_SDF,map<double,PVERT>& map_SDF);//SDF归一化
	/////////////////////////////////////////////////////GMM聚类
	static void GMM_seg(vector<double>& vec_SDF,vector<vector<PFACETTRI>>& Rbox,map<double,PFACETTRI>& map_SDF,int k,vector<double>& vec_u,vector<double>& vec_ct,vector<double>& vec_zmj);
	static double GMM_norm(double x,double u,double ct);//计算概率密度
	static double GMM_E(int zi,vector<double>& vec_u,vector<double>& vec_ct,vector<double>& vec_zmj,vector<double>& vec_SDF,vector<double>& vec_wj);//GMM聚类E步,返回wj
	static void GMM_M(vector<double>& vec_SDF,vector<double>& vec_wj,double& u,double& ct,double& mj);//跟新三个参数
	static int GMM_FacetoZi(vector<double>& vec_u,vector<double>& vec_ct,vector<double>& vec_zmj,double& theface);
	///////////////////////////////////////////点云聚类的两个
	static void segment_cloud(vec_PVERT& sdf_point,vector<vector<PFACETTRI>>& Rbox,vector<PFACETTRI>& m_vecPFacetTri);//种子生长分块
	static int Center_kmeans_cloud(vector<int>& cloud,vec_PVERT& sdf_point);//从一类里找中心
	static int cloud_center(vector<int>& vec_center,PVERT& xin,vec_PVERT& sdf_point);//判断属于哪个质心
	// 计算三角形面积
	static double AreaTri(PPOINT3D pPoint1, PPOINT3D pPoint2, PPOINT3D pPoint3);//计算三角形面积
	static double AreaTri(POINT3D Point1,POINT3D Point2,POINT3D Point3);
	static double AreaTri(PFACETTRI pFacTri);
	//点曲率
	static double  CalPointOnLineCurv(VECTOR3D &Tang, PHEDGE &pHEdge, POINT3D &p);   // 沿任意方向的边上点曲率
	static double CalculatePoint_Curv(PFACETTRI pFacet, POINT3D p, VECTOR3D vec);//三角形内部点上的曲率
	static void  CalculateCurvByTNT(PVERT pVex,VECTOR3D &t1,VECTOR3D &t2,double &K1,double &K2);//点的连续曲率
	 static void  VexTwoNeighVex(PVERT pVex, vec_PVERT& vecpVex);
	 static void  VexNeighVex(PVERT pVex, vec_PVERT& vecpVex);
	 static void  CalculateCurvature(double a, double b, double c, double p, double q, VECTOR3D& T1, VECTOR3D& T2,double &K1,double &K2);//离散曲率
	 static double  FacetTriArea(PFACETTRI pFacet, double& InnerAngle);
	 static double  TriArea(POINT3D p1, POINT3D p2, POINT3D p3);

	 /////////yT
	 static void CalculatePointNor(PFACETTRI pFacet, POINT3D p, VECTOR3D& pNormal);
	 static void  Facet_nPNeighFacet(PFACETTRI pFacet, vec_PFACETTRI& vecpFacet);
	 static bool PointInTri(PFACETTRI pFacet, POINT3D P);
	 static void CalculateNor(PFACETTRI pFacet, POINT3D p, VECTOR3D& pNormal);
	 static bool IntersectLinePlaneYT(PFACETTRI pFacTri, POINT3D LinePt, VECTOR3D LineNor, POINT3D& InterP);
	 static bool Line_TriIntersect(PFACETTRI pFacet, POINT3D LineP, VECTOR3D LineDir, POINT3D& IntersectP);
	 static bool PointOnLine(PFACETTRI pFacet, POINT3D P, PHEDGE& pHEdge);
	 static void  CalculateVexNor(PVERT  pVex, VECTOR3D& VexNormal);   // 顶点法矢
	 static void  FacetNeighVex(PFACETTRI pFacet, vec_PVERT& vecpVex);
	 static double SameSide(POINT3D A, POINT3D B, POINT3D C, POINT3D P);
	 static void  VexNeighFacet(PVERT pVex, vec_PFACETTRI& vecpFacet);
	 static VECTOR3D CalPointInTriNorByCore(PFACETTRI pFacet, POINT3D p);   // 面积权值法
	 static void  CalOnLineNorByFacNor(PHEDGE pHEdge, VECTOR3D& pNormal);  // 面片法矢平均
	 static bool IsHeight(POINT3D p0, POINT3D p1, POINT3D p2);//判断是否高波
protected:
private:
};