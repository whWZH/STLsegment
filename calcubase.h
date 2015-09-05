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
	CPlane(POINT3D point,VECTOR3D normal);  //�㷨ʽ����
	~CPlane(void);

	POINT3D   m_PointIn;       //ƽ���ϵ�һ����
	VECTOR3D  m_Normal;        //ƽ�淨����
	double    a,b,c,d;         //ƽ�淽�̣�ax+by+cz+d=0

	void GetParameter();       //�ɵ㷨ʽ��÷��̲���
	CPlane BuildPlane(POINT3D a, POINT3D b, POINT3D c);
	CPlane BuildPlane(VECTOR3D n, POINT3D p);  // pΪƽ����һ��
	CPlane BuildPlane(POINT3D p1, POINT3D p2, VECTOR3D NorF);
};
class Ccalcubase 
{
public:
	Ccalcubase(void);
	~Ccalcubase(void);
	//Octree* pOctree; 
	// �����߶���ƽ��Ľ���
	static BOOL IntersectLinePlane(PPOINT3D pStartP, PPOINT3D pEndP, CPlane Plane, POINT3D& InterP);
	static BOOL IntersectLinePlane(PHEDGE pHEdge, CPlane Plane, POINT3D& InterP);
	static BOOL IntersectLinePlane(POINT3D LinePt,VECTOR3D LineNor, PFACETTRI pFacTri, POINT3D& InterP);
	//������Ƭ��ƽ�潻�㣬2013.7.18
	static BOOL IntersectFacPlane(POINT3D PointPre,PHEDGE pHEdge,POINT3D StartP, CPlane Plane,PHEDGE& InterpHE, POINT3D& InterP);
	//��ռ���Ƿ�����Ƭ��,2014.1.4
	static BOOL IsPointInFac(POINT3D Pt,PFACETTRI pFacTri);
	//////////////////���·��
	static void Dijkstra_Point_ST_EN_new(PVERT ST,PVERT EN,map<PVERT,PVERT>& Path,vector<PVERT>& BLOOP,vec_PVERT& m_vecPVert);//����ģ����,ST��EN�����·��
	static BOOL IntersectTriangle(POINT3D orig,VECTOR3D dir,POINT3D v0,POINT3D v1,POINT3D v2,POINT3D &JPOT);//������
	static BOOL IntersectTriangle(POINT3D orig,VECTOR3D dir,POINT3D v0,POINT3D v1,POINT3D v2);
	static POINT3D Center(PHEDGE trigle);//�õ�һ��������Ƭ���е�
	static double dis(POINT3D A,POINT3D B);//�����A����B֮��ľ���
	static double calcuPOINTV(PVERT A);//�����İ�͹��
	static double calcuNoV(PVERT A);//������ƽ����͹��
	static double calcuNoV_new(PVERT A);//����һ�׵�İ�͹��
	static double calcuPOINTV(PFACETTRI A);//������Ƭ�İ�͹��
	static BOOL FindOneRing(PVERT& pVer,vec_PVERT& vecpVer);//����������һ�������
	static BOOL FindOnePH(PVERT A,vec_PHEDGE& vecpH);//����������һ��������
	static BOOL FindTwoRing(PVERT pVer,vec_PFACETTRI& vecpVerT);//���������Ķ���������
	static BOOL FindTwoRing(PVERT pVer,vec_PVERT& vecpVerT);//���������Ķ��������
	static BOOL FindTwoRing_new(PVERT pVer,vec_PVERT& vecpVerT);//���������Ķ��������
	static BOOL FindOneRFac(PVERT pVer,vec_PFACETTRI& vecpFac);////����������һ��������
	static BOOL findThreeRing(PVERT pVer,vec_PVERT& vecpVerR);///��������������������
	static BOOL FindPOneRFAC(PFACETTRI pFac,vec_PFACETTRI& vecpFacP);//��Ƭ��һ��������Ƭ
	static BOOL FindPOneRFAC_NEW(PFACETTRI pFac,vec_PFACETTRI& vecpFacP);//��Ƭ��һ��������Ƭ
	static void FindPtwoRFAC(PFACETTRI pFac, vec_PFACETTRI& vecpFacP);//��Ƭ�Ķ���������Ƭ
	static VECTOR3D CalcuVerNormal(PVERT pVer);//���㶥��ķ�ʸ,�Ե����������������ڼн�ΪȨֵ��������Ƭ��
	static VECTOR3D CalcuPntInFacNormal(POINT3D Pt,PFACETTRI Fac);//����������Ƭ������һ�㷨ʸ
	static VECTOR3D CalcuPntInLineNormal(POINT3D Pt,  PHEDGE pHEdge);//������Ͻ���ķ�ʸ
	static void CalcuVerPrinCurvature(PVERT pVer,double& K1,double& K2);//���㶥���������
	static double CalcuVerMeanCurvature(PVERT pVer);//���㶥���ƽ������
	static double CalcuVerGausCurvature(PVERT pVer);//���㶥��ĸ�˹����
	static double AreaTriMix(PPOINT3D pPoint1, PPOINT3D pPoint2, PPOINT3D pPoint3);//���㶥��һ��������������εĻ�����������������Ϊ��ָ�붥��[pPoint1]���ڵĲ��ֵ������
	static double AreaTriMixSum(PVERT pVer); //һ������������ܺ�
	static BOOL FindRHEd(PHEDGE A,vec_PFACETTRI& m_vecS);
	
	static void create_current_Frenet(MATRIX3D& currentM,VECTOR3D LineNor,POINT3D Center);//�����ֲ�����ϵ,����һ���ֲ�����ϵ
	static void create_BACK_Frenet(MATRIX3D& currentB,VECTOR3D LineNor,POINT3D Center);//�ֲ�����ϵ��ԭ����������
	static void create_current_matriX(MATRIX3D& matriX,double a);//��X��תa��
	static void create_current_matriZ(MATRIX3D& matriZ,double a);//������Z����ת������
	static void create_cone(vector<VECTOR3D>& m_vecVECTOR3D,VECTOR3D LineNor,POINT3D orig,double angA,double angB);//��������׶
	static void create_cone_new(vector<VECTOR3D>& m_vecVECTOR3D,VECTOR3D LineNorZ,VECTOR3D LineNor,POINT3D orig,double angA,int numB);//�������ο���׶
	static void cone_delel(vector<double>& vec_dis);//���ܵ����Ѱ
	static void create_cone_SDF(vector<VECTOR3D>& m_vecVECTOR3D,POINT3D orig,VECTOR3D liner,double angA,double angB);//����SDF׶
	
	////////////////////////////////////////////////////SDF����
	static int SDF_center(vector<double> &vec_center,double xi);//�����ж������ĸ�����
	static double SDF_new_center(vector<double> &vec_c);//ÿһ��������µ�����
	static void SDF_delet(vector<double>& vec_SDF);//��λ��ƽ����֮���
	static void SDF_nomal(vector<double>& vec_SDF,map<double,PFACETTRI>& map_SDF);//SDF��һ��
	static void SDF_nomal_old(vector<double>& vec_SDF,map<double,PFACETTRI>& map_SDF);//SDF��һ��
	static void SDF_nomal(vector<double>& vec_SDF,map<double,PVERT>& map_SDF);//SDF��һ��
	/////////////////////////////////////////////////////GMM����
	static void GMM_seg(vector<double>& vec_SDF,vector<vector<PFACETTRI>>& Rbox,map<double,PFACETTRI>& map_SDF,int k,vector<double>& vec_u,vector<double>& vec_ct,vector<double>& vec_zmj);
	static double GMM_norm(double x,double u,double ct);//��������ܶ�
	static double GMM_E(int zi,vector<double>& vec_u,vector<double>& vec_ct,vector<double>& vec_zmj,vector<double>& vec_SDF,vector<double>& vec_wj);//GMM����E��,����wj
	static void GMM_M(vector<double>& vec_SDF,vector<double>& vec_wj,double& u,double& ct,double& mj);//������������
	static int GMM_FacetoZi(vector<double>& vec_u,vector<double>& vec_ct,vector<double>& vec_zmj,double& theface);
	///////////////////////////////////////////���ƾ��������
	static void segment_cloud(vec_PVERT& sdf_point,vector<vector<PFACETTRI>>& Rbox,vector<PFACETTRI>& m_vecPFacetTri);//���������ֿ�
	static int Center_kmeans_cloud(vector<int>& cloud,vec_PVERT& sdf_point);//��һ����������
	static int cloud_center(vector<int>& vec_center,PVERT& xin,vec_PVERT& sdf_point);//�ж������ĸ�����
	// �������������
	static double AreaTri(PPOINT3D pPoint1, PPOINT3D pPoint2, PPOINT3D pPoint3);//�������������
	static double AreaTri(POINT3D Point1,POINT3D Point2,POINT3D Point3);
	static double AreaTri(PFACETTRI pFacTri);
	//������
	static double  CalPointOnLineCurv(VECTOR3D &Tang, PHEDGE &pHEdge, POINT3D &p);   // �����ⷽ��ı��ϵ�����
	static double CalculatePoint_Curv(PFACETTRI pFacet, POINT3D p, VECTOR3D vec);//�������ڲ����ϵ�����
	static void  CalculateCurvByTNT(PVERT pVex,VECTOR3D &t1,VECTOR3D &t2,double &K1,double &K2);//�����������
	 static void  VexTwoNeighVex(PVERT pVex, vec_PVERT& vecpVex);
	 static void  VexNeighVex(PVERT pVex, vec_PVERT& vecpVex);
	 static void  CalculateCurvature(double a, double b, double c, double p, double q, VECTOR3D& T1, VECTOR3D& T2,double &K1,double &K2);//��ɢ����
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
	 static void  CalculateVexNor(PVERT  pVex, VECTOR3D& VexNormal);   // ���㷨ʸ
	 static void  FacetNeighVex(PFACETTRI pFacet, vec_PVERT& vecpVex);
	 static double SameSide(POINT3D A, POINT3D B, POINT3D C, POINT3D P);
	 static void  VexNeighFacet(PVERT pVex, vec_PFACETTRI& vecpFacet);
	 static VECTOR3D CalPointInTriNorByCore(PFACETTRI pFacet, POINT3D p);   // ���Ȩֵ��
	 static void  CalOnLineNorByFacNor(PHEDGE pHEdge, VECTOR3D& pNormal);  // ��Ƭ��ʸƽ��
	 static bool IsHeight(POINT3D p0, POINT3D p1, POINT3D p2);//�ж��Ƿ�߲�
protected:
private:
};