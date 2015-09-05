#include "calcubase.h"
#include "solidPolygon.h"
#include "OpenGLDC.h"
#include "eigen/Eigen"
#include "OBB.h"
using namespace Eigen;
using Eigen::MatrixXd;
class CcutterPath :public CSolidPolygon
{
public:
	CcutterPath();
	~CcutterPath();
public:
	//////////////////初始刀具轨迹
	void update_direct(OBB& theBOX, VECTOR3D &Tdirection);//获得行距大方向
	void get_plane(vector<PVERT>& Bloop,double& A,double& B,double& C);//最小二乘拟合
	void CUT_line(OBB& theBOX,PPOLYGON &firstPLN,vector<vec_VECTOR3D>& CCnomal,vec_PHEDGE &vec_PH,vector<vec_PFACETTRI>& RBOX);//平面切分线
	void cut_partilize(vector<vector<PFACETTRI>>& Rbox, vec_PFACETTRI& m_vecPFacetTri);//重新分区

	inline void get_interPOT(PPOLYGON &firstPLN, PPOLYGON &nextPLN, vec_VECTOR3D &Cnomal, double &h, double &r, vec_PHEDGE &vec_PH, vec_PFACETTRI& facelist, VECTOR3D &Tdirection);//得到下一圈偏置点
	inline void dele_interPOT(PPOLYGON firstPLN, PPOLYGON &nextPLN, vec_PFACETTRI& facelist);//删除折线
	inline void ccpath_snooth(CPlane& cutPlane, vec_PFACETTRI& facelist, PPOLYGON &nextPLN, vec_VECTOR3D &Cnomal, vec_PHEDGE &vec_PH);//轨迹优化
	void change(PFACETTRI& theFac);//改变面片转向
	void get_CCpath(OBB& theBOX, vec_PPOLYPOLYGON& m_vecPPolyPolygons, double h, double r, vector<vec_PFACETTRI>& RBOX);//获得一个叶子分块的轨迹



	//PFACETTRI  OffsetPcc(PFACETTRI pFacet, VECTOR3D StepDir, POINT3D CurrentP, POINT3D& OffsetP, double ScallopHeight, double CutterR);   // 求偏置点
	//POINT3D ProjectionP_Sample(PFACETTRI pFacet, VECTOR3D StepDir, POINT3D CurrentP, CPlane& plane, double ScallopHeight, double CutterR,  double Interval);
	//bool CcutterPath::PlanewithEdge(PHEDGE pHEdge, CPlane Plane, POINT3D& p);
};
